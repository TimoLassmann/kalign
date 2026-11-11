#include "tldevel.h"
#include "tlrng.h"

#include "msa_struct.h"

#include "threadpool.h"

#include "task.h"

#include "aln_param.h"

#include "aln_struct.h"
#include "aln_mem.h"
#include "aln_setup.h"
#include "aln_controller.h"

#include "weave_alignment.h"
/* #include "weave_alignment.h" */
#include "sp_score.h"
#include "anchor_consistency.h"
#include "msa_consistency.h"
#include "aln_vsm.h"

#include <string.h>
#include <stdatomic.h>
#include <float.h>

#define ALN_RUN_IMPORT
#include "aln_run.h"

static int do_align(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id);
static int do_align_inline_refine(struct msa* msa, struct aln_tasks* t, struct aln_mem* m, int task_id, int n_trials);

/* ---------------------------------------------------------------------------
 * Non-recursive guide-tree walk.
 *
 * Every internal node (task) keeps a count of unfinished internal children.
 * Nodes whose children are all leaves are ready at the start.  A thread that
 * finishes a node decrements its parent's count and, if it reaches zero,
 * continues with the parent in the same loop.  No task ever waits for another
 * tree task, so the C stack use of a thread does not depend on tree depth.
 * Each node's alignment depends only on its two children, so results do not
 * depend on the execution order.
 * ------------------------------------------------------------------------- */

typedef atomic_int walk_cnt_t;

struct walk_ctx {
        struct msa *msa;
        struct aln_tasks *t;
        struct aln_param *ap;
        uint8_t *active;
        int n_trials;           /* < 0: plain alignment, >= 0: inline refinement */
        int *parent;            /* parent task index, -1 for the root */
        walk_cnt_t *pending;    /* unfinished internal children */
};

struct walk_arg {
        struct walk_ctx *ctx;
        int node;
};

/* Returns the new count. */
static int walk_dec_pending(struct walk_ctx *w, int p)
{
        int old;
        old = atomic_fetch_sub_explicit(&w->pending[p], 1, memory_order_acq_rel);
        return old - 1;
}

/* Align node c, then keep going up the tree while the parent becomes ready. */
static void walk_run_chain(struct walk_ctx *w, int c)
{
        struct msa *msa = w->msa;
        struct aln_tasks *t = w->t;
        struct aln_param *ap = w->ap;

        while(c >= 0){
                struct task *local_t = t->list[c];
                struct aln_mem *ml = NULL;

                alloc_aln_mem(&ml, 256);
                ml->ap = ap;
                ml->mode = ALN_MODE_FULL;
                ml->run_parallel = msa->run_parallel;
                ml->pool = ap->pool;
                if(w->n_trials < 0){
                        do_align(msa, t, ml, c);
                }else{
                        do_align_inline_refine(msa, t, ml, c, w->n_trials);
                }

                w->active[local_t->a] = 0;
                w->active[local_t->b] = 0;
                w->active[local_t->c] = 1;
                free_aln_mem(ml);

                int p = w->parent[c];
                if(p >= 0 && walk_dec_pending(w, p) == 0){
                        c = p;
                }else{
                        c = -1;
                }
        }
}

static void walk_task(void *arg)
{
        struct walk_arg *a = (struct walk_arg *)arg;
        walk_run_chain(a->ctx, a->node);
}

static int walk_tree(struct msa* msa, struct aln_tasks* t, struct aln_param* ap,
                     uint8_t* active, int n_trials)
{
        struct walk_ctx w;
        int n = t->n_tasks;
        int i;
        int *ready = NULL;
        int n_ready = 0;
        struct walk_arg *args = NULL;

        w.msa = msa;
        w.t = t;
        w.ap = ap;
        w.active = active;
        w.n_trials = n_trials;
        w.parent = NULL;
        w.pending = NULL;

        MMALLOC(w.parent, sizeof(int) * n);
        MMALLOC(w.pending, sizeof(walk_cnt_t) * n);
        MMALLOC(ready, sizeof(int) * n);
        for(i = 0; i < n; i++){
                w.parent[i] = -1;
                w.pending[i] = 0;
        }
        for(i = 0; i < n; i++){
                int ch[2];
                int j;
                ch[0] = t->list[i]->a - msa->numseq;
                ch[1] = t->list[i]->b - msa->numseq;
                for(j = 0; j < 2; j++){
                        if(ch[j] >= 0){
                                w.parent[ch[j]] = i;
                                w.pending[i]++;
                        }
                }
        }
        for(i = 0; i < n; i++){
                if(w.pending[i] == 0){
                        ready[n_ready++] = i;
                }
        }

        if(ap->pool){
                tp_group_t *g = tp_group_create(ap->pool);
                MMALLOC(args, sizeof(struct walk_arg) * n_ready);
                for(i = 0; i < n_ready; i++){
                        args[i].ctx = &w;
                        args[i].node = ready[i];
                        if(!g || tp_group_submit(g, walk_task, &args[i]) != 0){
                                walk_run_chain(&w, ready[i]);
                        }
                }
                if(g){
                        tp_group_wait(g);
                        tp_group_destroy(g);
                }
                MFREE(args);
        }else{
                for(i = 0; i < n_ready; i++){
                        walk_run_chain(&w, ready[i]);
                }
        }

        MFREE(ready);
        MFREE(w.pending);
        MFREE(w.parent);
        return OK;
ERROR:
        if(args){
                MFREE(args);
        }
        if(ready){
                MFREE(ready);
        }
        if(w.pending){
                MFREE(w.pending);
        }
        if(w.parent){
                MFREE(w.parent);
        }
        return FAIL;
}

int create_msa_tree(struct msa* msa, struct aln_param* ap,struct aln_tasks* t)
{ KPERF_PHASE_REGION(KP_PROGRESSIVE, KPH_PROGRESSIVE);
        int i;
        uint8_t* active = NULL;

        RUN(sort_tasks(t, TASK_ORDER_TREE));

        MMALLOC(active, sizeof(uint8_t)* msa->num_profiles);

        for(i = 0; i < msa->numseq;i++){
                active[i] = 1;
        }
        for(i = msa->numseq; i < msa->num_profiles;i++){
                active[i] = 0;
        }
        /* LOG_MSG("Setting threads to 1 for debugging!"); */
        /* ap->nthreads = 1; */
        msa->run_parallel = 1;
        if(ap->nthreads == 1){
                msa->run_parallel = 0;
        }

        RUN(walk_tree(msa, t, ap, active, -1));

        MFREE(active);
        return OK;
ERROR:
        if(active){
                MFREE(active);
        }
        return FAIL;
}


float compute_gap_scale(struct msa* msa, struct aln_param* ap, int a, int b)
{
        float ds = ap->dist_scale;
        if(ds <= 0.0f || msa->seq_distances == NULL){
                return 1.0f;
        }

        /* Compute mean normalized distance for all sequences in both clusters */
        float sum = 0.0f;
        int count = 0;
        int i;
        for(i = 0; i < msa->nsip[a]; i++){
                int si = msa->sip[a][i];
                if(si < msa->numseq){
                        sum += msa->seq_distances[si];
                        count++;
                }
        }
        for(i = 0; i < msa->nsip[b]; i++){
                int si = msa->sip[b][i];
                if(si < msa->numseq){
                        sum += msa->seq_distances[si];
                        count++;
                }
        }
        if(count == 0){
                return 1.0f;
        }
        float avg_div = sum / (float)count;

        /* Scale: 1.0 for similar sequences, decreasing for divergent.
           scale = max(0.3, 1.0 - dist_scale * avg_div)
           With dist_scale=0.5 and avg_div=1.0 (fully divergent): scale=0.5
           With dist_scale=0.5 and avg_div=0.0 (identical): scale=1.0 */
        float scale = 1.0f - ds * avg_div;
        if(scale < 0.3f) scale = 0.3f;
        if(scale > 1.0f) scale = 1.0f;
        return scale;
}

float compute_subm_offset(struct msa* msa, struct aln_param* ap, int a, int b)
{
        float amax = ap->vsm_amax;
        /* Legacy VSM only; KALIGN_VSM_MERGE is handled by vsm_merge_offset() */
        if(ap->vsm_mode != KALIGN_VSM_LEGACY || amax <= 0.0f || msa->seq_distances == NULL){
                return 0.0f;
        }

        /* Compute mean normalized distance for all sequences in both clusters */
        float sum = 0.0f;
        int count = 0;
        int i;
        for(i = 0; i < msa->nsip[a]; i++){
                int si = msa->sip[a][i];
                if(si < msa->numseq){
                        sum += msa->seq_distances[si];
                        count++;
                }
        }
        for(i = 0; i < msa->nsip[b]; i++){
                int si = msa->sip[b][i];
                if(si < msa->numseq){
                        sum += msa->seq_distances[si];
                        count++;
                }
        }
        if(count == 0){
                return 0.0f;
        }
        float avg_div = sum / (float)count;

        /* MAFFT VSM: a(d) = max(0, amax - d)
           Close sequences (small d) -> large offset -> more stringent scoring
           Distant sequences (large d) -> small/zero offset -> original scoring */
        float offset = amax - avg_div;
        if(offset < 0.0f) offset = 0.0f;

        return offset;
}

/* Leaf weight for make_profile_n — always 1.0 since sequence weighting
   is handled via balanced freq-count merging in update_n. */
static float leaf_weight(struct msa* msa, int node)
{
        (void)msa; (void)node;
        return 1.0f;
}

int do_align(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id)
{
        float* tmp = NULL;
        int a,b,c;
        int len_a;
        int len_b;
        int j,g;

        a = t->list[task_id]->a;
        b = t->list[task_id]->b;
        c = t->list[task_id]->c;

        /* Distance-dependent parameter scaling */
        struct aln_param* orig_ap = m->ap;
        struct aln_param scaled_ap;
        float* vsm_save_a = NULL;
        float* vsm_save_b = NULL;
        float gap_scale = compute_gap_scale(msa, m->ap, a, b);
        float subm_off = compute_subm_offset(msa, m->ap, a, b);
        if(gap_scale < 1.0f || subm_off > 0.0f){
                scaled_ap = *m->ap;  /* shallow copy — shares subm pointer */
                scaled_ap.gpo *= gap_scale;
                scaled_ap.gpe *= gap_scale;
                scaled_ap.tgpe *= gap_scale;
                scaled_ap.subm_offset = subm_off;
                m->ap = &scaled_ap;
        }

        if(msa->nsip[a] == 1){
                m->len_a = msa->sequences[a]->len;//  aln->sl[a];
                RUN(make_profile_n(m->ap, msa->sequences[a]->s,m->len_a,leaf_weight(msa,a),&t->profile[a]));
        }else{
                m->len_a = msa->plen[a];
                RUN(set_gap_penalties_n(t->profile[a],m->len_a,msa->nsip[b]));
        }

        if(msa->nsip[b] == 1){
                m->len_b = msa->sequences[b]->len;// aln->sl[b];
                RUN(make_profile_n(m->ap, msa->sequences[b]->s,m->len_b,leaf_weight(msa,b),&t->profile[b]));
        }else{
                m->len_b = msa->plen[b];
                RUN(set_gap_penalties_n(t->profile[b],m->len_b,msa->nsip[a]));
        }

        {
                float vsm_off = 0.0f;
                RUN(vsm_merge_offset(msa, t, m->ap, a, b, m->len_a, m->len_b, &vsm_off));
                if(vsm_off > 0.0f){
                        if(m->ap != &scaled_ap){
                                scaled_ap = *m->ap;
                                m->ap = &scaled_ap;
                        }
                        scaled_ap.subm_offset = vsm_off;
                        RUN(vsm_shift_profiles(t, a, b, m->len_a, m->len_b, vsm_off, &vsm_save_a, &vsm_save_b));
                }
        }
        RUN(init_alnmem(m));

        m->margin_sum = 0.0F;
        m->margin_count = 0;
        m->consistency = NULL;

        /* Compute consistency bonus for all merge types */
        {
                int dp_row_node, dp_col_node, dp_rows, dp_cols;
                int have_consistency = (msa->poar_consistency != NULL)
                        || (msa->consistency_table != NULL);

                if(have_consistency){
                        if(msa->nsip[a] == 1 && msa->nsip[b] == 1){
                                if(m->len_a < m->len_b){
                                        dp_row_node = a; dp_rows = m->len_a;
                                        dp_col_node = b; dp_cols = m->len_b;
                                }else{
                                        dp_row_node = b; dp_rows = m->len_b;
                                        dp_col_node = a; dp_cols = m->len_a;
                                }
                        }else if(msa->nsip[a] == 1){
                                dp_row_node = b; dp_rows = m->len_b;
                                dp_col_node = a; dp_cols = m->len_a;
                        }else if(msa->nsip[b] == 1){
                                dp_row_node = a; dp_rows = m->len_a;
                                dp_col_node = b; dp_cols = m->len_b;
                        }else{
                                if(m->len_a < m->len_b){
                                        dp_row_node = a; dp_rows = m->len_a;
                                        dp_col_node = b; dp_cols = m->len_b;
                                }else{
                                        dp_row_node = b; dp_rows = m->len_b;
                                        dp_col_node = a; dp_cols = m->len_a;
                                }
                        }

                        if(msa->poar_consistency != NULL){
                                RUN(poar_consistency_get_bonus(
                                        (struct poar_consistency_ctx*)msa->poar_consistency,
                                        msa, dp_row_node, dp_rows, dp_col_node, dp_cols,
                                        &m->consistency));
                        }else{
                                RUN(anchor_consistency_get_bonus_profile(
                                        (struct consistency_table*)msa->consistency_table,
                                        msa, dp_row_node, dp_rows, dp_col_node, dp_cols,
                                        &m->consistency));
                        }
                }
        }

        m->mode = ALN_MODE_FULL;
        if(msa->nsip[a] == 1){
                if(msa->nsip[b] == 1){
                        if(m->len_a < m->len_b){
                                m->seq1 = msa->sequences[a]->s;
                                m->seq2 = msa->sequences[b]->s;
                                /* LOG_MSG("%d %d", m->len_a, m->len_b); */
                                m->prof1 = NULL;
                                m->prof2 = NULL;
                                aln_runner(m);
                        }else{
                                len_b = m->len_b;
                                len_a = m->len_a;

                                m->enda = len_b;
                                m->endb = len_a;
                                m->len_a = len_b;
                                m->len_b = len_a;

                                m->seq1 = msa->sequences[b]->s;
                                m->seq2 = msa->sequences[a]->s;
                                m->prof1 = NULL;
                                m->prof2 = NULL;

                                aln_runner(m);

                                RUN(mirror_path_n(m,len_a,len_b));
                                m->len_a = len_a;
                                m->len_b = len_b;
                        }
                        /* m->seq1 = msa->sequences[a]->s; */
                        /* m->seq2 = msa->sequences[b]->s; */
                        /* m->prof1 = NULL; */
                        /* m->prof2 = NULL; */
                        /* aln_runner(m); */
                }else{
                        len_b = m->len_b;
                        len_a = m->len_a;

                        m->enda = len_b;
                        m->endb = len_a;
                        m->len_a = len_b;
                        m->len_b = len_a;

                        m->seq1 = NULL;
                        m->seq2 = msa->sequences[a]->s;
                        m->prof1 = t->profile[b];
                        m->prof2 = NULL;
                        m->sip = msa->nsip[b];

                        aln_runner(m);
                        RUN(mirror_path_n(m, len_a,len_b));
                        m->len_a = len_a;
                        m->len_b = len_b;
                }
        }else{
                if(msa->nsip[b] == 1){
                        m->seq1 = NULL;
                        m->seq2 = msa->sequences[b]->s;
                        m->prof1 = t->profile[a];
                        m->prof2 = NULL;
                        m->sip = msa->nsip[a];
                        aln_runner(m);
                }else{
                        if(m->len_a < m->len_b){
                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = t->profile[a];
                                m->prof2 = t->profile[b];
                                aln_runner(m);
                        }else{
                                len_b = m->len_b;
                                len_a = m->len_a;

                                m->enda = len_b;
                                m->endb = len_a;
                                m->len_a = len_b;
                                m->len_b = len_a;

                                m->seq1 = NULL;
                                m->seq2 = NULL;
                                m->prof1 = t->profile[b];
                                m->prof2 = t->profile[a];

                                aln_runner(m);

                                RUN(mirror_path_n(m,len_a,len_b));
                                m->len_a = len_a;
                                m->len_b = len_b;
                        }
                }
        }

        /* Store alignment confidence (average meetup margin) */
        if(m->margin_count > 0){
                t->list[task_id]->confidence = m->margin_sum / (float)m->margin_count;
        }else{
                t->list[task_id]->confidence = 0.0F;
        }

        RUN(add_gap_info_to_path_n(m)) ;
        /* LOG_MSG("Aligned %d and %d (len %d %d) -> path is of length: %d",a,b, m->len_a,m->len_b, 64*(m->path[0]+2)); */

        /* Free consistency bonus if allocated */
        if(m->consistency){
                sparse_bonus_free(m->consistency);
                m->consistency = NULL;
        }

        vsm_unshift_profiles(t, a, b, &vsm_save_a, &vsm_save_b);

        /* Restore original aln_param for profile update (unscaled base penalties) */
        m->ap = orig_ap;

        MMALLOC(tmp,sizeof(float)*64*(m->path[0]+2));

        if(task_id != t->n_tasks-1){
                update_n(t->profile[a],t->profile[b],tmp,m->ap,m->path,msa->nsip[a],msa->nsip[b]);
        }

        MFREE(t->profile[a]);
        MFREE(t->profile[b]);

        t->profile[c] = tmp;
        RUN(make_seq(msa,a,b,m->path));

        msa->plen[c] = m->path[0];

        msa->nsip[c] = msa->nsip[a] + msa->nsip[b];

        MREALLOC(msa->sip[c],sizeof(int)*(msa->nsip[a] + msa->nsip[b]));

        g = 0;
        for (j = msa->nsip[a];j--;){
                msa->sip[c][g] = msa->sip[a][j];
                g++;
        }
        for (j = msa->nsip[b];j--;){
                msa->sip[c][g] = msa->sip[b][j];
                g++;
        }

        return OK;
ERROR:
        /* Leave the profiles and parameters as they were before the merge. */
        if(m->consistency){
                sparse_bonus_free(m->consistency);
                m->consistency = NULL;
        }
        vsm_unshift_profiles(t, a, b, &vsm_save_a, &vsm_save_b);
        m->ap = orig_ap;
        return FAIL;
}

/* ---------------------------------------------------------------------------
 * Inline refinement: run N trials at each progressive merge, keep the best
 * by SP score. Single pass — no second replay pass needed.
 * ------------------------------------------------------------------------- */

int create_msa_tree_inline_refine(struct msa* msa, struct aln_param* ap,
                                  struct aln_tasks* t, int n_trials)
{ KPERF_PHASE_REGION(KP_PROGRESSIVE, KPH_PROGRESSIVE);
        int i;
        uint8_t* active = NULL;

        RUN(sort_tasks(t, TASK_ORDER_TREE));

        MMALLOC(active, sizeof(uint8_t) * msa->num_profiles);

        for(i = 0; i < msa->numseq; i++){
                active[i] = 1;
        }
        for(i = msa->numseq; i < msa->num_profiles; i++){
                active[i] = 0;
        }

        /* Tree-level fork-join is safe (sibling edges are independent).
           Multi-trial loop within each edge stays serial.
           Hirschberg fwd/bwd parallelism enabled via run_parallel. */
        msa->run_parallel = 1;
        if(ap->nthreads == 1){
                msa->run_parallel = 0;
        }

        RUN(walk_tree(msa, t, ap, active, n_trials));

        MFREE(active);
        return OK;
ERROR:
        if(active){
                MFREE(active);
        }
        return FAIL;
}

/* Run n_trials alignments for one edge, keep the best by SP score.
   Trial 0 is the deterministic baseline. Trials 1..n-1 use flip_threshold
   to explore alternative Hirschberg paths at uncertain meetup points. */
int do_align_inline_refine(struct msa* msa, struct aln_tasks* t,
                           struct aln_mem* m, int task_id, int n_trials)
{
        float* tmp = NULL;
        int* best_path = NULL;
        int a, b, c;
        int len_a, len_b;
        int j, g, k;
        float best_sp = -FLT_MAX;
        float avg_margin = 0.0F;

        a = t->list[task_id]->a;
        b = t->list[task_id]->b;
        c = t->list[task_id]->c;

        /* Distance-dependent parameter scaling */
        struct aln_param* orig_ap = m->ap;
        struct aln_param scaled_ap;
        float* vsm_save_a = NULL;
        float* vsm_save_b = NULL;
        float gap_scale = compute_gap_scale(msa, m->ap, a, b);
        float subm_off = compute_subm_offset(msa, m->ap, a, b);
        if(gap_scale < 1.0f || subm_off > 0.0f){
                scaled_ap = *m->ap;
                scaled_ap.gpo *= gap_scale;
                scaled_ap.gpe *= gap_scale;
                scaled_ap.tgpe *= gap_scale;
                scaled_ap.subm_offset = subm_off;
                m->ap = &scaled_ap;
        }

        /* Build profiles (same as do_align) */
        if(msa->nsip[a] == 1){
                m->len_a = msa->sequences[a]->len;
                RUN(make_profile_n(m->ap, msa->sequences[a]->s, m->len_a, leaf_weight(msa,a), &t->profile[a]));
        }else{
                m->len_a = msa->plen[a];
                RUN(set_gap_penalties_n(t->profile[a], m->len_a, msa->nsip[b]));
        }

        if(msa->nsip[b] == 1){
                m->len_b = msa->sequences[b]->len;
                RUN(make_profile_n(m->ap, msa->sequences[b]->s, m->len_b, leaf_weight(msa,b), &t->profile[b]));
        }else{
                m->len_b = msa->plen[b];
                RUN(set_gap_penalties_n(t->profile[b], m->len_b, msa->nsip[a]));
        }

        len_a = m->len_a;
        len_b = m->len_b;

        {
                float vsm_off = 0.0f;
                RUN(vsm_merge_offset(msa, t, m->ap, a, b, len_a, len_b, &vsm_off));
                if(vsm_off > 0.0f){
                        if(m->ap != &scaled_ap){
                                scaled_ap = *m->ap;
                                m->ap = &scaled_ap;
                        }
                        scaled_ap.subm_offset = vsm_off;
                        RUN(vsm_shift_profiles(t, a, b, len_a, len_b, vsm_off, &vsm_save_a, &vsm_save_b));
                }
        }

        RUN(init_alnmem(m));

        /* Compute consistency bonus for all merge types */
        m->consistency = NULL;
        {
                int dp_row_node, dp_col_node, dp_rows, dp_cols;
                int have_consistency = (msa->poar_consistency != NULL)
                        || (msa->consistency_table != NULL);

                if(have_consistency){
                        if(msa->nsip[a] == 1 && msa->nsip[b] == 1){
                                if(len_a < len_b){
                                        dp_row_node = a; dp_rows = len_a;
                                        dp_col_node = b; dp_cols = len_b;
                                }else{
                                        dp_row_node = b; dp_rows = len_b;
                                        dp_col_node = a; dp_cols = len_a;
                                }
                        }else if(msa->nsip[a] == 1){
                                dp_row_node = b; dp_rows = len_b;
                                dp_col_node = a; dp_cols = len_a;
                        }else if(msa->nsip[b] == 1){
                                dp_row_node = a; dp_rows = len_a;
                                dp_col_node = b; dp_cols = len_b;
                        }else{
                                if(len_a < len_b){
                                        dp_row_node = a; dp_rows = len_a;
                                        dp_col_node = b; dp_cols = len_b;
                                }else{
                                        dp_row_node = b; dp_rows = len_b;
                                        dp_col_node = a; dp_cols = len_a;
                                }
                        }

                        if(msa->poar_consistency != NULL){
                                RUN(poar_consistency_get_bonus(
                                        (struct poar_consistency_ctx*)msa->poar_consistency,
                                        msa, dp_row_node, dp_rows, dp_col_node, dp_cols,
                                        &m->consistency));
                        }else{
                                RUN(anchor_consistency_get_bonus_profile(
                                        (struct consistency_table*)msa->consistency_table,
                                        msa, dp_row_node, dp_rows, dp_col_node, dp_cols,
                                        &m->consistency));
                        }
                }
        }

        /* Allocate best_path buffer */
        MMALLOC(best_path, sizeof(int) * m->alloc_path_len);

        /* Multi-trial alignment */
        for(k = 0; k < n_trials; k++){
                float sp = 0.0F;

                /* Re-initialize DP state */
                {
                        int _g = MACRO_MAX(len_a, len_b) + 2;
                        int _i;
                        for(_i = 0; _i < _g; _i++){
                                m->path[_i] = -1;
                        }
                }
                m->starta = 0;
                m->startb = 0;
                m->enda = len_a;
                m->endb = len_b;
                m->len_a = len_a;
                m->len_b = len_b;
                m->f[0].a = 0.0F;
                m->f[0].ga = -FLT_MAX;
                m->f[0].gb = -FLT_MAX;
                m->b[0].a = 0.0F;
                m->b[0].ga = -FLT_MAX;
                m->b[0].gb = -FLT_MAX;
                m->margin_sum = 0.0F;
                m->margin_count = 0;

                if(k == 0){
                        m->flip_threshold = 0.0F;
                        m->flip_trial = 0;
                }else{
                        m->flip_threshold = avg_margin;
                        m->flip_trial = k;
                        m->flip_stride = n_trials - 1;
                        m->flip_counter = 0;
                }

                /* Dispatch alignment (handle seq/profile combinations + mirroring) */
                if(msa->nsip[a] == 1){
                        if(msa->nsip[b] == 1){
                                if(len_a < len_b){
                                        m->seq1 = msa->sequences[a]->s;
                                        m->seq2 = msa->sequences[b]->s;
                                        m->prof1 = NULL;
                                        m->prof2 = NULL;
                                        aln_runner(m);
                                }else{
                                        m->enda = len_b;
                                        m->endb = len_a;
                                        m->len_a = len_b;
                                        m->len_b = len_a;
                                        m->seq1 = msa->sequences[b]->s;
                                        m->seq2 = msa->sequences[a]->s;
                                        m->prof1 = NULL;
                                        m->prof2 = NULL;
                                        aln_runner(m);
                                        RUN(mirror_path_n(m, len_a, len_b));
                                        m->len_a = len_a;
                                        m->len_b = len_b;
                                }
                        }else{
                                m->enda = len_b;
                                m->endb = len_a;
                                m->len_a = len_b;
                                m->len_b = len_a;
                                m->seq1 = NULL;
                                m->seq2 = msa->sequences[a]->s;
                                m->prof1 = t->profile[b];
                                m->prof2 = NULL;
                                m->sip = msa->nsip[b];
                                aln_runner(m);
                                RUN(mirror_path_n(m, len_a, len_b));
                                m->len_a = len_a;
                                m->len_b = len_b;
                        }
                }else{
                        if(msa->nsip[b] == 1){
                                m->seq1 = NULL;
                                m->seq2 = msa->sequences[b]->s;
                                m->prof1 = t->profile[a];
                                m->prof2 = NULL;
                                m->sip = msa->nsip[a];
                                aln_runner(m);
                        }else{
                                if(len_a < len_b){
                                        m->seq1 = NULL;
                                        m->seq2 = NULL;
                                        m->prof1 = t->profile[a];
                                        m->prof2 = t->profile[b];
                                        aln_runner(m);
                                }else{
                                        m->enda = len_b;
                                        m->endb = len_a;
                                        m->len_a = len_b;
                                        m->len_b = len_a;
                                        m->seq1 = NULL;
                                        m->seq2 = NULL;
                                        m->prof1 = t->profile[b];
                                        m->prof2 = t->profile[a];
                                        aln_runner(m);
                                        RUN(mirror_path_n(m, len_a, len_b));
                                        m->len_a = len_a;
                                        m->len_b = len_b;
                                }
                        }
                }

                /* Convert raw Hirschberg path to 0/1/2 format with gap info */
                RUN(add_gap_info_to_path_n(m));

                /* Score this trial */
                RUN(compute_sp_score(msa, m->ap, m->path,
                                     msa->sip[a], msa->nsip[a],
                                     msa->sip[b], msa->nsip[b], &sp));

                if(sp > best_sp){
                        best_sp = sp;
                        memcpy(best_path, m->path,
                               sizeof(int) * (m->path[0] + 2));
                }

                /* After baseline, compute avg_margin for flip threshold */
                if(k == 0 && m->margin_count > 0){
                        avg_margin = m->margin_sum / (float)m->margin_count;
                }
        }

        /* Install best path */
        memcpy(m->path, best_path, sizeof(int) * (best_path[0] + 2));

        /* Free consistency bonus if allocated */
        if(m->consistency){
                sparse_bonus_free(m->consistency);
                m->consistency = NULL;
        }

        /* Store confidence */
        t->list[task_id]->confidence = best_sp;

        vsm_unshift_profiles(t, a, b, &vsm_save_a, &vsm_save_b);

        /* Restore original aln_param for profile update */
        m->ap = orig_ap;

        /* Merge profiles */
        MMALLOC(tmp, sizeof(float) * 64 * (m->path[0] + 2));
        if(task_id != t->n_tasks - 1){
                update_n(t->profile[a], t->profile[b], tmp, m->ap, m->path,
                         msa->nsip[a], msa->nsip[b]);
        }

        MFREE(t->profile[a]);
        MFREE(t->profile[b]);
        t->profile[a] = NULL;
        t->profile[b] = NULL;

        t->profile[c] = tmp;
        RUN(make_seq(msa, a, b, m->path));

        msa->plen[c] = m->path[0];
        msa->nsip[c] = msa->nsip[a] + msa->nsip[b];

        MREALLOC(msa->sip[c], sizeof(int) * (msa->nsip[a] + msa->nsip[b]));
        g = 0;
        for(j = msa->nsip[a]; j--;){
                msa->sip[c][g] = msa->sip[a][j];
                g++;
        }
        for(j = msa->nsip[b]; j--;){
                msa->sip[c][g] = msa->sip[b][j];
                g++;
        }

        MFREE(best_path);
        return OK;
ERROR:
        if(best_path){
                MFREE(best_path);
        }
        /* Leave the profiles and parameters as they were before the merge. */
        if(m->consistency){
                sparse_bonus_free(m->consistency);
                m->consistency = NULL;
        }
        vsm_unshift_profiles(t, a, b, &vsm_save_a, &vsm_save_b);
        m->ap = orig_ap;
        return FAIL;
}

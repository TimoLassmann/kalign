#include "tldevel.h"

#include <string.h>

#include "msa_struct.h"
#include "task.h"
#include "aln_param.h"
#include "aln_struct.h"
#include "aln_mem.h"
#include "aln_setup.h"
#include "aln_controller.h"
#include "weave_alignment.h"
#include "sp_score.h"
#include "aln_run.h"
#include "anchor_consistency.h"
#include "aln_vsm.h"

#define ALN_REFINE_IMPORT
#include "aln_refine.h"

#define REFINE_N_TRIALS 5

#include <stdatomic.h>
#include "threadpool.h"

struct refine_ctx {
        struct msa* msa;
        struct aln_param* ap;
        struct aln_tasks* t;
        int refine_mode;
        float threshold;
        int* parent;
        atomic_int* pending;
        atomic_int error;
};
struct refine_arg {
        struct refine_ctx* ctx;
        int node;
};
static void refine_run_chain(struct refine_ctx* c, int node);
static int refine_tree(struct msa* msa, struct aln_param* ap, struct aln_tasks* t,
                       int refine_mode, float threshold);

static int refine_or_replay(struct msa* msa, struct aln_param* ap, struct aln_tasks* t,
                            int refine_mode, float threshold, int task_id);
static int refine_edge(struct msa* msa, struct aln_param* ap, struct aln_tasks* t, int task_id);
static int replay_edge(struct msa* msa, struct aln_param* ap, struct aln_tasks* t, int task_id);
static int convert_raw_path(struct aln_mem* m);
static int compute_confidence_threshold(struct aln_tasks* t, float* threshold);

/* Leaf weight for make_profile_n — always 1.0 since sequence weighting
   is handled via balanced freq-count merging in update_n. */
static float leaf_weight(struct msa* msa, int node)
{
        (void)msa; (void)node;
        return 1.0f;
}

int refine_alignment(struct msa* msa, struct aln_param* ap, struct aln_tasks* t, int refine_mode)
{ KPERF_PHASE_REGION(KP_REFINE, KPH_REFINE);
        int i;
        float threshold = 0.0F;

        if(refine_mode == KALIGN_REFINE_NONE){
                return OK;
        }

        RUN(sort_tasks(t, TASK_ORDER_TREE));

        if(refine_mode == KALIGN_REFINE_CONFIDENT){
                RUN(compute_confidence_threshold(t, &threshold));
        }

        /* Reset alignment state: zero all gaps and restore leaf sip/nsip */
        RUN(clean_aln(msa));

        /* Clear stale profiles from the initial alignment */
        for(i = 0; i < msa->num_profiles; i++){
                if(t->profile[i]){
                        MFREE(t->profile[i]);
                        t->profile[i] = NULL;
                }
        }

        /* Process all edges bottom-up using progressive profile building.
           Each edge builds profiles from its children (identical to original
           progressive alignment), then refined edges run multi-trial DP.
           An edge depends only on its own subtree, so with a thread pool
           the two subtrees of every node are refined concurrently (as in
           create_msa_tree); the result does not depend on the order. */
        {
                for(i = 0; i < t->n_tasks; i++){
                        ASSERT(t->list[i]->c == msa->numseq + i, "Tasks not in tree order");
                }
                RUN(refine_tree(msa, ap, t, refine_mode, threshold));
        }

        return OK;
ERROR:
        return FAIL;
}

static int refine_or_replay(struct msa* msa, struct aln_param* ap, struct aln_tasks* t,
                            int refine_mode, float threshold, int task_id)
{
        int should_refine = 0;
        if(refine_mode == KALIGN_REFINE_ALL){
                should_refine = 1;
        }else if(refine_mode == KALIGN_REFINE_CONFIDENT){
                should_refine = (t->list[task_id]->confidence <= threshold);
        }
        if(should_refine){
                return refine_edge(msa, ap, t, task_id);
        }
        return replay_edge(msa, ap, t, task_id);
}

static void refine_task(void* arg)
{
        struct refine_arg* a = (struct refine_arg*)arg;
        refine_run_chain(a->ctx, a->node);
}

/* Refine node, then keep going up the tree for as long as this thread
   finishes the last unfinished child of the parent.  Nothing here waits
   on another tree task, so stack use does not depend on tree depth. */
static void refine_run_chain(struct refine_ctx* c, int node)
{
        while(node >= 0){
                int p;
                if(atomic_load(&c->error)){
                        return;
                }
                if(refine_or_replay(c->msa, c->ap, c->t, c->refine_mode, c->threshold, node) != OK){
                        atomic_store(&c->error, 1);
                        return;
                }
                p = c->parent[node];
                if(p < 0){
                        return;
                }
                if(atomic_fetch_sub_explicit(&c->pending[p], 1, memory_order_acq_rel) != 1){
                        return;
                }
                node = p;
        }
}

static int refine_tree(struct msa* msa, struct aln_param* ap, struct aln_tasks* t,
                       int refine_mode, float threshold)
{
        struct refine_ctx c;
        struct refine_arg* args = NULL;
        int* ready = NULL;
        int n = t->n_tasks;
        int n_ready = 0;
        int i;

        c.msa = msa;
        c.ap = ap;
        c.t = t;
        c.refine_mode = refine_mode;
        c.threshold = threshold;
        c.parent = NULL;
        c.pending = NULL;
        atomic_init(&c.error, 0);

        MMALLOC(c.parent, sizeof(int) * n);
        MMALLOC(c.pending, sizeof(atomic_int) * n);
        MMALLOC(ready, sizeof(int) * n);
        for(i = 0; i < n; i++){
                c.parent[i] = -1;
                atomic_init(&c.pending[i], 0);
        }
        for(i = 0; i < n; i++){
                int ch[2];
                int j;
                ch[0] = t->list[i]->a - msa->numseq;
                ch[1] = t->list[i]->b - msa->numseq;
                for(j = 0; j < 2; j++){
                        if(ch[j] >= 0){
                                c.parent[ch[j]] = i;
                                atomic_fetch_add(&c.pending[i], 1);
                        }
                }
        }
        for(i = 0; i < n; i++){
                if(atomic_load(&c.pending[i]) == 0){
                        ready[n_ready++] = i;
                }
        }
        if(ap->pool){
                tp_group_t* g = tp_group_create(ap->pool);
                MMALLOC(args, sizeof(struct refine_arg) * n_ready);
                for(i = 0; i < n_ready; i++){
                        args[i].ctx = &c;
                        args[i].node = ready[i];
                        if(!g || tp_group_submit(g, refine_task, &args[i]) != 0){
                                refine_run_chain(&c, ready[i]);
                        }
                }
                if(g){
                        tp_group_wait(g);
                        tp_group_destroy(g);
                }
                MFREE(args);
                args = NULL;
        }else{
                for(i = 0; i < n_ready; i++){
                        refine_run_chain(&c, ready[i]);
                }
        }
        if(atomic_load(&c.error)){
                ERROR_MSG("Refinement failed");
        }
        MFREE(ready);
        MFREE(c.pending);
        MFREE(c.parent);
        return OK;
ERROR:
        if(args){
                MFREE(args);
        }
        if(ready){
                MFREE(ready);
        }
        if(c.pending){
                MFREE(c.pending);
        }
        if(c.parent){
                MFREE(c.parent);
        }
        return FAIL;
}

/* Process an edge with multi-trial DP: builds progressive profiles
   (identical to do_align), runs K alignment trials (one deterministic
   baseline plus threshold-based alternatives), keeps the best by SP score. */
/* One alignment trial on ml: trial 0 is the deterministic baseline, trial
   k > 0 flips meetups whose margin is below avg_margin (every
   (n_trials-1)-th one, offset k).  Leaves the 0/1/2 path in ml->path and
   its sum-of-pairs score in *sp. */
static int run_trial(struct msa* msa, struct aln_param* ap, struct aln_mem* ml,
                     float* prof_a, float* prof_b, int a, int b, int len_a, int len_b,
                     int k, int n_trials, float avg_margin, float* sp)
{
        /* Re-initialize DP state before each trial. The Hirschberg
           recursion only writes matched positions to the path array;
           unmatched positions must be -1. After convert_raw_path swaps
           path/tmp_path, the buffer contains stale 0/1/2 values. */
        {
                int _g = MACRO_MAX(len_a, len_b) + 2;
                int _i;
                for(_i = 0; _i < _g; _i++){
                        ml->path[_i] = -1;
                }
        }
        ml->starta = 0;
        ml->startb = 0;
        ml->enda = len_a;
        ml->endb = len_b;
        ml->len_a = len_a;
        ml->len_b = len_b;
        ml->f[0].a = 0.0F;
        ml->f[0].ga = -FLT_MAX;
        ml->f[0].gb = -FLT_MAX;
        ml->b[0].a = 0.0F;
        ml->b[0].ga = -FLT_MAX;
        ml->b[0].gb = -FLT_MAX;
        ml->margin_sum = 0.0F;
        ml->margin_count = 0;

        if(k == 0){
                ml->flip_threshold = 0.0F;
                ml->flip_trial = 0;
        }else{
                ml->flip_threshold = avg_margin;
                ml->flip_trial = k;
                ml->flip_stride = n_trials - 1;
                ml->flip_counter = 0;
        }

        RUN(dispatch_alignment(msa, ml, prof_a, prof_b, a, b, len_a, len_b));

        /* Convert raw DP path to 0/1/2 format with gap info */
        RUN(convert_raw_path(ml));

        /* Score this candidate with sum-of-pairs */
        *sp = 0.0F;
        RUN(compute_sp_score(msa, ap, ml->path,
                             msa->sip[a], msa->nsip[a],
                             msa->sip[b], msa->nsip[b], sp));
        return OK;
ERROR:
        return FAIL;
}

struct trial_arg {
        struct msa* msa;
        struct aln_param* ap;
        struct aln_mem* ml;
        float* prof_a;
        float* prof_b;
        int a, b, len_a, len_b;
        int k, n_trials;
        float avg_margin;
        float sp;
        int ret;
};

static void trial_task(void* arg)
{
        struct trial_arg* x = (struct trial_arg*)arg;
        x->ret = run_trial(x->msa, x->ap, x->ml, x->prof_a, x->prof_b, x->a, x->b,
                           x->len_a, x->len_b, x->k, x->n_trials, x->avg_margin, &x->sp);
}

/* Trials 1..n_trials-1 concurrently, each on a copy of base's settings
   with its own DP memory (the consistency bonus is shared read-only).
   Then the serial loop's selection, in trial order: a trial replaces the
   best only if its score is strictly higher. */
static int run_trials_parallel(struct msa* msa, struct aln_param* ap, threadpool_t* pool,
                               struct aln_mem* base, float* prof_a, float* prof_b,
                               int a, int b, int len_a, int len_b, int n_trials, float avg_margin,
                               float* best_sp, float* best_margin_sum, int* best_margin_count,
                               int* best_path)
{
        struct trial_arg* x = NULL;
        tp_group_t* g = NULL;
        int k;
        const int n = n_trials - 1;

        MMALLOC(x, sizeof(struct trial_arg) * n);
        for(k = 0; k < n; k++){
                x[k].ml = NULL;
        }
        for(k = 0; k < n; k++){
                struct aln_mem* m = NULL;
                RUN(alloc_aln_mem(&m, 256));
                x[k].ml = m;
                m->pool = base->pool;
                m->run_parallel = base->run_parallel;
                m->ap = base->ap;
                m->mode = base->mode;
                m->len_a = len_a;
                m->len_b = len_b;
                RUN(init_alnmem(m));
                m->consistency = base->consistency;
                x[k].msa = msa;
                x[k].ap = ap;
                x[k].prof_a = prof_a;
                x[k].prof_b = prof_b;
                x[k].a = a;
                x[k].b = b;
                x[k].len_a = len_a;
                x[k].len_b = len_b;
                x[k].k = k + 1;
                x[k].n_trials = n_trials;
                x[k].avg_margin = avg_margin;
                x[k].sp = 0.0F;
                x[k].ret = FAIL;
        }

        g = tp_group_create(pool);
        for(k = 0; k < n; k++){
                if(!g || tp_group_submit(g, trial_task, &x[k]) != 0){
                        trial_task(&x[k]);
                }
        }
        if(g){
                tp_group_wait(g);
                tp_group_destroy(g);
        }

        for(k = 0; k < n; k++){
                ASSERT(x[k].ret == OK, "Refinement trial %d failed", k + 1);
                if(x[k].sp > *best_sp){
                        *best_sp = x[k].sp;
                        *best_margin_sum = x[k].ml->margin_sum;
                        *best_margin_count = x[k].ml->margin_count;
                        memcpy(best_path, x[k].ml->path,
                               sizeof(int) * (x[k].ml->path[0] + 2));
                }
        }
        for(k = 0; k < n; k++){
                x[k].ml->consistency = NULL;
                free_aln_mem(x[k].ml);
        }
        MFREE(x);
        return OK;
ERROR:
        if(x){
                for(k = 0; k < n; k++){
                        if(x[k].ml){
                                x[k].ml->consistency = NULL;
                                free_aln_mem(x[k].ml);
                        }
                }
                MFREE(x);
        }
        return FAIL;
}

int refine_edge(struct msa* msa, struct aln_param* ap, struct aln_tasks* t, int task_id)
{
        struct aln_mem* ml = NULL;
        float* tmp = NULL;
        int* best_path = NULL;
        int a, b, c;
        int len_a, len_b;
        int k, j, g;
        float best_sp = -FLT_MAX;
        float best_margin_sum = 0.0F;
        int best_margin_count = 0;
        float avg_margin = 0.0F;
        int n_trials = REFINE_N_TRIALS;

        a = t->list[task_id]->a;
        b = t->list[task_id]->b;
        c = t->list[task_id]->c;

        /* Distance-dependent parameter scaling (must match do_align) */
        struct aln_param* orig_ap = ap;
        struct aln_param scaled_ap;
        float* vsm_save_a = NULL;
        float* vsm_save_b = NULL;
        float gap_scale = compute_gap_scale(msa, ap, a, b);
        float subm_off = compute_subm_offset(msa, ap, a, b);
        if(gap_scale < 1.0f || subm_off > 0.0f){
                scaled_ap = *ap;
                scaled_ap.gpo *= gap_scale;
                scaled_ap.gpe *= gap_scale;
                scaled_ap.tgpe *= gap_scale;
                scaled_ap.subm_offset = subm_off;
                ap = &scaled_ap;
        }

        /* Build profiles progressively (same as do_align) */
        if(msa->nsip[a] == 1){
                len_a = msa->sequences[a]->len;
                RUN(make_profile_n(ap, msa->sequences[a]->s, len_a, leaf_weight(msa,a), &t->profile[a]));
        }else{
                len_a = msa->plen[a];
                RUN(set_gap_penalties_n(t->profile[a], len_a, msa->nsip[b]));
        }

        if(msa->nsip[b] == 1){
                len_b = msa->sequences[b]->len;
                RUN(make_profile_n(ap, msa->sequences[b]->s, len_b, leaf_weight(msa,b), &t->profile[b]));
        }else{
                len_b = msa->plen[b];
                RUN(set_gap_penalties_n(t->profile[b], len_b, msa->nsip[a]));
        }

        {
                float vsm_off = 0.0f;
                RUN(vsm_merge_offset(msa, t, ap, a, b, len_a, len_b, &vsm_off));
                if(vsm_off > 0.0f){
                        if(ap != &scaled_ap){
                                scaled_ap = *ap;
                                ap = &scaled_ap;
                        }
                        scaled_ap.subm_offset = vsm_off;
                        RUN(vsm_shift_profiles(t, a, b, len_a, len_b, vsm_off, &vsm_save_a, &vsm_save_b));
                }
        }

        /* Allocate alignment memory */
        RUN(alloc_aln_mem(&ml, 256));
        ml->pool = orig_ap->pool;
        ml->run_parallel = (orig_ap->pool != NULL);
        ml->ap = ap;
        ml->mode = ALN_MODE_FULL;
        ml->len_a = len_a;
        ml->len_b = len_b;
        ml->margin_sum = 0.0F;
        ml->margin_count = 0;
        RUN(init_alnmem(ml));

        /* Allocate best_path buffer */
        MMALLOC(best_path, sizeof(int) * ml->alloc_path_len);

        /* Compute consistency bonus for all merge types */
        ml->consistency = NULL;
        {
                struct consistency_table* ct = (struct consistency_table*)msa->consistency_table;
                if(ct != NULL){
                        int dp_row_node, dp_col_node, dp_rows, dp_cols;
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
                        RUN(anchor_consistency_get_bonus_profile(ct, msa,
                                dp_row_node, dp_rows, dp_col_node, dp_cols,
                                &ml->consistency));
                }
        }

        /* If adaptive budget, allocate flip_margins for baseline trial
           to record per-meetup margins for uncertainty analysis */
        if(ap->adaptive_budget){
                int est = MACRO_MIN(len_a, len_b) + 1;
                if(est < 64) est = 64;
                MMALLOC(ml->flip_margins, sizeof(float) * est);
                ml->flip_margin_alloc = est;
        }

        /* Multi-trial alignment: trial 0 is deterministic (baseline),
           trials 1..n_trials-1 use increasing flip thresholds.  Trials
           1..n_trials-1 depend only on avg_margin from trial 0, so with a
           thread pool they run concurrently, each on its own aln_mem; the
           best is then picked in trial order exactly as in the serial loop. */
        for(k = 0; k < n_trials; k++){
                float sp = 0.0F;

                if(k == 1 && orig_ap->pool != NULL && n_trials > 2){
                        RUN(run_trials_parallel(msa, ap, orig_ap->pool, ml, t->profile[a], t->profile[b],
                                                a, b, len_a, len_b, n_trials, avg_margin,
                                                &best_sp, &best_margin_sum, &best_margin_count, best_path));
                        break;
                }
                RUN(run_trial(msa, ap, ml, t->profile[a], t->profile[b], a, b, len_a, len_b,
                              k, n_trials, avg_margin, &sp));

                if(sp > best_sp){
                        best_sp = sp;
                        best_margin_sum = ml->margin_sum;
                        best_margin_count = ml->margin_count;
                        memcpy(best_path, ml->path,
                               sizeof(int) * (ml->path[0] + 2));
                }

                /* After baseline trial, compute avg_margin and adaptive budget */
                if(k == 0){
                        if(ml->margin_count > 0){
                                avg_margin = ml->margin_sum / (float)ml->margin_count;
                        }

                        /* Adaptive budget: count very uncertain meetups */
                        if(ap->adaptive_budget && ml->flip_margins && ml->margin_count > 0){
                                int n_very_uncertain = 0;
                                int m_i;
                                float vu_threshold = avg_margin * 0.25F;
                                for(m_i = 0; m_i < ml->margin_count; m_i++){
                                        if(ml->flip_margins[m_i] < vu_threshold){
                                                n_very_uncertain++;
                                        }
                                }
                                {
                                        float frac = (float)n_very_uncertain / (float)ml->margin_count;
                                        n_trials = 1 + (int)(7.0F * frac + 0.5F);
                                }
                        }

                        /* Free flip_margins after extracting stats */
                        if(ml->flip_margins){
                                MFREE(ml->flip_margins);
                                ml->flip_margins = NULL;
                        }
                }
        }

        /* Restore best path into ml */
        memcpy(ml->path, best_path,
               sizeof(int) * (best_path[0] + 2));

        /* Free consistency bonus */
        if(ml->consistency){
                sparse_bonus_free(ml->consistency);
                ml->consistency = NULL;
        }

        /* Update confidence from best trial */
        if(best_margin_count > 0){
                t->list[task_id]->confidence = best_margin_sum / (float)best_margin_count;
        }else{
                t->list[task_id]->confidence = 0.0F;
        }

        vsm_unshift_profiles(t, a, b, &vsm_save_a, &vsm_save_b);

        /* Restore original aln_param for profile update (unscaled base penalties) */
        ap = orig_ap;

        /* Merge profiles for downstream edges */
        MMALLOC(tmp, sizeof(float) * 64 * (ml->path[0] + 2));
        if(task_id != t->n_tasks - 1){
                update_n(t->profile[a], t->profile[b], tmp, ap, ml->path, msa->nsip[a], msa->nsip[b]);
        }
        t->profile[c] = tmp;

        /* Update gap arrays for all member sequences */
        RUN(make_seq(msa, a, b, ml->path));

        msa->plen[c] = ml->path[0];
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

        /* Free child profiles */
        MFREE(t->profile[a]);
        t->profile[a] = NULL;
        MFREE(t->profile[b]);
        t->profile[b] = NULL;

        /* Cleanup */
        MFREE(best_path);
        free_aln_mem(ml);
        return OK;
ERROR:
        if(best_path){
                MFREE(best_path);
        }
        /* Leave the profiles as they were before this edge. */
        if(ml && ml->consistency){
                sparse_bonus_free(ml->consistency);
                ml->consistency = NULL;
        }
        vsm_unshift_profiles(t, a, b, &vsm_save_a, &vsm_save_b);
        free_aln_mem(ml);
        return FAIL;
}

/* Process an edge using standard progressive profiles (identical to
   the original do_align in aln_run.c). Used for high-confidence edges
   in REFINE_CONFIDENT mode and as the default replay path. */
int replay_edge(struct msa* msa, struct aln_param* ap, struct aln_tasks* t, int task_id)
{
        struct aln_mem* ml = NULL;
        float* tmp = NULL;
        int a, b, c;
        int len_a, len_b;
        int j, g;

        a = t->list[task_id]->a;
        b = t->list[task_id]->b;
        c = t->list[task_id]->c;

        /* Distance-dependent parameter scaling (must match do_align) */
        struct aln_param* orig_ap = ap;
        struct aln_param scaled_ap;
        float* vsm_save_a = NULL;
        float* vsm_save_b = NULL;
        float gap_scale = compute_gap_scale(msa, ap, a, b);
        float subm_off = compute_subm_offset(msa, ap, a, b);
        if(gap_scale < 1.0f || subm_off > 0.0f){
                scaled_ap = *ap;
                scaled_ap.gpo *= gap_scale;
                scaled_ap.gpe *= gap_scale;
                scaled_ap.tgpe *= gap_scale;
                scaled_ap.subm_offset = subm_off;
                ap = &scaled_ap;
        }

        /* Build profiles (standard progressive approach) */
        if(msa->nsip[a] == 1){
                len_a = msa->sequences[a]->len;
                RUN(make_profile_n(ap, msa->sequences[a]->s, len_a, leaf_weight(msa,a), &t->profile[a]));
        }else{
                len_a = msa->plen[a];
                RUN(set_gap_penalties_n(t->profile[a], len_a, msa->nsip[b]));
        }

        if(msa->nsip[b] == 1){
                len_b = msa->sequences[b]->len;
                RUN(make_profile_n(ap, msa->sequences[b]->s, len_b, leaf_weight(msa,b), &t->profile[b]));
        }else{
                len_b = msa->plen[b];
                RUN(set_gap_penalties_n(t->profile[b], len_b, msa->nsip[a]));
        }

        {
                float vsm_off = 0.0f;
                RUN(vsm_merge_offset(msa, t, ap, a, b, len_a, len_b, &vsm_off));
                if(vsm_off > 0.0f){
                        if(ap != &scaled_ap){
                                scaled_ap = *ap;
                                ap = &scaled_ap;
                        }
                        scaled_ap.subm_offset = vsm_off;
                        RUN(vsm_shift_profiles(t, a, b, len_a, len_b, vsm_off, &vsm_save_a, &vsm_save_b));
                }
        }

        /* Allocate alignment memory and run alignment */
        RUN(alloc_aln_mem(&ml, 256));
        ml->pool = orig_ap->pool;
        ml->run_parallel = (orig_ap->pool != NULL);
        ml->ap = ap;
        ml->mode = ALN_MODE_FULL;
        ml->len_a = len_a;
        ml->len_b = len_b;
        ml->margin_sum = 0.0F;
        ml->margin_count = 0;
        RUN(init_alnmem(ml));

        /* Compute consistency bonus for all merge types */
        ml->consistency = NULL;
        {
                struct consistency_table* ct = (struct consistency_table*)msa->consistency_table;
                if(ct != NULL){
                        int dp_row_node, dp_col_node, dp_rows, dp_cols;
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
                        RUN(anchor_consistency_get_bonus_profile(ct, msa,
                                dp_row_node, dp_rows, dp_col_node, dp_cols,
                                &ml->consistency));
                }
        }

        RUN(dispatch_alignment(msa, ml, t->profile[a], t->profile[b], a, b, len_a, len_b));

        /* Free consistency bonus */
        if(ml->consistency){
                sparse_bonus_free(ml->consistency);
                ml->consistency = NULL;
        }

        /* Store alignment confidence */
        if(ml->margin_count > 0){
                t->list[task_id]->confidence = ml->margin_sum / (float)ml->margin_count;
        }else{
                t->list[task_id]->confidence = 0.0F;
        }

        RUN(convert_raw_path(ml));

        vsm_unshift_profiles(t, a, b, &vsm_save_a, &vsm_save_b);

        /* Restore original aln_param for profile update (unscaled base penalties) */
        ap = orig_ap;

        /* Merge profiles for downstream edges */
        MMALLOC(tmp, sizeof(float) * 64 * (ml->path[0] + 2));
        if(task_id != t->n_tasks - 1){
                update_n(t->profile[a], t->profile[b], tmp, ap, ml->path, msa->nsip[a], msa->nsip[b]);
        }
        t->profile[c] = tmp;

        /* Update gap arrays for all member sequences */
        RUN(make_seq(msa, a, b, ml->path));

        msa->plen[c] = ml->path[0];
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

        /* Free consumed profiles */
        MFREE(t->profile[a]);
        t->profile[a] = NULL;
        MFREE(t->profile[b]);
        t->profile[b] = NULL;

        free_aln_mem(ml);
        return OK;
ERROR:
        /* Leave the profiles as they were before this edge. */
        if(ml && ml->consistency){
                sparse_bonus_free(ml->consistency);
                ml->consistency = NULL;
        }
        vsm_unshift_profiles(t, a, b, &vsm_save_a, &vsm_save_b);
        free_aln_mem(ml);
        return FAIL;
}

/* Dispatch alignment based on group types (seq-seq, seq-profile, profile-profile).
   Handles the shorter-first convention and path mirroring. */
int dispatch_alignment(struct msa* msa, struct aln_mem* ml, float* prof_a, float* prof_b, int a, int b, int len_a, int len_b)
{
        if(msa->nsip[a] == 1){
                if(msa->nsip[b] == 1){
                        /* seq vs seq */
                        if(len_a < len_b){
                                ml->seq1 = msa->sequences[a]->s;
                                ml->seq2 = msa->sequences[b]->s;
                                ml->prof1 = NULL;
                                ml->prof2 = NULL;
                                aln_runner(ml);
                        }else{
                                ml->enda = len_b;
                                ml->endb = len_a;
                                ml->len_a = len_b;
                                ml->len_b = len_a;
                                ml->seq1 = msa->sequences[b]->s;
                                ml->seq2 = msa->sequences[a]->s;
                                ml->prof1 = NULL;
                                ml->prof2 = NULL;
                                aln_runner(ml);
                                RUN(mirror_path_n(ml, len_a, len_b));
                                ml->len_a = len_a;
                                ml->len_b = len_b;
                        }
                }else{
                        /* seq a vs profile b: profile must be prof1 */
                        ml->enda = len_b;
                        ml->endb = len_a;
                        ml->len_a = len_b;
                        ml->len_b = len_a;
                        ml->seq1 = NULL;
                        ml->seq2 = msa->sequences[a]->s;
                        ml->prof1 = prof_b;
                        ml->prof2 = NULL;
                        ml->sip = msa->nsip[b];
                        aln_runner(ml);
                        RUN(mirror_path_n(ml, len_a, len_b));
                        ml->len_a = len_a;
                        ml->len_b = len_b;
                }
        }else{
                if(msa->nsip[b] == 1){
                        /* profile a vs seq b */
                        ml->seq1 = NULL;
                        ml->seq2 = msa->sequences[b]->s;
                        ml->prof1 = prof_a;
                        ml->prof2 = NULL;
                        ml->sip = msa->nsip[a];
                        aln_runner(ml);
                }else{
                        /* profile vs profile */
                        if(len_a < len_b){
                                ml->seq1 = NULL;
                                ml->seq2 = NULL;
                                ml->prof1 = prof_a;
                                ml->prof2 = prof_b;
                                aln_runner(ml);
                        }else{
                                ml->enda = len_b;
                                ml->endb = len_a;
                                ml->len_a = len_b;
                                ml->len_b = len_a;
                                ml->seq1 = NULL;
                                ml->seq2 = NULL;
                                ml->prof1 = prof_b;
                                ml->prof2 = prof_a;
                                aln_runner(ml);
                                RUN(mirror_path_n(ml, len_a, len_b));
                                ml->len_a = len_a;
                                ml->len_b = len_b;
                        }
                }
        }
        return OK;
ERROR:
        return FAIL;
}

/* Convert raw Hirschberg path (path[i] = B position or -1) to the
   0/1/2 format with gap info bits that update_n and make_seq expect.

   This replaces add_gap_info_to_path_n for refinement trials. The
   original has a latent bug: it tracks the last B position via
   `b = path[i]`, which resets b to -1 on gap-in-B entries. When a
   later match skips B positions, the gap-in-A fill is skipped because
   `b == -1`. Stochastic Hirschberg sampling triggers this by producing
   interleaved gap patterns. This function tracks b_last correctly:
   only updated on actual matches (path[i] != -1). */
int convert_raw_path(struct aln_mem* m)
{ KPERF_REGION(KP_PATH_CONV);
        int* path = NULL;
        int* o_path = NULL;
        int* tmp_path = NULL;
        int i, j, a;
        int b_last;
        int len_a;
        int len_b;

        len_a = m->len_a;
        len_b = m->len_b;
        path = m->path;
        o_path = m->tmp_path;

        for(i = 0; i < len_a + len_b + 2; i++){
                o_path[i] = 0;
        }

        j = 1;
        b_last = 0;  /* last consumed B position (1-indexed), 0 = none */

        for(i = 1; i <= len_a; i++){
                if(path[i] == -1){
                        o_path[j] = 2;
                        j++;
                }else{
                        /* Fill gap-in-A for B positions b_last+1..path[i]-1 */
                        for(a = b_last + 1; a < path[i]; a++){
                                o_path[j] = 1;
                                j++;
                        }
                        o_path[j] = 0;
                        j++;
                        b_last = path[i];
                }
        }

        /* Trailing gap-in-A for remaining B positions */
        for(a = b_last + 1; a <= len_b; a++){
                o_path[j] = 1;
                j++;
        }

        o_path[0] = j - 1;
        o_path[j] = 3;

        /* Add gap info bits */
        i = 2;
        while(o_path[i] != 3){
                if ((o_path[i-1] & 3) && !(o_path[i] & 3)){
                        if(o_path[i-1] & 8){
                                o_path[i-1] += 8;
                        }else{
                                o_path[i-1] |= 16;
                        }
                }else if (!(o_path[i-1] & 3) && (o_path[i] & 3)){
                        o_path[i] |= 4;
                }else if ((o_path[i-1] & 1) && (o_path[i] & 1)){
                        o_path[i] |= 8;
                }else if ((o_path[i-1] & 2) && (o_path[i] & 2)){
                        o_path[i] |= 8;
                }
                i++;
        }

        /* Add terminal gap flags */
        i = 1;
        while(o_path[i] != 0){
                o_path[i] |= 32;
                i++;
        }
        i = o_path[0];
        while(o_path[i] != 0){
                o_path[i] |= 32;
                i--;
        }

        tmp_path = m->path;
        m->path = m->tmp_path;
        m->tmp_path = tmp_path;
        return OK;
}

int compute_confidence_threshold(struct aln_tasks* t, float* threshold)
{
        float* confidences = NULL;
        int n, i, j;
        float tmp;

        n = t->n_tasks;
        ASSERT(n > 0, "No tasks");

        MMALLOC(confidences, sizeof(float) * n);
        for(i = 0; i < n; i++){
                confidences[i] = t->list[i]->confidence;
        }

        /* Insertion sort for median */
        for(i = 1; i < n; i++){
                tmp = confidences[i];
                j = i - 1;
                while(j >= 0 && confidences[j] > tmp){
                        confidences[j + 1] = confidences[j];
                        j--;
                }
                confidences[j + 1] = tmp;
        }

        if(n % 2 == 0){
                *threshold = (confidences[n / 2 - 1] + confidences[n / 2]) / 2.0F;
        }else{
                *threshold = confidences[n / 2];
        }

        MFREE(confidences);
        return OK;
ERROR:
        if(confidences){
                MFREE(confidences);
        }
        return FAIL;
}

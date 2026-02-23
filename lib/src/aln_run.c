#include "tldevel.h"
#include "tlrng.h"

#include "msa_struct.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

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

#include <string.h>
#include <float.h>

#define ALN_RUN_IMPORT
#include "aln_run.h"

static void recursive_aln(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c);
/* static void recursive_aln_openMP(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c); */
/* static void recursive_aln_serial(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c); */

static int do_align(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id);
static int do_align_inline_refine(struct msa* msa, struct aln_tasks* t, struct aln_mem* m, int task_id, int n_trials);
static void recursive_aln_inline(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c, int n_trials);
/* static int do_align_serial(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id); */
/* static int do_score(struct msa* msa,struct aln_tasks* t,struct aln_mem* m, int task_id); */

/* static int SampleWithoutReplacement(struct rng_state* rng, int N, int n,int* samples); */
/* static int int_cmp(const void *a, const void *b); */

int create_msa_tree(struct msa* msa, struct aln_param* ap,struct aln_tasks* t)
{
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

#ifdef HAVE_OPENMP
#pragma omp parallel
#pragma omp single nowait
#endif
        recursive_aln(msa, t, ap, active, t->n_tasks-1);

        MFREE(active);
        return OK;
ERROR:
        if(active){
                MFREE(active);
        }
        return FAIL;
}


void recursive_aln(struct msa* msa, struct aln_tasks*t, struct aln_param* ap, uint8_t* active, int c)
{
        struct task* local_t = NULL;

        /* Follow left and right branch until I arrive at sequences / profiles
           ready to align.
        */
        int a;
        int b;
        local_t = t->list[c];

        a = local_t->a - msa->numseq;
        b = local_t->b - msa->numseq;

        if(!active[local_t->a] && local_t->a >= msa->numseq){
#ifdef HAVE_OPENMP
#pragma omp task shared(msa,t,ap,active) firstprivate(a)
#endif
                recursive_aln(msa, t, ap, active, a);
        }
        if(!active[local_t->b] && local_t->b >= msa->numseq){
#ifdef HAVE_OPENMP
#pragma omp task shared(msa,t,ap,active) firstprivate(b)
#endif
                recursive_aln(msa, t, ap, active, b);
        }
#ifdef HAVE_OPENMP
#pragma omp taskwait
#endif

        struct aln_mem* ml = NULL;

        alloc_aln_mem(&ml, 256);

        ml->ap = ap;
        ml->mode = ALN_MODE_FULL;
        do_align(msa,t,ml,c);

        active[local_t->a] = 0;
        active[local_t->b] = 0;
        active[local_t->c] = 1;
        /* LOG_MSG("Local: %d %d %d p:%d", local_t->a, local_t->b, local_t->c, local_t->p); */
        free_aln_mem(ml);
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
        if(amax <= 0.0f || msa->seq_distances == NULL){
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

        RUN(init_alnmem(m));

        m->margin_sum = 0.0F;
        m->margin_count = 0;
        m->consistency = NULL;
        m->consistency_stride = 0;

        /* Compute consistency bonus for all merge types */
        {
                struct consistency_table* ct = (struct consistency_table*)msa->consistency_table;
                if(ct != NULL){
                        int dp_row_node, dp_col_node, dp_rows, dp_cols;
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
                        RUN(anchor_consistency_get_bonus_profile(ct, msa,
                                dp_row_node, dp_rows, dp_col_node, dp_cols,
                                &m->consistency));
                        m->consistency_stride = dp_cols;
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
                MFREE(m->consistency);
                m->consistency = NULL;
                m->consistency_stride = 0;
        }

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
        return FAIL;
}

/* ---------------------------------------------------------------------------
 * Inline refinement: run N trials at each progressive merge, keep the best
 * by SP score. Single pass — no second replay pass needed.
 * ------------------------------------------------------------------------- */

int create_msa_tree_inline_refine(struct msa* msa, struct aln_param* ap,
                                  struct aln_tasks* t, int n_trials)
{
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

        /* Inline refine is sequential — multi-trial per edge isn't thread-safe */
        msa->run_parallel = 0;

        recursive_aln_inline(msa, t, ap, active, t->n_tasks - 1, n_trials);

        MFREE(active);
        return OK;
ERROR:
        if(active){
                MFREE(active);
        }
        return FAIL;
}

void recursive_aln_inline(struct msa* msa, struct aln_tasks* t,
                           struct aln_param* ap, uint8_t* active, int c,
                           int n_trials)
{
        struct task* local_t = NULL;
        int a, b;

        local_t = t->list[c];
        a = local_t->a - msa->numseq;
        b = local_t->b - msa->numseq;

        if(!active[local_t->a] && local_t->a >= msa->numseq){
                recursive_aln_inline(msa, t, ap, active, a, n_trials);
        }
        if(!active[local_t->b] && local_t->b >= msa->numseq){
                recursive_aln_inline(msa, t, ap, active, b, n_trials);
        }

        struct aln_mem* ml = NULL;
        alloc_aln_mem(&ml, 256);
        ml->ap = ap;
        ml->mode = ALN_MODE_FULL;

        do_align_inline_refine(msa, t, ml, c, n_trials);

        active[local_t->a] = 0;
        active[local_t->b] = 0;
        active[local_t->c] = 1;

        free_aln_mem(ml);
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

        RUN(init_alnmem(m));

        /* Compute consistency bonus for all merge types */
        m->consistency = NULL;
        m->consistency_stride = 0;
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
                                &m->consistency));
                        m->consistency_stride = dp_cols;
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
                MFREE(m->consistency);
                m->consistency = NULL;
                m->consistency_stride = 0;
        }

        /* Store confidence */
        t->list[task_id]->confidence = best_sp;

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
        return FAIL;
}

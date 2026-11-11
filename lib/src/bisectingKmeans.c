#include "tldevel.h"
#include "threadpool.h"

#ifdef HAVE_AVX2
#include <xmmintrin.h>
#include <mm_malloc.h>
#endif


#include "tlrng.h"

#include "bisectingKmeans.h"

#include "msa_struct.h"
/* #include "global.h" */
#include "task.h"
#include "sequence_distance.h"
#include "guide_dist.h"
#include "euclidean_dist.h"

/* #include "alignment.h" */
#include "pick_anchor.h"
#include "esl_stopwatch.h"

#include "sl_tree.h"   /* struct node; single-linkage and hybrid guide trees */

/* thread pool for the single-linkage trees (NULL: serial) */
#define SL_POOL(msa) ((void*)(msa)->pool)

struct kmeans_result{
        int* sl;
        int* sr;
        int nl;
        int nr;
        float score;
};

static struct kmeans_result* alloc_kmeans_result(int num_samples);
static void free_kmeans_results(struct kmeans_result* k);

struct node* upgma(float **dm,int* samples, int numseq);
static struct node* alloc_node(void);

static int label_internal(struct node *n, int label);
static void create_tasks(struct node*n, struct aln_tasks* t);
/* static void create_tasks(struct node*n, struct aln_tasks* t); */


/* static int bisecting_kmeans_serial(struct msa *msa, struct node **ret_n, float **dm, int *samples, int num_samples); */
static int bisecting_kmeans(struct msa* msa, struct node** ret_n,
                            const float * const * dm,
                            int* samples, int num_samples,
                            pair_dist_fn leaf_dist, void* leaf_ctx);
/* static int bisecting_kmeans_parallel(struct msa* msa, struct node** ret_n, float** dm,int* samples, int num_samples); */

/* One pending bisection: the samples to split and where to hang the result. */
struct bjob {
        struct node** ret_n;
        int* samples;
        int num_samples;
};
static int bisect_step(struct msa* msa, struct node** ret_n, const float * const * dm,
                       int* samples, int num_samples, pair_dist_fn leaf_dist, void* leaf_ctx,
                       struct bjob* out, int* nout);

static int split2(const float * const * dm,const int* samples, const int num_anchors,const int num_samples,const int seed_pick,struct kmeans_result** ret);

/* Wrapper args for threadpool task submission */
struct split2_task_arg {
        const float *const *dm;
        const int *samples;
        int num_anchors;
        int num_samples;
        int seed_idx;
        struct kmeans_result **res;
};

static void split2_task_fn(void *arg)
{
        struct split2_task_arg *a = (struct split2_task_arg *)arg;
        split2(a->dm, a->samples, a->num_anchors, a->num_samples, a->seed_idx, a->res);
}

static inline int cmp_floats(const float a, const float b);

inline int cmp_floats(const float a, const float b)
{
        const float epsilon = 1e-6; // Set a small epsilon value for tolerance
        if (fabsf(a - b) < epsilon) {
                return 0; // Numbers are equal
        } else if (a > b) {
                return 1; // First number is bigger
        } else {
                return -1; // Second number is bigger
        }
}

/* Default leaf-cluster distance: BPM on raw internal sequences */
static float** bpm_pair_dist(struct msa* msa, int* samples, int n, void* ctx)
{
        (void)ctx;
        return d_estimation(msa, samples, n, 1);
}

/* UPGMA leaf distances on the run's tree distance (tree_kmeans_dist = 1) */
static float** sl_pair_dist(struct msa* msa, int* samples, int n, void* ctx)
{
        (void)ctx;
        return sl_pair_matrix(msa, samples, n);
}

/* Anchor distances used by bisecting k-means when a pairwise-alignment guide
   distance is requested.  The top-level k-means is only reached when the
   input has at least KALIGN_KMEANS_UPGMA_THRESHOLD (floored at 3) sequences;
   below that the whole tree comes from the leaf distance, so the N x K
   alignments are skipped. */
static int guide_needs_anchor_dm(int numseq)
{
        int threshold = KALIGN_KMEANS_UPGMA_THRESHOLD;
        if(threshold < 3) threshold = 3;
        return numseq >= threshold;
}

static void free_anchor_dm(float** dm, int numseq)
{
        int i;
        if(!dm){
                return;
        }
        for(i = 0; i < numseq; i++){
                if(dm[i]){
#ifdef HAVE_AVX2
                        _mm_free(dm[i]);
#else
                        MFREE(dm[i]);
#endif
                }
        }
        MFREE(dm);
}

int build_tree_kmeans_noisy(struct msa* msa, struct aln_tasks** tasks,
                            uint64_t seed, float noise_sigma,
                            struct guide_dist_ctx* gctx)
{ KPERF_PHASE_REGION(KP_TREE, KPH_TREE);
        struct aln_tasks* t = NULL;
        struct node* root = NULL;
        float** dm = NULL;
        float** dm_tree = NULL;
        int* samples = NULL;
        int* anchors = NULL;
        int num_anchors;
        int numseq;
        int i;

        ASSERT(msa != NULL, "No alignment.");

        t = *tasks;
        if(!t){
                RUN(alloc_tasks(&t, msa->numseq));
        }
        numseq = msa->numseq;

        DECLARE_TIMER(timer);
        if(!msa->quiet){
                LOG_MSG("Calculating pairwise distances (noisy, seed=%lu)", (unsigned long)seed);
        }
        START_TIMER(timer);
        {
                float** anchor_bpm = NULL;
                RUNP(anchors = pick_anchor_dist(msa, &num_anchors, &anchor_bpm));
                dm = d_estimation_cached(msa, anchors, num_anchors, 0, anchor_bpm);
                pick_anchor_dist_free(anchor_bpm, numseq);
                if(!dm) ERROR_MSG("d_estimation failed");
        }
        if(gctx && guide_needs_anchor_dm(numseq)){
                RUNP(dm_tree = guide_dist_anchor_matrix(msa, gctx, anchors, num_anchors));
        }
        if(sl_kmeans_dist_override(msa)){
                free_anchor_dm(dm_tree, numseq);
                RUNP(dm_tree = sl_anchor_matrix(msa, anchors, num_anchors));
        }

        /* Add multiplicative Gaussian noise to distance matrix.
           The same factor perturbs the pairwise-alignment anchor distances
           (if used), so seq_distances stay identical to guide_dist = 0. */
        if(seed != 0 && noise_sigma > 0.0f){
                struct rng_state* rng = NULL;
                rng = init_rng(seed);
                for(i = 0; i < numseq; i++){
                        for(int j = 0; j < num_anchors; j++){
                                double noise = tl_random_gaussian(rng, 1.0, (double)noise_sigma);
                                if(noise < 0.1) noise = 0.1;
                                dm[i][j] *= (float)noise;
                                if(dm_tree){
                                        dm_tree[i][j] *= (float)noise;
                                }
                        }
                }
                free_rng(rng);
        }

        STOP_TIMER(timer);
        if(!msa->quiet){
                GET_TIMING(timer);
        }
        MFREE(anchors);

        MMALLOC(samples, sizeof(int) * numseq);
        for(i = 0; i < numseq; i++){
                samples[i] = i;
        }

        START_TIMER(timer);
        if(!msa->quiet){
                LOG_MSG("Building guide tree.");
        }

        if(sl_tree_mode(msa) == SL_MODE_PRIM || sl_tree_mode(msa) == SL_MODE_SLINK){
                RUN(sl_build(msa, sl_tree_mode(msa), noise_sigma, seed, SL_POOL(msa), &root));
        }else if(gctx){
                bisecting_kmeans(msa, &root,
                                 (const float * const *)(dm_tree ? dm_tree : dm),
                                 samples, numseq, sl_kmeans_dist_override(msa) ? sl_pair_dist : guide_dist_pair_matrix, gctx);
        }else{
                bisecting_kmeans(msa, &root, (const float * const *)(dm_tree ? dm_tree : dm), samples, numseq,
                                 sl_kmeans_dist_override(msa) ? sl_pair_dist : bpm_pair_dist, NULL);
        }
        ASSERT(root != NULL, "Guide tree construction failed.");

        STOP_TIMER(timer);
        if(!msa->quiet){
                GET_TIMING(timer);
        }
        label_internal(root, numseq);
        create_tasks(root, t);
        free_anchor_dm(dm_tree, numseq);
        dm_tree = NULL;

        /* Compute per-sequence normalized mean distance (noisy variant). */
        if(msa->seq_distances == NULL){
                MMALLOC(msa->seq_distances, sizeof(float) * numseq);
        }
        for(i = 0; i < numseq; i++){
                float sum = 0.0f;
                int j;
                for(j = 0; j < num_anchors; j++){
                        sum += dm[i][j];
                }
                float mean_dist = sum / (float)num_anchors;
                float seq_len = (float)msa->sequences[i]->len;
                msa->seq_distances[i] = (seq_len > 0.0f) ? mean_dist / seq_len : 0.0f;
        }

        MFREE(root);
        for(i = 0; i < numseq; i++){
#ifdef HAVE_AVX2
                _mm_free(dm[i]);
#else
                MFREE(dm[i]);
#endif
        }
        MFREE(dm);
        DESTROY_TIMER(timer);
        return OK;
ERROR:
        return FAIL;
}

int build_tree_kmeans(struct msa* msa, struct aln_tasks** tasks,
                      struct guide_dist_ctx* gctx)
{ KPERF_PHASE_REGION(KP_TREE, KPH_TREE);
        struct aln_tasks* t = NULL;
        struct node* root = NULL;
        float** dm = NULL;
        float** dm_tree = NULL;
        int* samples = NULL;
        int* anchors = NULL;
        int num_anchors;
        int numseq;

        int i;

        ASSERT(msa != NULL, "No alignment.");

        t = *tasks;
        if(!t){
                RUN(alloc_tasks(&t, msa->numseq));
        }
        numseq = msa->numseq;

        DECLARE_TIMER(timer);
        /* pick anchors . */
        if(!msa->quiet){
                LOG_MSG("Calculating pairwise distances");
        }
        START_TIMER(timer);
        {
                float** anchor_bpm = NULL;
                RUNP(anchors = pick_anchor_dist(msa, &num_anchors, &anchor_bpm));
                dm = d_estimation_cached(msa, anchors, num_anchors, 0, anchor_bpm);
                pick_anchor_dist_free(anchor_bpm, numseq);
                if(!dm) ERROR_MSG("d_estimation failed");
        }
        /* guide_dist = 1: tree from pairwise-alignment distances; the BPM
           matrix above is still used for seq_distances (VSM, gap scaling,
           consistency anchors), so only the guide tree changes. */
        if(gctx && guide_needs_anchor_dm(numseq)){
                RUNP(dm_tree = guide_dist_anchor_matrix(msa, gctx, anchors, num_anchors));
        }
        if(sl_kmeans_dist_override(msa)){
                free_anchor_dm(dm_tree, numseq);
                RUNP(dm_tree = sl_anchor_matrix(msa, anchors, num_anchors));
        }


        STOP_TIMER(timer);
        if(!msa->quiet){
                GET_TIMING(timer);
        }
        MFREE(anchors);

        MMALLOC(samples, sizeof(int)* numseq);
        for(i = 0; i < numseq;i++){
                samples[i] = i;
        }

        START_TIMER(timer);
        if(!msa->quiet){
                LOG_MSG("Building guide tree.");
        }

        /* if(n_threads == 1){ */
        /*         RUN(bisecting_kmeans_serial(msa,&root, dm, samples, numseq)); */
        /* }else{ */
        if(sl_tree_mode(msa) == SL_MODE_PRIM || sl_tree_mode(msa) == SL_MODE_SLINK){
                RUN(sl_build(msa, sl_tree_mode(msa), 0.0f, 0, SL_POOL(msa), &root));
        }else if(gctx){
                bisecting_kmeans(msa, &root,
                                 (const float * const *)(dm_tree ? dm_tree : dm),
                                 samples, numseq, sl_kmeans_dist_override(msa) ? sl_pair_dist : guide_dist_pair_matrix, gctx);
        }else{
                bisecting_kmeans(msa,&root, (const float * const *)(dm_tree ? dm_tree : dm), samples, numseq,
                                 sl_kmeans_dist_override(msa) ? sl_pair_dist : bpm_pair_dist, NULL);
        }
        /* } */
        ASSERT(root != NULL, "Guide tree construction failed.");

        STOP_TIMER(timer);
        if(!msa->quiet){
                GET_TIMING(timer);
        }
        label_internal(root, numseq);

        create_tasks(root, t);
        free_anchor_dm(dm_tree, numseq);
        dm_tree = NULL;

        /* Compute per-sequence normalized mean distance to anchors.
           Used for distance-dependent gap penalty scaling. */
        if(msa->seq_distances == NULL){
                MMALLOC(msa->seq_distances, sizeof(float) * numseq);
        }
        for(i = 0; i < numseq; i++){
                float sum = 0.0f;
                int j;
                for(j = 0; j < num_anchors; j++){
                        sum += dm[i][j];
                }
                float mean_dist = sum / (float)num_anchors;
                float seq_len = (float)msa->sequences[i]->len;
                msa->seq_distances[i] = (seq_len > 0.0f) ? mean_dist / seq_len : 0.0f;
        }

        MFREE(root);
        for(i = 0; i < numseq; i++){
#ifdef HAVE_AVX2
                _mm_free(dm[i]);
#else
                MFREE(dm[i]);
#endif
        }
        MFREE(dm);
        DESTROY_TIMER(timer);
        return OK;
ERROR:
        return FAIL;
}

/* Non-recursive driver for bisect_step.  Pending jobs live on a heap stack,
   so C stack use does not depend on the depth of the tree.  With a thread
   pool the left half of each split is handed to the pool (all tasks join one
   group owned by the caller; no task waits on another) and the right half is
   continued on the same thread.  The tree built does not depend on the
   order in which jobs run. */
struct bis_ctx {
        struct msa* msa;
        const float * const * dm;
        pair_dist_fn leaf_dist;
        void* leaf_ctx;
        tp_group_t* g;
};

struct bis_task {
        struct bis_ctx* ctx;
        struct bjob job;
};

static int bis_run(struct bis_ctx* c, struct bjob first, int* first_status);

static void bis_task_fn(void* arg)
{
        struct bis_task* bt = (struct bis_task*)arg;
        struct bis_ctx* c = bt->ctx;
        struct bjob j = bt->job;
        free(bt);
        bis_run(c, j, NULL);
}

static int bis_run(struct bis_ctx* c, struct bjob first, int* first_status)
{
        struct bjob* stack = NULL;
        int cap = 64;
        int top = 0;
        int ret = OK;

        stack = malloc(sizeof(struct bjob) * (size_t)cap);
        if(!stack){
                ERROR_MSG("out of memory");
        }
        stack[top++] = first;
        while(top > 0){
                struct bjob j = stack[--top];
                struct bjob out[2];
                int nout = 0;
                int r = bisect_step(c->msa, j.ret_n, c->dm, j.samples, j.num_samples,
                                    c->leaf_dist, c->leaf_ctx, out, &nout);
                if(first_status){
                        *first_status = r;
                        first_status = NULL;
                }
                if(r != OK){
                        ret = FAIL;
                }
                if(nout == 2){
                        int k = 1;
                        if(c->g){
                                struct bis_task* bt = malloc(sizeof(struct bis_task));
                                if(bt){
                                        bt->ctx = c;
                                        bt->job = out[0];
                                        if(tp_group_submit(c->g, bis_task_fn, bt) == 0){
                                                k = 1;
                                                /* left submitted: continue with the right half only */
                                                out[0] = out[1];
                                                nout = 1;
                                        }else{
                                                free(bt);
                                        }
                                }
                        }
                        if(top + nout > cap){
                                struct bjob* ns;
                                cap *= 2;
                                ns = realloc(stack, sizeof(struct bjob) * (size_t)cap);
                                if(!ns){
                                        ERROR_MSG("out of memory");
                                }
                                stack = ns;
                        }
                        /* push right first so that the left half is processed first */
                        for(k = nout - 1; k >= 0; k--){
                                stack[top++] = out[k];
                        }
                }
        }
        free(stack);
        return ret;
ERROR:
        if(stack){
                free(stack);
        }
        return FAIL;
}

int bisecting_kmeans(struct msa* msa, struct node** ret_n, const float * const * dm,
                     int* samples, int num_samples, pair_dist_fn leaf_dist, void* leaf_ctx)
{
        struct bis_ctx c;
        struct bjob root;
        int status = OK;

        c.msa = msa;
        c.dm = dm;
        c.leaf_dist = leaf_dist;
        c.leaf_ctx = leaf_ctx;
        root.ret_n = ret_n;
        root.samples = samples;
        root.num_samples = num_samples;
        c.g = NULL;
        if(msa->pool){
                c.g = tp_group_create(msa->pool);
        }
        bis_run(&c, root, &status);
        if(c.g){
                tp_group_wait(c.g);
                tp_group_destroy(c.g);
        }
        return status;
}

int bisect_step(struct msa* msa, struct node** ret_n, const float * const * dm,
                int* samples, int num_samples, pair_dist_fn leaf_dist, void* leaf_ctx,
                struct bjob* out, int* nout)
{ KPERF_REGION(KP_TREE_KMEANS);
        struct kmeans_result* res_tmp = NULL;
        struct kmeans_result* best = NULL;
        /* struct kmeans_result** res = NULL; */
        struct node* n = NULL;
        int num_anchors = 0;

        int i,j;
        int tries = 40;
        /* int t_iter; */
        /* int r; */
        int* sl = NULL;
        int* sr = NULL;
        int num_l,num_r;

        /* LOG_MSG("num_samples: %d", num_samples); */

        *nout = 0;

        /* Base cases: 0 or 1 samples cannot be split further */
        if(num_samples <= 0){
                return OK;
        }
        if(num_samples == 1){
                n = alloc_node();
                n->id = samples[0];
                *ret_n = n;
                MFREE(samples);
                return OK;
        }

        num_anchors = MACRO_MIN(32, msa->numseq);

        /* K-means needs at least 3 samples for a non-degenerate bisection.
           Floor the threshold to 3 so UPGMA handles tiny inputs safely. */
        {
                int threshold = sl_kmeans_leaf(msa);
                if(threshold < 3) threshold = 3;
                if(num_samples < threshold && sl_kmeans_leaf_is_sl(msa)){
                        /* hybrid tree: the group's subtree by single linkage */
                        float nz = msa->tree.in_realign ? 0.0f : msa->tree.sl_noise;
                        RUN(sl_build_subset(msa, samples, num_samples, SL_MODE_PRIM, nz,
                                            nz > 0.0f ? msa->tree.sl_seed : 0, SL_POOL(msa), &n));
                        *ret_n = n;
                        MFREE(samples);
                        return OK;
                }
                if(num_samples < threshold){
                        float** dm_local = NULL;
                        RUNP(dm_local = leaf_dist(msa, samples, num_samples, leaf_ctx));
                        n = upgma(dm_local, samples, num_samples);
                        *ret_n = n;
                        gfree(dm_local);
                        MFREE(samples);
                        return OK;
                }
        }

        /* else if(num_samples < 1000){ */
        /*         RUN(bisecting_kmeans_serial(msa, &n, dm, samples, num_samples)); */
        /*         *ret_n = n; */
        /*         return OK; */
        /* } */


        best = NULL;
        res_tmp = NULL;
        struct kmeans_result* res[4];

        /* MMALLOC(res, sizeof(struct kmeans_result*) * 4); */
        for(i = 0; i < 4;i++){
                res[i] = NULL;
        }
        tries = MACRO_MIN(tries, num_samples);
        int step = num_samples / tries;
        int change = 0;
        for(i = 0;i < tries;i += 4){
                change = 0;

                {
                        struct split2_task_arg args[4];
                        for (int t = 0; t < 4; t++) {
                                args[t].dm = dm;
                                args[t].samples = samples;
                                args[t].num_anchors = num_anchors;
                                args[t].num_samples = num_samples;
                                args[t].seed_idx = ((i + t) * step) % num_samples;
                                args[t].res = &res[t];
                        }
                        tp_group_t *g = tp_group_create(msa->pool);
                        for (int t = 0; t < 4; t++)
                                tp_group_submit(g, split2_task_fn, &args[t]);
                        tp_group_wait(g);
                        tp_group_destroy(g);
                }

                for(j = 0; j < 4;j++){
                        if(!best){
                                change++;
                                best = res[j];
                                res[j] = NULL;
                        }else{
                                if(best->score > res[j]->score){
                                        res_tmp = best;
                                        best = res[j];
                                        res[j] = res_tmp;
                                        /* LOG_MSG("Better!!! %f %f", res_tmp->score,best->score); */

                                        change++;
                                }
                        }
                }
                if(!change){
                        break;
                }
        }

        sl = best->sl;
        sr = best->sr;

        num_l = best->nl;
        num_r = best->nr;

        /* free_kmeans_results(res[0]); */
        /* free_kmeans_results(res[1]); */
        /* free_kmeans_results(res[2]); */
        /* free_kmeans_results(res[3]); */

        for(i = 0; i < 4;i++){
                free_kmeans_results(res[i]);
        }
        /* MFREE(res); */
        MFREE(best);

        MFREE(samples);
        n = alloc_node();

        /* The two halves are bisected by the caller (bisecting_kmeans),
           so this function does not recurse. */
        out[0].ret_n = &n->left;
        out[0].samples = sl;
        out[0].num_samples = num_l;
        out[1].ret_n = &n->right;
        out[1].samples = sr;
        out[1].num_samples = num_r;
        *nout = 2;
        *ret_n = n;
        return OK;
ERROR:
        return FAIL;
}

/* int bisecting_kmeans_serial(struct msa* msa, struct node** ret_n, float** dm,int* samples, int num_samples) */
/* { */
/*         struct kmeans_result* res_tmp = NULL; */
/*         struct kmeans_result* best = NULL; */
/*         struct kmeans_result** res_ptr = NULL; */
/*         int num_anchors = 0; */
/*         struct node* n = NULL; */
/*         int i,j; */
/*         int tries = 40; */
/*         /\* int t_iter; *\/ */
/*         /\* int r; *\/ */
/*         int* sl = NULL; */
/*         int* sr = NULL; */
/*         int num_l,num_r; */

/*         num_anchors = MACRO_MIN(32, msa->numseq); */

/*         if(num_samples < msa->kmeans_upgma_threshold){ */
/*                 float** dm = NULL; */
/*                 RUNP(dm = d_estimation(msa, samples, num_samples,1));// anchors, num_anchors,1)); */
/*                 n = upgma(dm,samples, num_samples); */
/*                 *ret_n = n; */
/*                 gfree(dm); */
/*                 MFREE(samples); */
/*                 return OK; */
/*         } */

/*         best = NULL; */
/*         res_tmp = NULL; */

/*         MMALLOC(res_ptr, sizeof(struct kmeans_result*) * 4); */
/*         for(i = 0; i < 4;i++){ */
/*                 res_ptr[i] = NULL; */
/*         } */
/*         tries = MACRO_MIN(tries, num_samples); */
/*         int step = num_samples / tries; */
/*         int change = 0; */
/*         for(i = 0;i < tries;i += 4){ */
/*                 change = 0; */
/*                 for(j = 0; j < 4;j++){ */
/*                         split(dm,samples,num_anchors, num_samples, (i+ j)*step, &res_ptr[j]); */
/*                 } */

/*                 for(j = 0; j < 4;j++){ */
/*                         if(!best){ */
/*                                 change++; */
/*                                 best = res_ptr[j]; */
/*                                 res_ptr[j] = NULL; */
/*                         }else{ */
/*                                 if(best->score > res_ptr[j]->score){ */
/*                                         res_tmp = best; */
/*                                         best = res_ptr[j]; */
/*                                         res_ptr[j] = res_tmp; */
/*                                         /\* LOG_MSG("Better!!! %f %f", res_tmp->score,best->score); *\/ */

/*                                         change++; */
/*                                 } */
/*                         } */
/*                 } */
/*                 if(!change){ */
/*                         break; */
/*                 } */

/*         } */
/*         sl = best->sl; */
/*         sr = best->sr; */

/*         num_l = best->nl; */
/*         num_r = best->nr; */

/*         for(i = 0; i < 4;i++){ */
/*                 free_kmeans_results(res_ptr[i]); */
/*         } */
/*         MFREE(res_ptr); */
/*         MFREE(best); */

/*         MFREE(samples); */
/*         n = alloc_node(); */

/*         bisecting_kmeans_serial(msa,&n->left , dm, sl, num_l); */
/*         bisecting_kmeans_serial(msa,&n->right, dm, sr, num_r); */
/*         *ret_n = n; */
/*         return OK; */
/* ERROR: */
/*         return FAIL; */
/* } */

int split2(const float * const * dm,const int* samples, const int num_anchors,const int num_samples,const int seed_pick,struct kmeans_result** ret)
{
        struct kmeans_result* res = NULL;
        int* sl = NULL;
        int* sr = NULL;
        int num_l,num_r;
        float* w = NULL;
        float* wl = NULL;
        float* wr = NULL;
        float* cl = NULL;
        float* cr = NULL;
        float dl = 0.0F;
        float dr = 0.0F;
        float score;
        int num_var;
        int i;
        int s;
        int j;

        num_var = num_anchors / 8;
        if( num_anchors%8){
                num_var++;
        }
        num_var = num_var << 3;




#ifdef HAVE_AVX2
        wr = _mm_malloc(sizeof(float) * num_var,32);
        wl = _mm_malloc(sizeof(float) * num_var,32);
        cr = _mm_malloc(sizeof(float) * num_var,32);
        cl = _mm_malloc(sizeof(float) * num_var,32);
        w = _mm_malloc(sizeof(float) * num_var,32);
#else
        MMALLOC(wr,sizeof(float) * num_var);
        MMALLOC(wl,sizeof(float) * num_var);
        MMALLOC(cr,sizeof(float) * num_var);
        MMALLOC(cl,sizeof(float) * num_var);
        MMALLOC(w,sizeof(float) * num_var);
#endif

        if(*ret){
                res = *ret;
        }else{
                RUNP(res = alloc_kmeans_result(num_samples));
        }

        res->score = FLT_MAX;

        sl = res->sl;
        sr = res->sr;


        for(i = 0; i < num_var;i++){
                w[i] = 0.0F;
                wr[i] = 0.0F;
                wl[i] = 0.0F;
                cr[i] = 0.0F;
                cl[i] = 0.0F;
        }
        for(i = 0; i < num_samples;i++){
                s = samples[i];
                for(j = 0; j < num_anchors;j++){
                        w[j] += dm[s][j];
                }
        }

        for(j = 0; j < num_anchors;j++){
                w[j] /= (float)num_samples;
        }
        //r = tl_random_int(rng  , num_samples);
        //r = sel[t_iter];

        s = samples[seed_pick];
        /* LOG_MSG("Selected %d\n",s); */
        for(j = 0; j < num_anchors;j++){
                cl[j] = dm[s][j];
        }

        for(j = 0; j < num_anchors;j++){
                cr[j] = w[j] - (cl[j] - w[j]);

        }

#ifdef HAVE_AVX2
        _mm_free(w);
#else
        MFREE(w);
#endif

        w = NULL;
        for(int stop = 0; stop < 500; stop++){
                num_l = 0;
                num_r = 0;

                for(i = 0; i < num_anchors;i++){
                        wr[i] = 0.0F;
                        wl[i] = 0.0F;
                }
                score = 0.0f;
                for(i = 0; i < num_samples;i++){
                        s = samples[i];
#ifdef HAVE_AVX2
                        edist_256(dm[s], cl, num_anchors, &dl);
                        edist_256(dm[s], cr, num_anchors, &dr);
#else
                        edist_serial(dm[s], cl, num_anchors, &dl);
                        edist_serial(dm[s], cr, num_anchors, &dr);
#endif
                        score += MACRO_MIN(dl,dr);


                        int res = cmp_floats(dr,dl);
                        if(res == -1){
                                w = wr;
                                sr[num_r] = s;
                                num_r++;
                        }else if (res == 1){
                                w = wl;
                                sl[num_l] = s;
                                num_l++;
                        }else{
                                if(i & 1){
                                        w = wr;
                                        sr[num_r] = s;
                                        num_r++;
                                }else{
                                        w = wl;
                                        sl[num_l] = s;
                                        num_l++;
                                }
                        }
                        for(j = 0; j < num_anchors;j++){
                                w[j] += dm[s][j];
                        }
                }
                if(num_l == 0 || num_r == 0){
                        score = 0.0F;
                        num_l = 0;
                        num_r = 0;

                        for(i = 0; i < num_samples/2;i++){
                                sl[num_l] = samples[i];
                                num_l++;
                        }
                        for(i = num_samples/2; i < num_samples;i++){
                                sr[num_r] = samples[i];
                                num_r++;
                        }
                        break;
                }

                for(j = 0; j < num_anchors;j++){
                        wl[j] /= (float)num_l;
                        wr[j] /= (float)num_r;
                }

                s = 0;

                for(j = 0; j < num_anchors;j++){
                        int res = cmp_floats(wl[j],cl[j]);
                        if(res != 0){
                                s = 1;
                                break;
                        }
                        res = cmp_floats(wr[j],cr[j]);
                        if(res != 0){
                                s = 1;
                                break;
                        }
                }

                if(s){
                        w = cl;
                        cl = wl;
                        wl = w;

                        w = cr;
                        cr = wr;
                        wr = w;
                }else{
                        break;
                }
        }

#ifdef HAVE_AVX2
        _mm_free(wr);
        _mm_free(wl);
        _mm_free(cr);
        _mm_free(cl);
#else
        MFREE(wr);
        MFREE(wl);
        MFREE(cr);
        MFREE(cl);
#endif

        res->nl =  num_l;
        res->nr =  num_r;
        res->score = score;
        *ret = res;
        return OK;
ERROR:
        return FAIL;
}


struct node* upgma(float **dm,int* samples, int numseq)
{
        struct node** tree = NULL;
        struct node* tmp = NULL;

        int i,j;
        int *as = NULL;

        float max;
        int node_a = 0;
        int node_b = 0;
        int cnode = numseq;
        int numprofiles;


        numprofiles = (numseq << 1) - 1;

        MMALLOC(as,sizeof(int)*numseq);
        for (i = numseq; i--;){
                as[i] = i+1;
        }


        MMALLOC(tree,sizeof(struct node*)*numseq);
        for (i = 0;i < numseq;i++){
                tree[i] = NULL;
                tree[i] = alloc_node();
                tree[i]->id = samples[i];
        }

        while (cnode != numprofiles){
                max = FLT_MAX;
                for (i = 0;i < numseq-1; i++){
                        if (as[i]){
                                for ( j = i + 1;j < numseq;j++){
                                        if (as[j]){
                                                if (dm[i][j] < max){
                                                        max = dm[i][j];
                                                        node_a = i;
                                                        node_b = j;
                                                }
                                        }
                                }
                        }
                }
                tmp = NULL;
                tmp = alloc_node();
                tmp->left = tree[node_a];
                tmp->right = tree[node_b];


                tree[node_a] = tmp;
                tree[node_b] = NULL;

                /*deactivate  sequences to be joined*/
                as[node_a] = cnode+1;
                as[node_b] = 0;


                cnode++;

                /*calculate new distances*/
                for (j = numseq;j--;){
                        if (j != node_b){
                                dm[node_a][j] = (dm[node_a][j] + dm[node_b][j])*0.5F + 0.001F;
                        }
                        //fprintf(stdout,"\n");
                }
                dm[node_a][node_a] = 0.0F;
                for (j = numseq;j--;){
                        dm[j][node_a] = dm[node_a][j];
                }
        }
        tmp = tree[node_a];
        MFREE(tree);
        MFREE(as);
        return tmp;
ERROR:
        return NULL;
}

struct node* alloc_node(void)
{
        struct node* n = NULL;
        MMALLOC(n, sizeof(struct node));
        n->left = NULL;
        n->right = NULL;
        n->id = -1;
        return n;
ERROR:
        return NULL;
}

struct tree_frame {
        struct node* n;
        int st;
};

static void tree_stack_grow(void** stack, int* cap, size_t elem)
{
        void* ns;
        *cap *= 2;
        ns = realloc(*stack, elem * (size_t)*cap);
        if(!ns){
                fprintf(stderr, "out of memory in tree walk\n");
                abort();
        }
        *stack = ns;
}

/* Post-order labelling (left, right, self), with an explicit stack. */
int label_internal(struct node*n, int label)
{
        struct tree_frame* stack = NULL;
        int cap = 64;
        int top = 0;

        stack = malloc(sizeof(struct tree_frame) * (size_t)cap);
        if(!stack){
                fprintf(stderr, "out of memory in tree walk\n");
                abort();
        }
        stack[top].n = n;
        stack[top].st = 0;
        top++;
        while(top > 0){
                struct tree_frame* f = &stack[top - 1];
                struct node* child = NULL;
                if(f->st == 0){
                        f->st = 1;
                        child = f->n->left;
                }else if(f->st == 1){
                        f->st = 2;
                        child = f->n->right;
                }else{
                        if(f->n->id == -1){
                                f->n->id = label;
                                label++;
                        }
                        top--;
                        continue;
                }
                if(child){
                        if(top == cap){
                                tree_stack_grow((void**)&stack, &cap, sizeof(struct tree_frame));
                        }
                        stack[top].n = child;
                        stack[top].st = 0;
                        top++;
                }
        }
        free(stack);
        return label;
}

/* Pre-order task emission, with an explicit stack.  As before, the two
   children of a node that has both are freed once visited. */
void create_tasks(struct node*n, struct aln_tasks* t)
{
        struct tree_frame* stack = NULL;   /* st: 1 = free node after visiting */
        int cap = 64;
        int top = 0;

        stack = malloc(sizeof(struct tree_frame) * (size_t)cap);
        if(!stack){
                fprintf(stderr, "out of memory in tree walk\n");
                abort();
        }
        stack[top].n = n;
        stack[top].st = 0;
        top++;
        while(top > 0){
                struct node* x = stack[--top].n;
                int free_it = stack[top].st;
                struct node* l = x->left;
                struct node* r = x->right;
                if(l && r){
                        struct task* task;
                        task = t->list[t->n_tasks];
                        task->a = l->id;
                        task->b = r->id;
                        task->c = x->id;
                        t->n_tasks++;
                }
                if(top + 2 > cap){
                        tree_stack_grow((void**)&stack, &cap, sizeof(struct tree_frame));
                }
                /* push right first so that the left subtree is emitted first */
                if(r){
                        stack[top].n = r;
                        stack[top].st = (l != NULL);
                        top++;
                }
                if(l){
                        stack[top].n = l;
                        stack[top].st = (r != NULL);
                        top++;
                }
                if(free_it){
                        MFREE(x);
                }
        }
        free(stack);
}


struct kmeans_result* alloc_kmeans_result(int num_samples)
{
        struct kmeans_result* k = NULL;
        ASSERT(num_samples != 0, "No samples???");

        MMALLOC(k, sizeof(struct kmeans_result));

        k->nl = 0;
        k->nr = 0;
        k->sl = NULL;
        k->sr = NULL;
        MMALLOC(k->sl, sizeof(int) * num_samples);
        MMALLOC(k->sr, sizeof(int) * num_samples);
        k->score = FLT_MAX;
        return k;
ERROR:
        free_kmeans_results(k);
        return NULL;
}

void free_kmeans_results(struct kmeans_result* k)
{
        if(k){
                if(k->sl){
                        MFREE(k->sl);
                }
                if(k->sr){
                        MFREE(k->sr);
                }
                MFREE(k);
        }
}

/* Build guide tree from a pre-computed N×K distance matrix using bisecting
   k-means.  Used by the realign loop where distances are computed from aligned
   sequences (identity distances, already 0..1) rather than BPM on raw sequences.

   The caller owns dm and is responsible for freeing it.  This function does NOT
   free dm — unlike build_tree_kmeans which computes and frees its own dm. */
int build_tree_kmeans_from_dm(struct msa* msa, struct aln_tasks** tasks,
                               float** dm, int num_anchors,
                               pair_dist_fn leaf_dist)
{ KPERF_PHASE_REGION(KP_TREE, KPH_TREE);
        struct aln_tasks* t = NULL;
        struct node* root = NULL;
        int* samples = NULL;
        int numseq;
        int i;

        ASSERT(msa != NULL, "No alignment.");
        ASSERT(dm != NULL, "No distance matrix.");
        ASSERT(num_anchors > 0, "num_anchors must be > 0");

        t = *tasks;
        if(!t){
                RUN(alloc_tasks(&t, msa->numseq));
        }
        numseq = msa->numseq;

        MMALLOC(samples, sizeof(int) * numseq);
        for(i = 0; i < numseq; i++){
                samples[i] = i;
        }

        /* Bisecting k-means — same parallel algorithm as initial tree.
           Note: bisecting_kmeans() takes ownership of samples and frees it. */
        bisecting_kmeans(msa, &root, (const float * const *)dm, samples, numseq, leaf_dist, NULL);
        samples = NULL;  /* owned and freed by bisecting_kmeans */
        ASSERT(root != NULL, "Bisecting k-means tree construction failed.");

        label_internal(root, numseq);
        create_tasks(root, t);

        /* Compute per-sequence mean distance to anchors.
           Identity distances are already 0..1, so NO length normalization
           (unlike BPM distances in build_tree_kmeans which divides by seq_len). */
        if(msa->seq_distances == NULL){
                MMALLOC(msa->seq_distances, sizeof(float) * numseq);
        }
        for(i = 0; i < numseq; i++){
                float sum = 0.0f;
                int j;
                for(j = 0; j < num_anchors; j++){
                        sum += dm[i][j];
                }
                msa->seq_distances[i] = sum / (float)num_anchors;
        }

        MFREE(root);
        *tasks = t;
        return OK;
ERROR:
        if(samples) MFREE(samples);
        return FAIL;
}

int build_tree_from_pairwise(struct msa* msa, struct aln_tasks** tasks, float** dm)
{ KPERF_PHASE_REGION(KP_TREE, KPH_TREE);
        struct aln_tasks* t = NULL;
        struct node* root = NULL;
        int* samples = NULL;
        int numseq;
        int i, j;

        ASSERT(msa != NULL, "No alignment.");
        ASSERT(dm != NULL, "No distance matrix.");

        t = *tasks;
        if(!t){
                RUN(alloc_tasks(&t, msa->numseq));
        }
        numseq = msa->numseq;

        /* Compute per-sequence mean pairwise distance (for gap/VSM scaling).
           Must be done BEFORE upgma() which modifies dm in place. */
        if(msa->seq_distances == NULL){
                MMALLOC(msa->seq_distances, sizeof(float) * numseq);
        }
        for(i = 0; i < numseq; i++){
                float sum = 0.0f;
                for(j = 0; j < numseq; j++){
                        if(j != i) sum += dm[i][j];
                }
                msa->seq_distances[i] = (numseq > 1) ? sum / (float)(numseq - 1) : 0.0f;
        }

        MMALLOC(samples, sizeof(int) * numseq);
        for(i = 0; i < numseq; i++){
                samples[i] = i;
        }

        /* Build UPGMA tree from pairwise distances.
           Note: upgma() modifies dm in-place. */
        root = upgma(dm, samples, numseq);
        ASSERT(root != NULL, "UPGMA tree construction failed.");

        label_internal(root, numseq);
        create_tasks(root, t);

        MFREE(samples);
        MFREE(root);
        *tasks = t;
        return OK;
ERROR:
        if(samples) MFREE(samples);
        return FAIL;
}

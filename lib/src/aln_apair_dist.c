#include "tldevel.h"
#include "msa_struct.h"

#ifdef USE_THREADPOOL
#include "threadpool.h"
#endif

#define ALN_APAIR_DIST_IMPORT
#include "aln_apair_dist.h"

static float pairwise_identity_dist(const char* a, const char* b, int alnlen);

#ifdef USE_THREADPOOL
struct apair_ctx {
        struct msa_seq** seqs;
        float** dm;
        int n;
        int alnlen;
};

static void apair_row_fn(int row_start, int row_end, void *arg)
{
        struct apair_ctx *c = (struct apair_ctx *)arg;
        for (int i = row_start; i < row_end; i++) {
                const char *seq_i = c->seqs[i]->seq;
                for (int j = i + 1; j < c->n; j++) {
                        float d = pairwise_identity_dist(seq_i, c->seqs[j]->seq, c->alnlen);
                        c->dm[i][j] = d;
                        c->dm[j][i] = d;
                }
        }
}
#endif

int compute_aln_pairwise_dist(struct msa* msa, float*** dm_ptr)
{
        float** dm = NULL;
        int n;
        int i, j;

        ASSERT(msa != NULL, "No MSA");
        ASSERT(msa->aligned == ALN_STATUS_FINAL, "MSA must be finalized");

        n = msa->numseq;

        MMALLOC(dm, sizeof(float*) * n);
        for(i = 0; i < n; i++){
                dm[i] = NULL;
        }
        for(i = 0; i < n; i++){
                MMALLOC(dm[i], sizeof(float) * n);
                dm[i][i] = 0.0f;
        }

#ifdef USE_THREADPOOL
        if(msa->pool && n >= KALIGN_DIST_MIN_SEQS){
                struct apair_ctx ctx = { msa->sequences, dm, n, msa->alnlen };
                tp_parallel_for_chunked(msa->pool, 0, n - 1, KALIGN_PFOR_MIN_CHUNK, apair_row_fn, &ctx);
        }else{
#endif
        for(i = 0; i < n - 1; i++){
                const char* seq_i = msa->sequences[i]->seq;
                for(j = i + 1; j < n; j++){
                        float d = pairwise_identity_dist(seq_i, msa->sequences[j]->seq,
                                                         msa->alnlen);
                        dm[i][j] = d;
                        dm[j][i] = d;
                }
        }
#ifdef USE_THREADPOOL
        }
#endif

        *dm_ptr = dm;
        return OK;
ERROR:
        if(dm){
                for(i = 0; i < n; i++){
                        if(dm[i]) MFREE(dm[i]);
                }
                MFREE(dm);
        }
        return FAIL;
}

void free_aln_dm(float** dm, int n)
{
        int i;
        if(dm == NULL) return;
        for(i = 0; i < n; i++){
                if(dm[i]) MFREE(dm[i]);
        }
        MFREE(dm);
}

/* Distance = 1.0 - identity.
   Only counts columns where both sequences have a residue (no gap). */
float pairwise_identity_dist(const char* a, const char* b, int alnlen)
{
        int matches = 0;
        int aligned = 0;
        int i;

        for(i = 0; i < alnlen; i++){
                if(a[i] != '-' && b[i] != '-'){
                        aligned++;
                        if(a[i] == b[i]){
                                matches++;
                        }
                }
        }

        if(aligned == 0){
                return 1.0f;
        }
        return 1.0f - (float)matches / (float)aligned;
}

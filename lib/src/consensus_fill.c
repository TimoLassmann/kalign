#include "tldevel.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "msa_struct.h"
#include "msa_op.h"
#include "aln_wrap.h"

#include "kalign/kalign.h"
#include "kalign/kalign_config.h"

#define CONSENSUS_FILL_IMPORT
#include "consensus_fill.h"

/* Realign the residue strings frag[0..n-1] (lengths flen[]) and write
   the aligned block (width *out_w) to blk[i] (caller frees).  Empty
   fragments come back as all-gap rows. */
static int align_segment(char** frag, int* flen, int n,
                         const struct msa* parent,
                         const struct kalign_run_config* cfg,
                         int n_threads,
                         char** blk, int* out_w)
{
        struct msa* sub = NULL;
        char** in = NULL;
        int* in_len = NULL;
        int* map = NULL;
        int k = 0;
        int max_len = 0;
        int w;
        int i;

        MMALLOC(in, sizeof(char*) * n);
        MMALLOC(in_len, sizeof(int) * n);
        MMALLOC(map, sizeof(int) * n);
        for(i = 0; i < n; i++){
                blk[i] = NULL;
                if(flen[i] > 0){
                        in[k] = frag[i];
                        in_len[k] = flen[i];
                        map[k] = i;
                        k++;
                        if(flen[i] > max_len) max_len = flen[i];
                }
        }

        if(k <= 1 || max_len <= 1){
                /* Nothing to align: one fragment, or every fragment is a
                   single residue (one column). */
                w = max_len;
                for(i = 0; i < n; i++){
                        MMALLOC(blk[i], w + 1);
                        memset(blk[i], '-', w);
                        blk[i][w] = '\0';
                        if(flen[i] > 0){
                                memcpy(blk[i], frag[i], flen[i]);
                        }
                }
        }else{
                RUN(kalign_arr_to_msa(in, in_len, k, &sub));
                sub->quiet = 1;
                sub->biotype = parent->biotype;
                /* The progressive aligner submits work to msa->pool
                   unconditionally when built with the threadpool. */
                sub->pool = parent->pool;
                RUN(kalign_single_run(sub, cfg, n_threads));
                w = sub->alnlen;
                for(i = 0; i < n; i++){
                        MMALLOC(blk[i], w + 1);
                        memset(blk[i], '-', w);
                        blk[i][w] = '\0';
                }
                /* kalign_single_run restores input order (msa_sort_rank) */
                for(i = 0; i < k; i++){
                        ASSERT((int)strlen(sub->sequences[i]->seq) == w,
                               "Segment row %d has length %d, expected %d",
                               i, (int)strlen(sub->sequences[i]->seq), w);
                        memcpy(blk[map[i]], sub->sequences[i]->seq, w);
                }
                kalign_free_msa(sub);
                sub = NULL;
        }

        *out_w = w;
        MFREE(in);
        MFREE(in_len);
        MFREE(map);
        return OK;
ERROR:
        if(sub) kalign_free_msa(sub);
        for(i = 0; i < n; i++){
                if(blk[i]){
                        MFREE(blk[i]);
                        blk[i] = NULL;
                }
        }
        if(in) MFREE(in);
        if(in_len) MFREE(in_len);
        if(map) MFREE(map);
        return FAIL;
}

static int cmp_int(const void* a, const void* b)
{
        int x = *(const int*)a;
        int y = *(const int*)b;
        return (x > y) - (x < y);
}

/* Project a guide alignment onto the segment: residue r0[i]..r0[i]+flen[i]-1
   of sequence i goes to the column it occupies in the guide, and guide
   columns holding none of the segment residues are dropped.  Because each
   sequence contributes a contiguous stretch, the projection is a valid
   alignment of the fragments. */
static int project_segment(int** g_r2c, int* r0, char** frag, int* flen, int n,
                           char** blk, int* out_w)
{
        int* cols = NULL;
        int total = 0;
        int m = 0;
        int i, r;

        for(i = 0; i < n; i++){
                blk[i] = NULL;
                total += flen[i];
        }
        MMALLOC(cols, sizeof(int) * (total + 1));
        for(i = 0; i < n; i++){
                for(r = 0; r < flen[i]; r++){
                        cols[m++] = g_r2c[i][r0[i] + r];
                }
        }
        qsort(cols, m, sizeof(int), cmp_int);
        {
                int u = 0;
                for(r = 0; r < m; r++){
                        if(u == 0 || cols[u - 1] != cols[r]){
                                cols[u++] = cols[r];
                        }
                }
                m = u;
        }
        for(i = 0; i < n; i++){
                MMALLOC(blk[i], m + 1);
                memset(blk[i], '-', m);
                blk[i][m] = '\0';
                for(r = 0; r < flen[i]; r++){
                        int key = g_r2c[i][r0[i] + r];
                        int* hit = bsearch(&key, cols, m, sizeof(int), cmp_int);
                        ASSERT(hit != NULL, "Guide column not found");
                        blk[i][hit - cols] = frag[i][r];
                }
        }
        MFREE(cols);
        *out_w = m;
        return OK;
ERROR:
        for(i = 0; i < n; i++){
                if(blk[i]){
                        MFREE(blk[i]);
                        blk[i] = NULL;
                }
        }
        if(cols) MFREE(cols);
        return FAIL;
}

/* Append src (length len) to the growing row *dst (capacity *cap). */
static int row_append(char** dst, int* used, int* cap, const char* src, int len)
{
        if(*used + len + 1 > *cap){
                int nc = *cap;
                while(*used + len + 1 > nc){
                        nc *= 2;
                }
                MREALLOC(*dst, nc);
                *cap = nc;
        }
        memcpy(*dst + *used, src, len);
        *used += len;
        (*dst)[*used] = '\0';
        return OK;
ERROR:
        return FAIL;
}

int consensus_fill_in(struct msa* cons,
                      const struct kalign_run_config* cfg,
                      const struct msa* guide,
                      int mode,
                      float anchor_frac,
                      int n_threads)
{
        struct kalign_run_config seg_cfg;
        int numseq;
        int alnlen;
        int thr;
        int* count = NULL;
        char** rows = NULL;
        int* used = NULL;
        int* cap = NULL;
        char** frag = NULL;
        int* flen = NULL;
        char** blk = NULL;
        int** g_r2c = NULL;
        int* rpos = NULL;
        int i, c;

        ASSERT(cons != NULL, "No consensus alignment");
        ASSERT(cfg != NULL, "No run config");
        ASSERT(mode == KALIGN_FILL_REALIGN || mode == KALIGN_FILL_PROJECT,
               "Unknown fill-in mode %d", mode);
        ASSERT(mode != KALIGN_FILL_PROJECT || guide != NULL,
               "Projection fill-in needs a guide alignment");

        numseq = cons->numseq;
        alnlen = cons->alnlen;

        /* Plain progressive run: no refinement, no tree perturbation,
           no anchor consistency (segments are short and local). */
        seg_cfg = *cfg;
        seg_cfg.refine = KALIGN_REFINE_NONE;
        seg_cfg.adaptive_budget = 0;
        seg_cfg.realign = 0;
        seg_cfg.tree_seed = 0;
        seg_cfg.tree_noise = 0.0f;
        seg_cfg.consistency_anchors = 0;

        thr = (int)ceilf(anchor_frac * (float)numseq);
        if(thr < 2) thr = 2;

        MMALLOC(count, sizeof(int) * (alnlen + 1));
        for(c = 0; c < alnlen; c++){
                count[c] = 0;
        }
        for(i = 0; i < numseq; i++){
                const char* s = cons->sequences[i]->seq;
                ASSERT((int)strlen(s) == alnlen, "Consensus row %d is not alnlen wide", i);
                for(c = 0; c < alnlen; c++){
                        if(s[c] != '-'){
                                count[c]++;
                        }
                }
        }

        MMALLOC(rpos, sizeof(int) * numseq);
        for(i = 0; i < numseq; i++){
                rpos[i] = 0;
        }
        if(mode == KALIGN_FILL_PROJECT){
                ASSERT(guide->numseq == numseq, "Guide has %d sequences, consensus %d",
                       guide->numseq, numseq);
                MMALLOC(g_r2c, sizeof(int*) * numseq);
                for(i = 0; i < numseq; i++){
                        g_r2c[i] = NULL;
                }
                for(i = 0; i < numseq; i++){
                        const char* g = guide->sequences[i]->seq;
                        int glen = (int)strlen(g);
                        int r = 0;
                        MMALLOC(g_r2c[i], sizeof(int) * (glen + 1));
                        for(c = 0; c < glen; c++){
                                if(g[c] != '-'){
                                        g_r2c[i][r++] = c;
                                }
                        }
                }
        }

        MMALLOC(rows, sizeof(char*) * numseq);
        MMALLOC(used, sizeof(int) * numseq);
        MMALLOC(cap, sizeof(int) * numseq);
        MMALLOC(frag, sizeof(char*) * numseq);
        MMALLOC(flen, sizeof(int) * numseq);
        MMALLOC(blk, sizeof(char*) * numseq);
        for(i = 0; i < numseq; i++){
                rows[i] = NULL;
                frag[i] = NULL;
                blk[i] = NULL;
        }
        for(i = 0; i < numseq; i++){
                cap[i] = alnlen + 1;
                used[i] = 0;
                MMALLOC(rows[i], cap[i]);
                rows[i][0] = '\0';
                MMALLOC(frag[i], alnlen + 1);
        }

        c = 0;
        while(c < alnlen){
                if(count[c] >= thr){
                        /* Anchor column: copy unchanged */
                        for(i = 0; i < numseq; i++){
                                RUN(row_append(&rows[i], &used[i], &cap[i],
                                               cons->sequences[i]->seq + c, 1));
                                if(cons->sequences[i]->seq[c] != '-'){
                                        rpos[i]++;
                                }
                        }
                        c++;
                        continue;
                }
                /* Segment [c, e) of non-anchor columns */
                int e = c;
                while(e < alnlen && count[e] < thr){
                        e++;
                }
                for(i = 0; i < numseq; i++){
                        const char* s = cons->sequences[i]->seq;
                        int n = 0;
                        for(int p = c; p < e; p++){
                                if(s[p] != '-'){
                                        frag[i][n++] = s[p];
                                }
                        }
                        frag[i][n] = '\0';
                        flen[i] = n;
                }
                int w = 0;
                if(mode == KALIGN_FILL_PROJECT){
                        RUN(project_segment(g_r2c, rpos, frag, flen, numseq, blk, &w));
                }else{
                        RUN(align_segment(frag, flen, numseq, cons,
                                          &seg_cfg, n_threads, blk, &w));
                }
                for(i = 0; i < numseq; i++){
                        rpos[i] += flen[i];
                }
                for(i = 0; i < numseq; i++){
                        RUN(row_append(&rows[i], &used[i], &cap[i], blk[i], w));
                        MFREE(blk[i]);
                        blk[i] = NULL;
                }
                c = e;
        }

        /* All rows must share one width */
        for(i = 1; i < numseq; i++){
                ASSERT(used[i] == used[0], "Filled rows differ in width (%d vs %d)",
                       used[i], used[0]);
        }

        for(i = 0; i < numseq; i++){
                MFREE(cons->sequences[i]->seq);
                cons->sequences[i]->seq = rows[i];
                cons->sequences[i]->len = used[i];
                cons->sequences[i]->alloc_len = cap[i];
                rows[i] = NULL;
        }
        cons->alnlen = used[0];
        cons->aligned = ALN_STATUS_FINAL;

        for(i = 0; i < numseq; i++){
                MFREE(frag[i]);
        }
        MFREE(frag);
        MFREE(flen);
        MFREE(blk);
        MFREE(rows);
        MFREE(used);
        MFREE(cap);
        MFREE(count);
        MFREE(rpos);
        if(g_r2c){
                for(i = 0; i < numseq; i++){
                        MFREE(g_r2c[i]);
                }
                MFREE(g_r2c);
        }
        return OK;
ERROR:
        if(rows){
                for(i = 0; i < numseq; i++){
                        if(rows[i]) MFREE(rows[i]);
                }
                MFREE(rows);
        }
        if(frag){
                for(i = 0; i < numseq; i++){
                        if(frag[i]) MFREE(frag[i]);
                }
                MFREE(frag);
        }
        if(blk){
                for(i = 0; i < numseq; i++){
                        if(blk[i]) MFREE(blk[i]);
                }
                MFREE(blk);
        }
        if(flen) MFREE(flen);
        if(used) MFREE(used);
        if(cap) MFREE(cap);
        if(count) MFREE(count);
        if(rpos) MFREE(rpos);
        if(g_r2c){
                for(i = 0; i < numseq; i++){
                        if(g_r2c[i]) MFREE(g_r2c[i]);
                }
                MFREE(g_r2c);
        }
        return FAIL;
}

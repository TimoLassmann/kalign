#include <math.h>
#include <string.h>

#include "tldevel.h"
#include "probmsa_logmath.h"
#include "probmsa_emission.h"

#define PROBMSA_PAIRHMM_IMPORT
#include "probmsa_pairhmm.h"

static double** alloc_2d(int rows, int cols);
static void free_2d(double** m, int rows);

/* Position-specific gap-open: delta_eff = delta_base * cs / (conc + cs).
   High concentration → small delta → hard to gap.
   Low concentration → delta → δ_base → easy to gap.
   Returns 0 if position-specific modulation is disabled (conc_scale <= 0). */
static inline double pos_delta(double conc_a, double conc_b,
                                double delta_base, double conc_scale)
{
        if(conc_scale <= 0.0){
                return delta_base;
        }
        double conc = 0.5 * (conc_a + conc_b);
        double d = delta_base * conc_scale / (conc + conc_scale);
        /* Clamp to valid range */
        if(d >= 0.499) d = 0.499;
        if(d < 1e-15) d = 1e-15;
        return d;
}

int probmsa_dp_alloc(struct probmsa_dp** out, int max_a, int max_b,
                      int alpha_size)
{
        struct probmsa_dp* dp = NULL;
        int rows = max_a + 1;
        int cols = max_b + 1;

        MMALLOC(dp, sizeof(struct probmsa_dp));
        memset(dp, 0, sizeof(struct probmsa_dp));

        dp->M  = alloc_2d(rows, cols);
        dp->X  = alloc_2d(rows, cols);
        dp->Y  = alloc_2d(rows, cols);
        dp->bM = alloc_2d(rows, cols);
        dp->bX = alloc_2d(rows, cols);
        dp->bY = alloc_2d(rows, cols);

        if(!dp->M || !dp->X || !dp->Y || !dp->bM || !dp->bX || !dp->bY){
                ERROR_MSG("Failed to allocate DP matrices");
        }

        /* 5-state matrices are NULL for 3-state model */
        dp->Xl = NULL;  dp->Yl = NULL;
        dp->bXl = NULL; dp->bYl = NULL;

        dp->alpha_size = alpha_size;
        dp->alloc_a = max_a;
        dp->alloc_b = max_b;
        dp->mea_threshold = 0.0;
        dp->mea_gpo = 0.0;
        dp->mea_gpe = 0.0;
        dp->use_dotproduct = 0;
        dp->delta_base = 0.0;
        dp->epsilon = 0.0;
        dp->conc_scale = 0.0;
        dp->use_5state = 0;

        *out = dp;
        return OK;
ERROR:
        probmsa_dp_free(dp);
        return FAIL;
}

int probmsa_dp_set_transitions(struct probmsa_dp* dp,
                                double delta, double epsilon,
                                double tau)
{
        ASSERT(dp != NULL, "NULL dp");
        ASSERT(delta > 0.0 && delta < 0.5, "delta must be in (0, 0.5)");
        ASSERT(epsilon > 0.0 && epsilon < 1.0, "epsilon must be in (0, 1)");

        dp->delta_base = delta;
        dp->epsilon = epsilon;

        dp->MM = log(1.0 - 2.0 * delta);
        dp->MX = log(delta);
        dp->MY = log(delta);
        dp->XM = log(1.0 - epsilon);
        dp->XX = log(epsilon);
        dp->YM = log(1.0 - epsilon);
        dp->YY = log(epsilon);

        if(tau > 0.0){
                dp->SX = log(tau);
                dp->SY = log(tau);
                dp->semiglobal = 1;
        }else{
                dp->SX = dp->MX;  /* same as internal delta */
                dp->SY = dp->MY;
                dp->semiglobal = 0;
        }

        return OK;
ERROR:
        return FAIL;
}

void probmsa_dp_free(struct probmsa_dp* dp)
{
        if(!dp) return;
        int rows = dp->alloc_a + 1;
        free_2d(dp->M,  rows);
        free_2d(dp->X,  rows);
        free_2d(dp->Y,  rows);
        free_2d(dp->bM, rows);
        free_2d(dp->bX, rows);
        free_2d(dp->bY, rows);
        /* 5-state extra matrices (NULL-safe) */
        free_2d(dp->Xl,  rows);
        free_2d(dp->Yl,  rows);
        free_2d(dp->bXl, rows);
        free_2d(dp->bYl, rows);
        MFREE(dp);
}

/* Forward algorithm with position-specific gap transitions.
   When use_dotproduct=1, passes NULL for joint matrix in emit_match.
   When conc_scale > 0, computes position-specific delta from profile
   concentrations at each M cell. */
int probmsa_forward(struct probmsa_dp* dp,
                     const struct probmsa_profile* a,
                     const struct probmsa_profile* b,
                     double* score)
{
        int i, j;
        int len_a = a->length;
        int len_b = b->length;
        double** M = dp->M;
        double** X = dp->X;
        double** Y = dp->Y;

        /* Global transitions (used for X→M, X→X, Y→M, Y→Y) */
        const double tXM = dp->XM;
        const double tXX = dp->XX;
        const double tYM = dp->YM;
        const double tYY = dp->YY;
        const int sg = dp->semiglobal;
        const int dotprod = dp->use_dotproduct;
        const int posspec = (dp->conc_scale > 0.0);

        /* For position-specific delta */
        const double delta_base = dp->delta_base;
        const double conc_scale = dp->conc_scale;

        /* For global fallback */
        const double tMM = dp->MM;
        const double tMX = dp->MX;
        const double tMY = dp->MY;

        const double (*joint_ptr)[PROBMSA_ALPHA_MAX] = dotprod ? NULL : dp->joint;
        const double* bg_ptr = dotprod ? NULL : dp->background;

        /* Initialize all boundary to -INF */
        for(i = 0; i <= len_a; i++){
                M[i][0] = -INFINITY;
                X[i][0] = -INFINITY;
                Y[i][0] = -INFINITY;
        }
        for(j = 0; j <= len_b; j++){
                M[0][j] = -INFINITY;
                X[0][j] = -INFINITY;
                Y[0][j] = -INFINITY;
        }

        /* Leading gaps in A (first row, Y state) */
        for(j = 1; j <= len_b; j++){
                double emit_b = probmsa_emit_gap(b, j - 1, bg_ptr, dp->alpha_size);
                if(j == 1){
                        Y[0][j] = dp->SY + emit_b;
                }else{
                        Y[0][j] = Y[0][j-1] + tYY + emit_b;
                }
        }

        /* Leading gaps in B (first column, X state) */
        for(i = 1; i <= len_a; i++){
                double emit_a = probmsa_emit_gap(a, i - 1, bg_ptr, dp->alpha_size);
                if(i == 1){
                        X[i][0] = dp->SX + emit_a;
                }else{
                        X[i][0] = X[i-1][0] + tXX + emit_a;
                }
        }

        /* Main DP */
        for(i = 1; i <= len_a; i++){
                double emit_a = probmsa_emit_gap(a, i - 1, bg_ptr, dp->alpha_size);

                for(j = 1; j <= len_b; j++){
                        double emit_m = probmsa_emit_match(a, i - 1, b, j - 1,
                                                            joint_ptr, dp->alpha_size);
                        double emit_b = probmsa_emit_gap(b, j - 1, bg_ptr, dp->alpha_size);

                        /* Position-specific transitions from M cells */
                        double t_mm_here, t_mx_here, t_my_here;

                        if(posspec){
                                /* M→M transition from M[i-1][j-1]: uses delta at (i-1,j-1) */
                                if(i >= 2 && j >= 2){
                                        double d = pos_delta(a->cols[i-2].concentration,
                                                              b->cols[j-2].concentration,
                                                              delta_base, conc_scale);
                                        t_mm_here = log(1.0 - 2.0 * d);
                                }else{
                                        t_mm_here = tMM;
                                }

                                /* M→X transition from M[i-1][j]: uses delta at (i-1,j) */
                                if(i >= 2){
                                        double d = pos_delta(a->cols[i-2].concentration,
                                                              b->cols[j-1].concentration,
                                                              delta_base, conc_scale);
                                        t_mx_here = sg && j == len_b ? dp->SX : log(d);
                                }else{
                                        t_mx_here = sg && j == len_b ? dp->SX : tMX;
                                }

                                /* M→Y transition from M[i][j-1]: uses delta at (i,j-1) */
                                if(j >= 2){
                                        double d = pos_delta(a->cols[i-1].concentration,
                                                              b->cols[j-2].concentration,
                                                              delta_base, conc_scale);
                                        t_my_here = sg && i == len_a ? dp->SY : log(d);
                                }else{
                                        t_my_here = sg && i == len_a ? dp->SY : tMY;
                                }
                        }else{
                                /* Global transitions */
                                t_mx_here = (sg && j == len_b) ? dp->SX : tMX;
                                t_my_here = (sg && i == len_a) ? dp->SY : tMY;
                                t_mm_here = tMM;
                        }

                        /* Match state */
                        if(i == 1 && j == 1){
                                M[i][j] = emit_m;
                        }else{
                                M[i][j] = LOGSUM3(M[i-1][j-1] + t_mm_here,
                                                  X[i-1][j-1] + tXM,
                                                  Y[i-1][j-1] + tYM) + emit_m;
                        }

                        /* X state: gap in B (advance A) */
                        X[i][j] = logsum(X[i-1][j] + tXX,
                                         M[i-1][j] + t_mx_here) + emit_a;

                        /* Y state: gap in A (advance B) */
                        Y[i][j] = logsum(Y[i][j-1] + tYY,
                                         M[i][j-1] + t_my_here) + emit_b;
                }
        }

        *score = LOGSUM3(M[len_a][len_b],
                         X[len_a][len_b],
                         Y[len_a][len_b]);
        return OK;
}

/* Backward algorithm with position-specific gap transitions. */
int probmsa_backward(struct probmsa_dp* dp,
                      const struct probmsa_profile* a,
                      const struct probmsa_profile* b,
                      double* score)
{
        int i, j;
        int len_a = a->length;
        int len_b = b->length;
        double** M = dp->bM;
        double** X = dp->bX;
        double** Y = dp->bY;

        const double tXM = dp->XM;
        const double tXX = dp->XX;
        const double tYM = dp->YM;
        const double tYY = dp->YY;
        const int sg = dp->semiglobal;
        const int dotprod = dp->use_dotproduct;
        const int posspec = (dp->conc_scale > 0.0);
        const double delta_base = dp->delta_base;
        const double conc_scale = dp->conc_scale;
        const double tMM = dp->MM;
        const double tMX = dp->MX;
        const double tMY = dp->MY;

        const double (*joint_ptr)[PROBMSA_ALPHA_MAX] = dotprod ? NULL : dp->joint;
        const double* bg_ptr = dotprod ? NULL : dp->background;

        /* Initialize all to -INF */
        for(i = 0; i <= len_a; i++){
                for(j = 0; j <= len_b; j++){
                        M[i][j] = -INFINITY;
                        X[i][j] = -INFINITY;
                        Y[i][j] = -INFINITY;
                }
        }

        /* Last cell */
        M[len_a][len_b] = probmsa_emit_match(a, len_a - 1, b, len_b - 1,
                                               joint_ptr, dp->alpha_size);
        X[len_a][len_b] = probmsa_emit_gap(a, len_a - 1, bg_ptr, dp->alpha_size);
        Y[len_a][len_b] = probmsa_emit_gap(b, len_b - 1, bg_ptr, dp->alpha_size);

        /* Last row (i = len_a): trailing gaps in A */
        for(j = len_b - 1; j >= 1; j--){
                double emit_m = probmsa_emit_match(a, len_a - 1, b, j - 1,
                                                    joint_ptr, dp->alpha_size);
                double emit_y = probmsa_emit_gap(b, j - 1, bg_ptr, dp->alpha_size);

                /* M→Y transition from M[len_a][j]: uses delta at (len_a, j) */
                double t_my;
                if(sg){
                        t_my = dp->SY;
                }else if(posspec){
                        double d = pos_delta(a->cols[len_a-1].concentration,
                                              b->cols[j-1].concentration,
                                              delta_base, conc_scale);
                        t_my = log(d);
                }else{
                        t_my = tMY;
                }

                M[len_a][j] = emit_m + t_my + Y[len_a][j+1];
                X[len_a][j] = -INFINITY;
                Y[len_a][j] = emit_y + tYY + Y[len_a][j+1];
        }

        /* Last column (j = len_b): trailing gaps in B */
        for(i = len_a - 1; i >= 1; i--){
                double emit_m = probmsa_emit_match(a, i - 1, b, len_b - 1,
                                                    joint_ptr, dp->alpha_size);
                double emit_x = probmsa_emit_gap(a, i - 1, bg_ptr, dp->alpha_size);

                /* M→X transition from M[i][len_b]: uses delta at (i, len_b) */
                double t_mx;
                if(sg){
                        t_mx = dp->SX;
                }else if(posspec){
                        double d = pos_delta(a->cols[i-1].concentration,
                                              b->cols[len_b-1].concentration,
                                              delta_base, conc_scale);
                        t_mx = log(d);
                }else{
                        t_mx = tMX;
                }

                M[i][len_b] = emit_m + t_mx + X[i+1][len_b];
                X[i][len_b] = emit_x + tXX + X[i+1][len_b];
                Y[i][len_b] = -INFINITY;
        }

        /* Interior cells */
        for(i = len_a - 1; i >= 1; i--){
                double emit_x = probmsa_emit_gap(a, i - 1, bg_ptr, dp->alpha_size);

                for(j = len_b - 1; j >= 1; j--){
                        double emit_m = probmsa_emit_match(a, i - 1, b, j - 1,
                                                            joint_ptr, dp->alpha_size);
                        double emit_y = probmsa_emit_gap(b, j - 1, bg_ptr, dp->alpha_size);

                        /* Outgoing transitions from M[i][j] */
                        double t_mm_out, t_mx_out, t_my_out;
                        if(posspec){
                                double d = pos_delta(a->cols[i-1].concentration,
                                                      b->cols[j-1].concentration,
                                                      delta_base, conc_scale);
                                t_mm_out = log(1.0 - 2.0 * d);
                                t_mx_out = log(d);
                                t_my_out = log(d);
                        }else{
                                t_mm_out = tMM;
                                t_mx_out = tMX;
                                t_my_out = tMY;
                        }

                        M[i][j] = LOGSUM3(M[i+1][j+1] + t_mm_out,
                                          X[i+1][j]   + t_mx_out,
                                          Y[i][j+1]   + t_my_out) + emit_m;

                        X[i][j] = logsum(M[i+1][j+1] + tXM,
                                         X[i+1][j]   + tXX) + emit_x;

                        Y[i][j] = logsum(M[i+1][j+1] + tYM,
                                         Y[i][j+1]   + tYY) + emit_y;
                }

                /* Column j=0: only X state valid (leading gaps in B) */
                {
                        double emit_x2 = probmsa_emit_gap(a, i - 1, bg_ptr, dp->alpha_size);
                        X[i][0] = logsum(M[i+1][1] + tXM,
                                         X[i+1][0] + tXX) + emit_x2;
                }
        }

        /* Row i=0: only Y state valid (leading gaps in A) */
        for(j = len_b - 1; j >= 1; j--){
                double emit_y = probmsa_emit_gap(b, j - 1, bg_ptr, dp->alpha_size);
                Y[0][j] = logsum(M[1][j+1] + tYM,
                                 Y[0][j+1] + tYY) + emit_y;
        }

        /* Corner (0, len_b) */
        M[0][len_b] = -INFINITY;
        X[0][len_b] = -INFINITY;
        Y[0][len_b] = -INFINITY;

        /* Backward score */
        {
                double sm = 0.0;
                *score = LOGSUM3(sm  + M[1][1],
                                 dp->SX + X[1][0],
                                 dp->SY + Y[0][1]);
        }
        return OK;
}

/* Compute posterior probabilities.
   Overwrites forward M with log P(Match at i,j | data). */
int probmsa_posterior(struct probmsa_dp* dp,
                       const struct probmsa_profile* a,
                       const struct probmsa_profile* b,
                       int len_a, int len_b, double total)
{
        int i, j;
        double** fM = dp->M;
        double** bM = dp->bM;
        const double (*joint_ptr)[PROBMSA_ALPHA_MAX] = dp->use_dotproduct ? NULL : dp->joint;

        /* Set boundary cells to -INF */
        fM[0][0] = -INFINITY;
        for(j = 1; j <= len_b; j++){
                fM[0][j] = -INFINITY;
        }

        for(i = 1; i <= len_a; i++){
                fM[i][0] = -INFINITY;
                for(j = 1; j <= len_b; j++){
                        double emit = probmsa_emit_match(
                                a, i - 1, b, j - 1,
                                joint_ptr, dp->alpha_size);
                        fM[i][j] = fM[i][j] + bM[i][j] - emit - total;
                }
        }

        return OK;
}

/* MEA traceback with affine gap penalties, producing kalign-convention path.
   Before MEA DP, converts match posteriors from log-space to probability
   and saves them in dp->bY for later use in profile merging.
   Uses three MEA states: Sm (match), Ga (gap-in-A), Gb (gap-in-B).
     Sm[i][j] = max(Sm,Ga,Gb)[i-1][j-1] + (post(i,j) - threshold)
     Ga[i][j] = max(Sm[i][j-1] - gpo, Ga[i][j-1] - gpe, Gb[i][j-1] - gpo)
     Gb[i][j] = max(Sm[i-1][j] - gpo, Ga[i-1][j] - gpo, Gb[i-1][j] - gpe)
   Matrices used: M=Sm, X=Ga, bX=Gb scores; bM=Sm tb, Y=Ga tb, bY posteriors.
   Path format: path[0]=aln_len, path[1..N]=ops (0/1/2), path[N+1]=3. */
int probmsa_mea_traceback(struct probmsa_dp* dp, int len_a, int len_b,
                            int** path)
{
        double** Sm = dp->M;    /* posteriors → overwritten with match-state MEA scores */
        double** Ga = dp->X;    /* gap-in-A state scores (advance B) */
        double** Gb = dp->bX;   /* gap-in-B state scores (advance A) */
        double** tb_m = dp->bM; /* traceback for Sm: 0=Sm, 1=Ga, 2=Gb */
        double** tb_a = dp->Y;  /* traceback for Ga: 0=Sm, 1=Ga, 2=Gb */
        double** tb_b = dp->bY; /* saved posteriors (phase 1) + traceback for Gb (phase 2) */
        int* p = NULL;
        int i, j, c;
        int state;
        double threshold = dp->mea_threshold;
        double gpo = dp->mea_gpo;
        double gpe = dp->mea_gpe;

        /* Phase 1: Convert posteriors to probabilities, save in temp storage */
        double** post_save = NULL;
        post_save = (double**)calloc(len_a + 1, sizeof(double*));
        if(!post_save) goto ERROR;
        for(i = 0; i <= len_a; i++){
                post_save[i] = (double*)malloc(sizeof(double) * (len_b + 1));
                if(!post_save[i]) goto ERROR;
        }

        for(i = 0; i <= len_a; i++){
                for(j = 0; j <= len_b; j++){
                        double post = log2prob(Sm[i][j]);
                        post_save[i][j] = post;
                }
        }

        /* Phase 2: Three-state affine MEA DP */
        for(i = 0; i <= len_a; i++){
                Sm[i][0] = -INFINITY;
                Ga[i][0] = -INFINITY;
                Gb[i][0] = -INFINITY;
                tb_m[i][0] = 0.0;
                tb_a[i][0] = 0.0;
                tb_b[i][0] = 0.0;
        }
        for(j = 0; j <= len_b; j++){
                Sm[0][j] = -INFINITY;
                Ga[0][j] = -INFINITY;
                Gb[0][j] = -INFINITY;
                tb_m[0][j] = 0.0;
                tb_a[0][j] = 0.0;
                tb_b[0][j] = 0.0;
        }
        Sm[0][0] = 0.0;

        /* First row: can only advance B (gap-in-A) */
        for(j = 1; j <= len_b; j++){
                Sm[0][j] = -INFINITY;
                Gb[0][j] = -INFINITY;
                if(j == 1){
                        Ga[0][j] = Sm[0][0] - gpo;
                        tb_a[0][j] = 0.0;
                }else{
                        double from_ga = Ga[0][j-1] - gpe;
                        Ga[0][j] = from_ga;
                        tb_a[0][j] = 1.0;
                }
        }

        /* First column: can only advance A (gap-in-B) */
        for(i = 1; i <= len_a; i++){
                Sm[i][0] = -INFINITY;
                Ga[i][0] = -INFINITY;
                if(i == 1){
                        Gb[i][0] = Sm[0][0] - gpo;
                        tb_b[i][0] = 0.0;
                }else{
                        double from_gb = Gb[i-1][0] - gpe;
                        Gb[i][0] = from_gb;
                        tb_b[i][0] = 2.0;
                }
        }

        /* Interior cells */
        for(i = 1; i <= len_a; i++){
                for(j = 1; j <= len_b; j++){
                        double post = post_save[i][j] - threshold;

                        /* Sm[i][j] = max(Sm,Ga,Gb)[i-1][j-1] + post */
                        {
                                double v0 = Sm[i-1][j-1];
                                double v1 = Ga[i-1][j-1];
                                double v2 = Gb[i-1][j-1];
                                double best = v0; double from = 0.0;
                                if(v1 > best){ best = v1; from = 1.0; }
                                if(v2 > best){ best = v2; from = 2.0; }
                                Sm[i][j] = best + post;
                                tb_m[i][j] = from;
                        }

                        /* Ga[i][j]: gap-in-A = advance B only */
                        {
                                double from_sm = Sm[i][j-1] - gpo;
                                double from_ga = Ga[i][j-1] - gpe;
                                double from_gb = Gb[i][j-1] - gpo;
                                double best = from_sm; double from = 0.0;
                                if(from_ga > best){ best = from_ga; from = 1.0; }
                                if(from_gb > best){ best = from_gb; from = 2.0; }
                                Ga[i][j] = best;
                                tb_a[i][j] = from;
                        }

                        /* Gb[i][j]: gap-in-B = advance A only */
                        {
                                double from_sm = Sm[i-1][j] - gpo;
                                double from_ga = Ga[i-1][j] - gpo;
                                double from_gb = Gb[i-1][j] - gpe;
                                double best = from_sm; double from = 0.0;
                                if(from_ga > best){ best = from_ga; from = 1.0; }
                                if(from_gb > best){ best = from_gb; from = 2.0; }
                                Gb[i][j] = best;
                                tb_b[i][j] = from;
                        }
                }
        }

        /* Phase 3: Traceback */
        MMALLOC(p, sizeof(int) * (len_a + len_b + 3));

        i = len_a;
        j = len_b;
        {
                double best = Sm[i][j]; state = 0;
                if(Ga[i][j] > best){ best = Ga[i][j]; state = 1; }
                if(Gb[i][j] > best){ best = Gb[i][j]; state = 2; }
        }

        c = 1;
        while(i + j > 0){
                if(state == 0){
                        int prev = (int)tb_m[i][j];
                        p[c] = 0;
                        i--;
                        j--;
                        state = prev;
                }else if(state == 1){
                        int prev = (int)tb_a[i][j];
                        p[c] = 1;
                        j--;
                        state = prev;
                }else{
                        int prev = (int)tb_b[i][j];
                        p[c] = 2;
                        i--;
                        state = prev;
                }
                c++;
        }

        int aln_len = c - 1;
        p[0] = aln_len;

        /* Reverse the path */
        for(i = 1, j = aln_len; i < j; i++, j--){
                state = p[i];
                p[i] = p[j];
                p[j] = state;
        }

        p[aln_len + 1] = 3;

        /* Copy posteriors to dp->Y for profile merge */
        for(i = 0; i <= len_a; i++){
                for(j = 0; j <= len_b; j++){
                        dp->Y[i][j] = post_save[i][j];
                }
        }
        for(i = 0; i <= len_a; i++) free(post_save[i]);
        free(post_save);
        post_save = NULL;

        *path = p;
        return OK;
ERROR:
        if(post_save){
                for(i = 0; i <= len_a; i++){
                        if(post_save[i]) free(post_save[i]);
                }
                free(post_save);
        }
        if(p) MFREE(p);
        return FAIL;
}

/* ========================================================================
 * 5-state pair-HMM: M, Xs (short gap-A), Ys (short gap-B),
 *                      Xl (long gap-A),  Yl (long gap-B)
 *
 * Short gaps: frequent (high δs), low self-transition (low εs) → mean length ~2
 * Long gaps:  rare (low δl), high self-transition (high εl) → mean length ~20
 *
 * Reuses M, X=Xs, Y=Ys, bM, bX=bXs, bY=bYs from the struct.
 * Adds Xl, Yl, bXl, bYl for long gap states.
 * MEA traceback is reused unchanged (operates on match posteriors only).
 * ======================================================================== */

int probmsa_dp_alloc_5state(struct probmsa_dp** out, int max_a, int max_b,
                              int alpha_size)
{
        struct probmsa_dp* dp = NULL;
        int rows = max_a + 1;
        int cols = max_b + 1;

        MMALLOC(dp, sizeof(struct probmsa_dp));
        memset(dp, 0, sizeof(struct probmsa_dp));

        /* Shared with 3-state: M, Xs=X, Ys=Y */
        dp->M  = alloc_2d(rows, cols);
        dp->X  = alloc_2d(rows, cols);
        dp->Y  = alloc_2d(rows, cols);
        dp->bM = alloc_2d(rows, cols);
        dp->bX = alloc_2d(rows, cols);
        dp->bY = alloc_2d(rows, cols);

        /* 5-state additional: Xl, Yl */
        dp->Xl  = alloc_2d(rows, cols);
        dp->Yl  = alloc_2d(rows, cols);
        dp->bXl = alloc_2d(rows, cols);
        dp->bYl = alloc_2d(rows, cols);

        if(!dp->M || !dp->X || !dp->Y || !dp->bM || !dp->bX || !dp->bY ||
           !dp->Xl || !dp->Yl || !dp->bXl || !dp->bYl){
                ERROR_MSG("Failed to allocate 5-state DP matrices");
        }

        dp->alpha_size = alpha_size;
        dp->alloc_a = max_a;
        dp->alloc_b = max_b;
        dp->mea_threshold = 0.0;
        dp->mea_gpo = 0.0;
        dp->mea_gpe = 0.0;
        dp->use_dotproduct = 0;
        dp->delta_base = 0.0;
        dp->epsilon = 0.0;
        dp->conc_scale = 0.0;
        dp->use_5state = 1;

        *out = dp;
        return OK;
ERROR:
        probmsa_dp_free(dp);
        return FAIL;
}

int probmsa_dp_set_transitions_5state(struct probmsa_dp* dp,
                                        double delta_s, double epsilon_s,
                                        double delta_l, double epsilon_l)
{
        ASSERT(dp != NULL, "NULL dp");
        ASSERT(delta_s > 0.0 && delta_s < 0.5, "delta_s must be in (0, 0.5)");
        ASSERT(delta_l > 0.0 && delta_l < 0.5, "delta_l must be in (0, 0.5)");
        ASSERT(2.0 * delta_s + 2.0 * delta_l < 1.0,
               "2*delta_s + 2*delta_l must be < 1");
        ASSERT(epsilon_s > 0.0 && epsilon_s < 1.0,
               "epsilon_s must be in (0, 1)");
        ASSERT(epsilon_l > 0.0 && epsilon_l < 1.0,
               "epsilon_l must be in (0, 1)");

        dp->t5_MM  = log(1.0 - 2.0 * delta_s - 2.0 * delta_l);
        dp->t5_MXs = log(delta_s);
        dp->t5_MYs = log(delta_s);
        dp->t5_MXl = log(delta_l);
        dp->t5_MYl = log(delta_l);

        dp->t5_XsM  = log(1.0 - epsilon_s);
        dp->t5_XsXs = log(epsilon_s);
        dp->t5_YsM  = log(1.0 - epsilon_s);
        dp->t5_YsYs = log(epsilon_s);

        dp->t5_XlM  = log(1.0 - epsilon_l);
        dp->t5_XlXl = log(epsilon_l);
        dp->t5_YlM  = log(1.0 - epsilon_l);
        dp->t5_YlYl = log(epsilon_l);

        /* Start transitions: global alignment, same as internal gap-open */
        dp->t5_SXs = dp->t5_MXs;
        dp->t5_SYs = dp->t5_MYs;
        dp->t5_SXl = dp->t5_MXl;
        dp->t5_SYl = dp->t5_MYl;

        dp->use_5state = 1;

        return OK;
ERROR:
        return FAIL;
}

/* 5-state forward algorithm.
   M=match, X=Xs (short gap-A), Y=Ys (short gap-B), Xl (long gap-A), Yl (long gap-B).
   Global alignment: all gaps penalized. */
int probmsa_forward_5state(struct probmsa_dp* dp,
                             const struct probmsa_profile* a,
                             const struct probmsa_profile* b,
                             double* score)
{
        int i, j;
        int len_a = a->length;
        int len_b = b->length;
        double** M  = dp->M;
        double** Xs = dp->X;
        double** Ys = dp->Y;
        double** Xl = dp->Xl;
        double** Yl = dp->Yl;
        const int dotprod = dp->use_dotproduct;

        const double (*joint_ptr)[PROBMSA_ALPHA_MAX] = dotprod ? NULL : dp->joint;
        const double* bg_ptr = dotprod ? NULL : dp->background;

        /* Transitions */
        const double tMM   = dp->t5_MM;
        const double tMXs  = dp->t5_MXs;
        const double tMYs  = dp->t5_MYs;
        const double tMXl  = dp->t5_MXl;
        const double tMYl  = dp->t5_MYl;
        const double tXsM  = dp->t5_XsM;
        const double tXsXs = dp->t5_XsXs;
        const double tYsM  = dp->t5_YsM;
        const double tYsYs = dp->t5_YsYs;
        const double tXlM  = dp->t5_XlM;
        const double tXlXl = dp->t5_XlXl;
        const double tYlM  = dp->t5_YlM;
        const double tYlYl = dp->t5_YlYl;

        /* Initialize all boundaries to -INF */
        for(i = 0; i <= len_a; i++){
                M[i][0]  = -INFINITY;
                Xs[i][0] = -INFINITY;
                Ys[i][0] = -INFINITY;
                Xl[i][0] = -INFINITY;
                Yl[i][0] = -INFINITY;
        }
        for(j = 0; j <= len_b; j++){
                M[0][j]  = -INFINITY;
                Xs[0][j] = -INFINITY;
                Ys[0][j] = -INFINITY;
                Xl[0][j] = -INFINITY;
                Yl[0][j] = -INFINITY;
        }

        /* Leading gaps in A (first row: Y states advance B) */
        for(j = 1; j <= len_b; j++){
                double emit_b = probmsa_emit_gap(b, j - 1, bg_ptr, dp->alpha_size);
                if(j == 1){
                        Ys[0][j] = dp->t5_SYs + emit_b;
                        Yl[0][j] = dp->t5_SYl + emit_b;
                }else{
                        Ys[0][j] = Ys[0][j-1] + tYsYs + emit_b;
                        Yl[0][j] = Yl[0][j-1] + tYlYl + emit_b;
                }
        }

        /* Leading gaps in B (first column: X states advance A) */
        for(i = 1; i <= len_a; i++){
                double emit_a = probmsa_emit_gap(a, i - 1, bg_ptr, dp->alpha_size);
                if(i == 1){
                        Xs[i][0] = dp->t5_SXs + emit_a;
                        Xl[i][0] = dp->t5_SXl + emit_a;
                }else{
                        Xs[i][0] = Xs[i-1][0] + tXsXs + emit_a;
                        Xl[i][0] = Xl[i-1][0] + tXlXl + emit_a;
                }
        }

        /* Main DP */
        for(i = 1; i <= len_a; i++){
                double emit_a = probmsa_emit_gap(a, i - 1, bg_ptr, dp->alpha_size);

                for(j = 1; j <= len_b; j++){
                        double emit_m = probmsa_emit_match(a, i - 1, b, j - 1,
                                                            joint_ptr, dp->alpha_size);
                        double emit_b = probmsa_emit_gap(b, j - 1, bg_ptr, dp->alpha_size);

                        /* Match: from M, Xs, Ys, Xl, Yl at [i-1][j-1] */
                        if(i == 1 && j == 1){
                                M[i][j] = emit_m;
                        }else{
                                M[i][j] = LOGSUM5(
                                        M[i-1][j-1]  + tMM,
                                        Xs[i-1][j-1] + tXsM,
                                        Ys[i-1][j-1] + tYsM,
                                        Xl[i-1][j-1] + tXlM,
                                        Yl[i-1][j-1] + tYlM) + emit_m;
                        }

                        /* Short gap-in-B (advance A): from M or Xs at [i-1][j] */
                        Xs[i][j] = logsum(M[i-1][j]  + tMXs,
                                          Xs[i-1][j] + tXsXs) + emit_a;

                        /* Short gap-in-A (advance B): from M or Ys at [i][j-1] */
                        Ys[i][j] = logsum(M[i][j-1]  + tMYs,
                                          Ys[i][j-1] + tYsYs) + emit_b;

                        /* Long gap-in-B (advance A): from M or Xl at [i-1][j] */
                        Xl[i][j] = logsum(M[i-1][j]  + tMXl,
                                          Xl[i-1][j] + tXlXl) + emit_a;

                        /* Long gap-in-A (advance B): from M or Yl at [i][j-1] */
                        Yl[i][j] = logsum(M[i][j-1]  + tMYl,
                                          Yl[i][j-1] + tYlYl) + emit_b;
                }
        }

        *score = LOGSUM5(M[len_a][len_b],
                         Xs[len_a][len_b],
                         Ys[len_a][len_b],
                         Xl[len_a][len_b],
                         Yl[len_a][len_b]);
        return OK;
}

/* 5-state backward algorithm. */
int probmsa_backward_5state(struct probmsa_dp* dp,
                              const struct probmsa_profile* a,
                              const struct probmsa_profile* b,
                              double* score)
{
        int i, j;
        int len_a = a->length;
        int len_b = b->length;
        double** bM  = dp->bM;
        double** bXs = dp->bX;
        double** bYs = dp->bY;
        double** bXl = dp->bXl;
        double** bYl = dp->bYl;
        const int dotprod = dp->use_dotproduct;

        const double (*joint_ptr)[PROBMSA_ALPHA_MAX] = dotprod ? NULL : dp->joint;
        const double* bg_ptr = dotprod ? NULL : dp->background;

        const double tMM   = dp->t5_MM;
        const double tMXs  = dp->t5_MXs;
        const double tMYs  = dp->t5_MYs;
        const double tMXl  = dp->t5_MXl;
        const double tMYl  = dp->t5_MYl;
        const double tXsM  = dp->t5_XsM;
        const double tXsXs = dp->t5_XsXs;
        const double tYsM  = dp->t5_YsM;
        const double tYsYs = dp->t5_YsYs;
        const double tXlM  = dp->t5_XlM;
        const double tXlXl = dp->t5_XlXl;
        const double tYlM  = dp->t5_YlM;
        const double tYlYl = dp->t5_YlYl;

        /* Initialize all to -INF */
        for(i = 0; i <= len_a; i++){
                for(j = 0; j <= len_b; j++){
                        bM[i][j]  = -INFINITY;
                        bXs[i][j] = -INFINITY;
                        bYs[i][j] = -INFINITY;
                        bXl[i][j] = -INFINITY;
                        bYl[i][j] = -INFINITY;
                }
        }

        /* Last cell (len_a, len_b) */
        bM[len_a][len_b]  = probmsa_emit_match(a, len_a - 1, b, len_b - 1,
                                                  joint_ptr, dp->alpha_size);
        bXs[len_a][len_b] = probmsa_emit_gap(a, len_a - 1, bg_ptr, dp->alpha_size);
        bYs[len_a][len_b] = probmsa_emit_gap(b, len_b - 1, bg_ptr, dp->alpha_size);
        bXl[len_a][len_b] = probmsa_emit_gap(a, len_a - 1, bg_ptr, dp->alpha_size);
        bYl[len_a][len_b] = probmsa_emit_gap(b, len_b - 1, bg_ptr, dp->alpha_size);

        /* Last row (i=len_a): trailing gaps in A (only Y states can advance B) */
        for(j = len_b - 1; j >= 1; j--){
                double emit_m = probmsa_emit_match(a, len_a - 1, b, j - 1,
                                                    joint_ptr, dp->alpha_size);
                double emit_ys = probmsa_emit_gap(b, j - 1, bg_ptr, dp->alpha_size);

                /* M at last row: can only go to Ys or Yl (advance B) */
                bM[len_a][j] = emit_m + logsum(tMYs + bYs[len_a][j+1],
                                                tMYl + bYl[len_a][j+1]);

                /* X states impossible at last row (can't advance A) */
                bXs[len_a][j] = -INFINITY;
                bXl[len_a][j] = -INFINITY;

                /* Y states: can only self-loop (Y→M needs i+1 which is out of bounds) */
                bYs[len_a][j] = emit_ys + tYsYs + bYs[len_a][j+1];
                bYl[len_a][j] = emit_ys + tYlYl + bYl[len_a][j+1];
        }

        /* Last column (j=len_b): trailing gaps in B (only X states can advance A) */
        for(i = len_a - 1; i >= 1; i--){
                double emit_m = probmsa_emit_match(a, i - 1, b, len_b - 1,
                                                    joint_ptr, dp->alpha_size);
                double emit_xs = probmsa_emit_gap(a, i - 1, bg_ptr, dp->alpha_size);

                /* M at last column: can only go to Xs or Xl (advance A) */
                bM[i][len_b] = emit_m + logsum(tMXs + bXs[i+1][len_b],
                                                tMXl + bXl[i+1][len_b]);

                /* X states: can only self-loop (X→M needs j+1 which is out of bounds) */
                bXs[i][len_b] = emit_xs + tXsXs + bXs[i+1][len_b];
                bXl[i][len_b] = emit_xs + tXlXl + bXl[i+1][len_b];

                /* Y states impossible at last column */
                bYs[i][len_b] = -INFINITY;
                bYl[i][len_b] = -INFINITY;
        }

        /* Interior cells */
        for(i = len_a - 1; i >= 1; i--){
                double emit_x = probmsa_emit_gap(a, i - 1, bg_ptr, dp->alpha_size);

                for(j = len_b - 1; j >= 1; j--){
                        double emit_m = probmsa_emit_match(a, i - 1, b, j - 1,
                                                            joint_ptr, dp->alpha_size);
                        double emit_y = probmsa_emit_gap(b, j - 1, bg_ptr, dp->alpha_size);

                        /* M: outgoing to M, Xs, Ys, Xl, Yl */
                        bM[i][j] = emit_m + LOGSUM5(
                                tMM  + bM[i+1][j+1],
                                tMXs + bXs[i+1][j],
                                tMYs + bYs[i][j+1],
                                tMXl + bXl[i+1][j],
                                tMYl + bYl[i][j+1]);

                        /* Xs: outgoing to M or Xs */
                        bXs[i][j] = emit_x + logsum(tXsM  + bM[i+1][j+1],
                                                     tXsXs + bXs[i+1][j]);

                        /* Ys: outgoing to M or Ys */
                        bYs[i][j] = emit_y + logsum(tYsM  + bM[i+1][j+1],
                                                     tYsYs + bYs[i][j+1]);

                        /* Xl: outgoing to M or Xl */
                        bXl[i][j] = emit_x + logsum(tXlM  + bM[i+1][j+1],
                                                     tXlXl + bXl[i+1][j]);

                        /* Yl: outgoing to M or Yl */
                        bYl[i][j] = emit_y + logsum(tYlM  + bM[i+1][j+1],
                                                     tYlYl + bYl[i][j+1]);
                }

                /* Column j=0: leading gaps in B (only X states valid) */
                {
                        double emit_x2 = probmsa_emit_gap(a, i - 1, bg_ptr, dp->alpha_size);
                        bXs[i][0] = emit_x2 + logsum(tXsM  + bM[i+1][1],
                                                       tXsXs + bXs[i+1][0]);
                        bXl[i][0] = emit_x2 + logsum(tXlM  + bM[i+1][1],
                                                       tXlXl + bXl[i+1][0]);
                }
        }

        /* Row i=0: leading gaps in A (only Y states valid) */
        for(j = len_b - 1; j >= 1; j--){
                double emit_y = probmsa_emit_gap(b, j - 1, bg_ptr, dp->alpha_size);
                bYs[0][j] = emit_y + logsum(tYsM  + bM[1][j+1],
                                             tYsYs + bYs[0][j+1]);
                bYl[0][j] = emit_y + logsum(tYlM  + bM[1][j+1],
                                             tYlYl + bYl[0][j+1]);
        }

        /* Corners that must be -INF */
        bM[0][len_b] = -INFINITY;
        bXs[0][len_b] = -INFINITY;
        bYs[0][len_b] = -INFINITY;
        bXl[0][len_b] = -INFINITY;
        bYl[0][len_b] = -INFINITY;

        /* Backward score: sum over start states */
        {
                double sm = 0.0;  /* P(start→M) = 1 in log */
                *score = LOGSUM5(sm            + bM[1][1],
                                 dp->t5_SXs   + bXs[1][0],
                                 dp->t5_SYs   + bYs[0][1],
                                 dp->t5_SXl   + bXl[1][0],
                                 dp->t5_SYl   + bYl[0][1]);
        }
        return OK;
}

/* 5-state posterior computation.
   Overwrites forward M with log P(Match at i,j | data). */
int probmsa_posterior_5state(struct probmsa_dp* dp,
                               const struct probmsa_profile* a,
                               const struct probmsa_profile* b,
                               int len_a, int len_b, double total)
{
        int i, j;
        double** fM = dp->M;
        double** bM = dp->bM;
        const double (*joint_ptr)[PROBMSA_ALPHA_MAX] = dp->use_dotproduct ? NULL : dp->joint;

        /* Boundary cells = -INF */
        fM[0][0] = -INFINITY;
        for(j = 1; j <= len_b; j++){
                fM[0][j] = -INFINITY;
        }

        for(i = 1; i <= len_a; i++){
                fM[i][0] = -INFINITY;
                for(j = 1; j <= len_b; j++){
                        double emit = probmsa_emit_match(
                                a, i - 1, b, j - 1,
                                joint_ptr, dp->alpha_size);
                        fM[i][j] = fM[i][j] + bM[i][j] - emit - total;
                }
        }

        return OK;
}

/* --- Internal helpers --- */

static double** alloc_2d(int rows, int cols)
{
        double** m = NULL;
        int i;

        m = (double**)malloc(sizeof(double*) * rows);
        if(!m) return NULL;

        for(i = 0; i < rows; i++){
                m[i] = (double*)malloc(sizeof(double) * cols);
                if(!m[i]){
                        int k;
                        for(k = 0; k < i; k++) free(m[k]);
                        free(m);
                        return NULL;
                }
                for(int j = 0; j < cols; j++){
                        m[i][j] = -INFINITY;
                }
        }
        return m;
}

static void free_2d(double** m, int rows)
{
        if(!m) return;
        int i;
        for(i = 0; i < rows; i++){
                if(m[i]) free(m[i]);
        }
        free(m);
}

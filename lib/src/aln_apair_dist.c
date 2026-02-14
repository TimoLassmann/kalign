/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019, 2024, 2025 Timo Lassmann

    This file is part of kalign.

    Kalign is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#include "tldevel.h"
#include "msa_struct.h"

#define ALN_APAIR_DIST_IMPORT
#include "aln_apair_dist.h"

static float pairwise_identity_dist(const char* a, const char* b, int alnlen);

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

        for(i = 0; i < n - 1; i++){
                const char* seq_i = msa->sequences[i]->seq;
                for(j = i + 1; j < n; j++){
                        float d = pairwise_identity_dist(seq_i, msa->sequences[j]->seq,
                                                         msa->alnlen);
                        dm[i][j] = d;
                        dm[j][i] = d;
                }
        }

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

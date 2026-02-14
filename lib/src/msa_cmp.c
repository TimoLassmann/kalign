#include "tldevel.h"

#include <ctype.h>
#include "msa_struct.h"
#include "msa_check.h"
#include "msa_op.h"

#define MSA_CMP_IMPORT
#include "msa_cmp.h"

struct cmp_stats {
        uint64_t ref_total_aligned_pairs;
        uint64_t ref_total_gap_pairs;
        uint64_t identical_aligned;
        uint64_t identical_gaps;


        uint64_t test_total_aligned_pairs;
        uint64_t test_total_gap_pairs;
};

struct detailed_pair_stats {
        int64_t ref_scored_pairs;
        int64_t test_pairs;
        int64_t common_scored; /* matches in scored columns → for recall */
        int64_t common_all;    /* all matches → for precision */
};


static int compare_pair(char *seq1A, char *seq2A, char *seq1B, char *seq2B,
                        int len_a, int len_b, struct cmp_stats *stat);

static int compare_pair_detailed(char *seq1A, char *seq2A, char *seq1B, char *seq2B,
                                 int len_a, int len_b, int *scored_cols,
                                 struct detailed_pair_stats *stat);

int kalign_msa_compare(struct msa *r, struct msa *t,  float *score)
{
        struct cmp_stats* stat = NULL;
        ASSERT(r != NULL, "No reference alignment");
        ASSERT(t != NULL, "No test alignment");



        if(r->aligned == ALN_STATUS_ALIGNED){
                finalise_alignment(r);
        }

        if(t->aligned == ALN_STATUS_ALIGNED){
                finalise_alignment(t);
        }
        RUN(kalign_check_msa(r,1));
        RUN(kalign_check_msa(t,1));

        kalign_sort_msa(r);
        kalign_sort_msa(t);

        MMALLOC(stat, sizeof(struct cmp_stats));
        stat->identical_gaps = 0;

        stat->identical_aligned = 0;
        stat->ref_total_gap_pairs = 0;
        stat->ref_total_aligned_pairs = 0;

        stat->test_total_gap_pairs = 0;
        stat->test_total_aligned_pairs = 0;
        /* float s = 0.0; */
        /* float c = 0.0; */

        for(int i = 0; i < r->numseq;i++){
                for(int j = i + 1; j < r->numseq;j++){
                        compare_pair(r->sequences[i]->seq,
                                     r->sequences[j]->seq,
                                     t->sequences[i]->seq,
                                     t->sequences[j]->seq,
                                     r->alnlen,
                                     t->alnlen,
                                     stat
                                );
                }
        }


        /* LOG_MSG("%ld %ld %ld %ld %ld %ld", */
        /*         stat->ref_total_aligned_pairs, */
        /*         stat->ref_total_gap_pairs, */
        /*         stat->test_total_aligned_pairs, */
        /*         stat->test_total_gap_pairs, */
        /*         stat->identical_aligned, */
        /*         stat->identical_gaps); */
        double a;
        double b;

        /* Pairs of aligned (and unaligned) residues in common with reference
         divided by number of pairs in reference */
        a = (double) (stat->identical_aligned+ stat->identical_gaps );
        b = (double) (stat->ref_total_aligned_pairs + stat->ref_total_gap_pairs);
        /* LOG_MSG("Score: %f (%f / %f)", 100.0 * a  / b, a,b); */

        /* Pairs of aligned residues in common with reference
         divided by number of pairs in reference */
        a = (double) (stat->identical_aligned);
        b = (double) (stat->ref_total_aligned_pairs);
        /* LOG_MSG("Score: %f (%f / %f)", 100.0 * a  / b, a,b); */

        /* Pairs of aligned (and unaligned) residues in common with reference ,
           divided by average number if pairs  */
        a = (double) (stat->identical_aligned+ stat->identical_gaps );
        b = (double)(stat->ref_total_aligned_pairs + stat->ref_total_gap_pairs + stat->test_total_aligned_pairs + stat->test_total_gap_pairs)/ 2.0;
        /* LOG_MSG("Score: %f (%f / %f)", 100.0 * a  / b, a,b); */

        /* Pairs of aligned residues in common with reference ,
           divided by average number if pairs  */
        a = (double) (stat->identical_aligned );
        b = (double)(stat->ref_total_aligned_pairs   + stat->test_total_aligned_pairs) / 2.0;
        /* LOG_MSG("Score: %f (%f / %f)", 100.0 * a  / b, a,b); */




        /* LOG_MSG("ORG Score: %f", s); */
        a = (double) (stat->identical_aligned+ stat->identical_gaps );
        b = (double) (stat->ref_total_aligned_pairs + stat->ref_total_gap_pairs);
        *score = 100.0 * a  / b;
        /* exit(0); */
        MFREE(stat);
        return OK;
ERROR:
        return FAIL;
}

int compare_pair(char* seq1A, char* seq2A, char* seq1B, char* seq2B, int len_a, int len_b,struct cmp_stats* stat)
{
        int* codes1_A = NULL;
        int* codes2_A = NULL;

        int* codes1_B = NULL;
        int* codes2_B = NULL;

        int p1;
        int p2;
        int g1;
        int g2;

        if(len_a == 0 || len_b == 0){
                return OK;
        }

        MMALLOC(codes1_A, sizeof(int) * len_a);
        MMALLOC(codes2_A, sizeof(int) * len_a);

        MMALLOC(codes1_B, sizeof(int) * len_b);
        MMALLOC(codes2_B, sizeof(int) * len_b);
        /* process aln A */
        p1 = -1;
        p2 = -1;
        for(int i = 0; i < len_a;i++){
                g1 = 0;
                g2 = 0;
                if(isalpha((int) seq1A[i])){
                        p1++;
                }else{
                        g1 = 1;
                }
                if(isalpha((int) seq2A[i])){
                        p2++;
                }else{
                        g2 = 1;
                }
                if(!g1 && ! g2){
                        stat->ref_total_aligned_pairs++;
                        codes1_A[p1] = p2;
                        stat->ref_total_aligned_pairs++;
                        codes2_A[p2] = p1;
                }else if(!g1 && g2){
                        stat->ref_total_gap_pairs++;
                        codes1_A[p1] = -1;
                }else if(g1 && !g2){
                        stat->ref_total_gap_pairs++;
                        codes2_A[p2] = -1;
                }
        }

        /* process aln B  */
        p1 = -1;
        p2 = -1;
        for(int i = 0; i < len_b;i++){
                g1 = 0;
                g2 = 0;
                if(isalpha((int) seq1B[i])){
                        p1++;
                }else{
                        g1 = 1;
                }
                if(isalpha((int) seq2B[i])){
                        p2++;
                }else{
                        g2 = 1;
                }
                if(!g1 && ! g2){
                        stat->test_total_aligned_pairs++;
                        codes1_B[p1] = p2;
                        stat->test_total_aligned_pairs++;
                        codes2_B[p2] = p1;
                }else if(!g1 && g2){
                         stat->test_total_gap_pairs++;
                        codes1_B[p1] = -1;
                }else if(g1 && !g2){
                        stat->test_total_gap_pairs++;
                        codes2_B[p2] = -1;
                }
        }
        /* LOG_MSG("P1: %d P2: %d len b: %d", p1, p2,len_b); */
        /* fprintf(stdout,"%s\n%s\n",seq1A,seq2A); */
        /* fprintf(stdout,"%s\n%s\n",seq1B, seq2B); */
        /* fprintf(stdout,"\n"); */

        for(int i = 0; i <= p1;i++){
                if(codes1_A[i] != -1){
                        if(codes1_A[i] == codes1_B[i]){
                                stat->identical_aligned++;
                        }
                }else{
                        if(codes1_A[i] == codes1_B[i]){
                                stat->identical_gaps++;
                        }
                }
                /* if(codes1_A[i] != codes1_B[i]){ */
                /*         d += 1.0; */
                /*         /\* LOG_MSG("%d %d %d", i, codes1_A[i], codes1_B[i]); *\/ */
                /* } */
        }
        for(int i = 0; i <= p2;i++){
                if(codes2_A[i] != -1){
                        if(codes2_A[i] == codes2_B[i]){
                                stat->identical_aligned++;
                        }
                }else{
                        if(codes2_A[i] == codes2_B[i]){
                                stat->identical_gaps++;
                        }
                }


                /* if(codes2_A[i] != codes2_B[i]){ */
                /*         d += 1.0; */
                /*         /\* LOG_MSG("%d %d %d (second seq)", i, codes1_A[i], codes1_B[i]); *\/ */
                /* } */
        }

        /* s = d /(float)( p1 + p2); */
        /* LOG_MSG("score: %f", s); */
        /* *score = s; */
        MFREE(codes1_A);
        MFREE(codes1_B);
        MFREE(codes2_A);
        MFREE(codes2_B);
        return OK;
ERROR:
        MFREE(codes1_A);
        MFREE(codes1_B);
        MFREE(codes2_A);
        MFREE(codes2_B);
        return FAIL;
}


/* Shared helper: compute POAR + TC scores given a pre-built column mask. */
static int compare_with_mask_helper(struct msa *r, struct msa *t,
                                    int *scored_cols, struct poar_score *out)
{
        struct detailed_pair_stats dstat;
        int tc_correct = 0;
        int tc_total = 0;

        /* Accumulate POAR stats over all sequence pairs */
        dstat.ref_scored_pairs = 0;
        dstat.test_pairs = 0;
        dstat.common_scored = 0;
        dstat.common_all = 0;

        for(int i = 0; i < r->numseq; i++){
                for(int j = i + 1; j < r->numseq; j++){
                        compare_pair_detailed(
                                r->sequences[i]->seq,
                                r->sequences[j]->seq,
                                t->sequences[i]->seq,
                                t->sequences[j]->seq,
                                r->alnlen,
                                t->alnlen,
                                scored_cols,
                                &dstat
                        );
                }
        }

        /* TC score: for each scored column in reference, check if all non-gap
           residues map to the same column in the test alignment.
           Build residue-to-test-column map for each sequence. */
        {
                int** res_to_tcol = NULL;
                MMALLOC(res_to_tcol, sizeof(int*) * t->numseq);
                for(int s = 0; s < t->numseq; s++){
                        res_to_tcol[s] = NULL;
                        MMALLOC(res_to_tcol[s], sizeof(int) * (t->sequences[s]->len + 1));
                }

                /* Fill residue→test_column map */
                for(int s = 0; s < t->numseq; s++){
                        int p = 0;
                        for(int c = 0; c < t->alnlen; c++){
                                if(isalpha((int)t->sequences[s]->seq[c])){
                                        res_to_tcol[s][p] = c;
                                        p++;
                                }
                        }
                }

                /* Check each scored reference column */
                for(int c = 0; c < r->alnlen; c++){
                        if(!scored_cols[c]){
                                continue;
                        }
                        int first_tcol = -1;
                        int all_same = 1;
                        int nres = 0;
                        /* Count non-gap residues at this column */
                        for(int s = 0; s < r->numseq; s++){
                                if(isalpha((int)r->sequences[s]->seq[c])){
                                        nres++;
                                }
                        }
                        if(nres < 2){
                                /* Need at least 2 residues for a "pair" column */
                                continue;
                        }

                        tc_total++;
                        /* Check if all residues map to same test column */
                        {
                                int pos[r->numseq];
                                for(int s = 0; s < r->numseq; s++){
                                        pos[s] = 0;
                                }
                                for(int cc = 0; cc < c; cc++){
                                        for(int s = 0; s < r->numseq; s++){
                                                if(isalpha((int)r->sequences[s]->seq[cc])){
                                                        pos[s]++;
                                                }
                                        }
                                }
                                for(int s = 0; s < r->numseq; s++){
                                        if(isalpha((int)r->sequences[s]->seq[c])){
                                                int tcol = res_to_tcol[s][pos[s]];
                                                if(first_tcol < 0){
                                                        first_tcol = tcol;
                                                }else if(tcol != first_tcol){
                                                        all_same = 0;
                                                        break;
                                                }
                                        }
                                }
                        }
                        if(all_same){
                                tc_correct++;
                        }
                }

                for(int s = 0; s < t->numseq; s++){
                        MFREE(res_to_tcol[s]);
                }
                MFREE(res_to_tcol);
        }

        /* Fill output */
        out->ref_pairs = dstat.ref_scored_pairs;
        out->test_pairs = dstat.test_pairs;
        out->common = dstat.common_scored;

        if(dstat.ref_scored_pairs > 0){
                out->recall = (double)dstat.common_scored / (double)dstat.ref_scored_pairs;
        }else{
                out->recall = 0.0;
        }
        if(dstat.test_pairs > 0){
                out->precision = (double)dstat.common_all / (double)dstat.test_pairs;
        }else{
                out->precision = 0.0;
        }
        if(out->recall + out->precision > 0.0){
                out->f1 = 2.0 * out->recall * out->precision / (out->recall + out->precision);
        }else{
                out->f1 = 0.0;
        }
        if(tc_total > 0){
                out->tc = (double)tc_correct / (double)tc_total;
        }else{
                out->tc = 0.0;
        }

        return OK;
ERROR:
        return FAIL;
}

int kalign_msa_compare_detailed(struct msa *r, struct msa *t,
                                float max_gap_frac, struct poar_score *out)
{
        int* scored_cols = NULL;

        ASSERT(r != NULL, "No reference alignment");
        ASSERT(t != NULL, "No test alignment");
        ASSERT(out != NULL, "No output struct");

        if(r->aligned == ALN_STATUS_ALIGNED){
                finalise_alignment(r);
        }
        if(t->aligned == ALN_STATUS_ALIGNED){
                finalise_alignment(t);
        }
        RUN(kalign_check_msa(r, 1));
        RUN(kalign_check_msa(t, 1));

        kalign_sort_msa(r);
        kalign_sort_msa(t);

        /* Build scored column mask from reference alignment */
        MMALLOC(scored_cols, sizeof(int) * r->alnlen);
        for(int c = 0; c < r->alnlen; c++){
                if(max_gap_frac < 0.0f){
                        scored_cols[c] = 1;
                }else{
                        int ngaps = 0;
                        for(int s = 0; s < r->numseq; s++){
                                if(!isalpha((int)r->sequences[s]->seq[c])){
                                        ngaps++;
                                }
                        }
                        float gf = (float)ngaps / (float)r->numseq;
                        scored_cols[c] = (gf <= max_gap_frac) ? 1 : 0;
                }
        }

        RUN(compare_with_mask_helper(r, t, scored_cols, out));

        MFREE(scored_cols);
        return OK;
ERROR:
        MFREE(scored_cols);
        return FAIL;
}

int kalign_msa_compare_with_mask(struct msa *r, struct msa *t,
                                 int *scored_cols, int n_cols,
                                 struct poar_score *out)
{
        ASSERT(r != NULL, "No reference alignment");
        ASSERT(t != NULL, "No test alignment");
        ASSERT(scored_cols != NULL, "No column mask");
        ASSERT(out != NULL, "No output struct");

        if(r->aligned == ALN_STATUS_ALIGNED){
                finalise_alignment(r);
        }
        if(t->aligned == ALN_STATUS_ALIGNED){
                finalise_alignment(t);
        }
        RUN(kalign_check_msa(r, 1));
        RUN(kalign_check_msa(t, 1));

        kalign_sort_msa(r);
        kalign_sort_msa(t);

        ASSERT(n_cols == r->alnlen,
               "Mask length (%d) != reference alignment length (%d)",
               n_cols, r->alnlen);

        RUN(compare_with_mask_helper(r, t, scored_cols, out));

        return OK;
ERROR:
        return FAIL;
}


int compare_pair_detailed(char* seq1A, char* seq2A, char* seq1B, char* seq2B,
                          int len_a, int len_b, int* scored_cols,
                          struct detailed_pair_stats* stat)
{
        int* codes1_A = NULL;
        int* codes2_A = NULL;
        int* codes1_B = NULL;
        int* codes2_B = NULL;
        int* in_scored1 = NULL;
        int* in_scored2 = NULL;

        int p1, p2, g1, g2;

        if(len_a == 0 || len_b == 0){
                return OK;
        }

        MMALLOC(codes1_A, sizeof(int) * len_a);
        MMALLOC(codes2_A, sizeof(int) * len_a);
        MMALLOC(codes1_B, sizeof(int) * len_b);
        MMALLOC(codes2_B, sizeof(int) * len_b);
        MMALLOC(in_scored1, sizeof(int) * len_a);
        MMALLOC(in_scored2, sizeof(int) * len_a);

        /* Process reference alignment */
        p1 = -1;
        p2 = -1;
        for(int i = 0; i < len_a; i++){
                g1 = 0;
                g2 = 0;
                if(isalpha((int)seq1A[i])){
                        p1++;
                        in_scored1[p1] = 0;
                }else{
                        g1 = 1;
                }
                if(isalpha((int)seq2A[i])){
                        p2++;
                        in_scored2[p2] = 0;
                }else{
                        g2 = 1;
                }
                if(!g1 && !g2){
                        /* Both residues aligned */
                        codes1_A[p1] = p2;
                        codes2_A[p2] = p1;
                        if(scored_cols[i]){
                                stat->ref_scored_pairs += 2; /* both directions */
                                in_scored1[p1] = 1;
                                in_scored2[p2] = 1;
                        }
                }else if(!g1 && g2){
                        codes1_A[p1] = -1;
                }else if(g1 && !g2){
                        codes2_A[p2] = -1;
                }
        }

        /* Process test alignment */
        p1 = -1;
        p2 = -1;
        for(int i = 0; i < len_b; i++){
                g1 = 0;
                g2 = 0;
                if(isalpha((int)seq1B[i])){
                        p1++;
                }else{
                        g1 = 1;
                }
                if(isalpha((int)seq2B[i])){
                        p2++;
                }else{
                        g2 = 1;
                }
                if(!g1 && !g2){
                        stat->test_pairs += 2; /* both directions */
                        codes1_B[p1] = p2;
                        codes2_B[p2] = p1;
                }else if(!g1 && g2){
                        codes1_B[p1] = -1;
                }else if(g1 && !g2){
                        codes2_B[p2] = -1;
                }
        }

        /* Compare: count scored matches (for recall) and all matches (for precision) */
        for(int i = 0; i <= p1; i++){
                if(codes1_A[i] >= 0 && codes1_A[i] == codes1_B[i]){
                        stat->common_all++;
                        if(in_scored1[i]){
                                stat->common_scored++;
                        }
                }
        }
        for(int i = 0; i <= p2; i++){
                if(codes2_A[i] >= 0 && codes2_A[i] == codes2_B[i]){
                        stat->common_all++;
                        if(in_scored2[i]){
                                stat->common_scored++;
                        }
                }
        }

        MFREE(codes1_A);
        MFREE(codes1_B);
        MFREE(codes2_A);
        MFREE(codes2_B);
        MFREE(in_scored1);
        MFREE(in_scored2);
        return OK;
ERROR:
        MFREE(codes1_A);
        MFREE(codes1_B);
        MFREE(codes2_A);
        MFREE(codes2_B);
        MFREE(in_scored1);
        MFREE(in_scored2);
        return FAIL;
}

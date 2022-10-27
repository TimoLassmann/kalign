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



static int compare_pair(char *seq1A, char *seq2A, char *seq1B, char *seq2B,
                        int len_a, int len_b, struct cmp_stats *stat);

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

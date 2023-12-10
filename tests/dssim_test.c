#include "tldevel.h"
#include "tlrng.h"

#include "kalign/kalign.h"

#include "msa_struct.h"
#include "msa_sort.h"
#include "msa_op.h"
#include "dssim.h"

int msa_restore_ord_seq_order(struct msa* m);
int msa_remove_gaps(struct msa *m);
/* int test_consistency(int num_tests, int numseq, int seed); */
int test_consistency(int num_tests, int numseq,int dna,int seed);

void progress_bar(int c, int t, int w);

int main(void)
{
        LOG_MSG("DSSim - a simplistic sequence simulator.");

        LOG_MSG("Protein alignments - 20 seq");
        RUN(test_consistency(1000,20,0,42));
        LOG_MSG("Protein alignments - 1000 seq");
        RUN(test_consistency(10,1000,0,42));
        LOG_MSG("Protein alignments - 10000 seq");
        RUN(test_consistency(2,2000,0,42));

        LOG_MSG("DNA     alignments - 20 seq");
        RUN(test_consistency(1000,20,1,42));
        LOG_MSG("DNA     alignments - 1000 seq");
        RUN(test_consistency(10,1000,1,42));
        LOG_MSG("DNA     alignments - 10000 seq");
        RUN(test_consistency(2,2000,1,42));

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int test_consistency(int num_tests, int numseq,int dna,int seed)
{
        struct msa* m = NULL;
        struct msa* m2 = NULL;
        struct rng_state* rng = NULL;


        RUNP(rng = init_rng(seed));

        for(int i = 0; i < num_tests;i++){
                int local_seed = tl_random_int(rng, 1000000);
                int t1  = tl_random_int(rng, 8) + 1;
                int t2 = tl_random_int(rng, 8) + 1;
                while(t1 == t2){
                        t2 = tl_random_int(rng, 8) + 1;
                }
                float score = 0.0;

                dssim_get_fasta(&m, numseq, 10,dna, 45, local_seed);

                m->quiet = 1;

                msa_cpy(&m2, m);
                msa_shuffle_seq(m, rng);
                msa_shuffle_seq(m2, rng);
                kalign_run(m, t1, KALIGN_TYPE_UNDEFINED,0.0,0.0,0.0);
                kalign_run(m2, t2, KALIGN_TYPE_UNDEFINED,0.0,0.0,0.0);

                kalign_msa_compare(m, m2, &score);
                if(score != 100.0f){
                        LOG_MSG("Testing %d : %d %d %f", i , t1 ,t2,score);
                        kalign_write_msa(m, NULL, "msf");
                }
                kalign_free_msa(m);
                kalign_free_msa(m2);
                m = NULL;
                m2 = NULL;
                progress_bar(i, num_tests, 50);
        }
        progress_bar(num_tests, num_tests, 50);
        fprintf(stdout,"\n");
        free_rng(rng);
        return OK;
ERROR:
        return FAIL;
}

int msa_remove_gaps(struct msa *m)
{
        int aln_len = m->alnlen;
        for(int i = 0; i < m->numseq;i++){
                struct msa_seq* s = m->sequences[i];

                int c = 0;
                for(int j = 0; j < aln_len;j++){
                        if(s->seq[j] != '-'){
                                s->seq[c] = s->seq[j];
                                c++;
                        }
                }
                s->len = c;
                m->sequences[i]->rank = atoi(m->sequences[i]->name);
        }
        dealign_msa(m);
        return OK;
}


int msa_restore_ord_seq_order(struct msa* m)
{
        for(int i = 0; i < m->numseq;i++){
                m->sequences[i]->rank = atoi(m->sequences[i]->name);

        }
        RUN(msa_sort_rank(m));
        return OK;
ERROR:
        return FAIL;
}

void progress_bar(int c, int t, int w)

{
        double progress = (double) c / (double) t;
        int l = (int)round(progress * (double) w);

        fprintf(stdout,"\r[");
        for(int i = 0; i < w;i++){
                if(i < l){
                        fprintf(stdout,"=");
                }else{
                        fprintf(stdout," ");
                }
        }
        fprintf(stdout,"] %3d%%", (int)round(progress * 100.0));
        fflush(stdout);
}

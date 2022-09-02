/* #include "tldevel.h" */
#include "libkalign.h"

/* tldevel  */

#include "mod_tldevel.c"
/* Interface  */

#include "mod_interface.c"
/* Alignment tasks  */

#include "task.c"

/* MSa  */

#include "msa.c"

/* MSA io */

#include "mod_msaio.c"

/* Tree */

#include "mod_tree.c"
/* Alphabet  */
/* #include "alphabet.h" */
/* #include "alphabet.c" */
/* Alignment module  */

#include "mod_aln.c"

#include "alphabet.c"
/*  */

/* #include "alphabet.h" */
/* #include "global.h" */
/* #include "msa.h" */

/* #include "aln_param.h" */
/* #include "aln_run.h" */
/* #include "task.h" */


/* #include "bisectingKmeans.h" */

/* #include "weave_alignment.h" */


#ifdef HAVE_OPENMP
#include <omp.h>
#endif

static int to_msa(char** input_sequences, int* len, int numseq,struct msa** multiple_aln);
int echo_aln(struct msa* msa);

int kalign(char **seq, int *len,int numseq, char ***aligned, int *out_aln_len)
{
        struct msa* msa = NULL;
        struct aln_param* ap = NULL;
        struct aln_tasks* tasks = NULL;
        struct parameters* param = NULL;
        char** out = NULL;
        RUNP(param = init_param());

        ASSERT(seq != NULL, "No sequences");
        /* get sequences into aln struct  */

        RUN(to_msa(seq, len, numseq,&msa));

        RUN(detect_alphabet(msa));
        RUN(detect_aligned(msa));

        RUN(set_sip_nsip(msa));
        if(!msa->quiet){
                LOG_MSG("Detected: %d sequences.", msa->numseq);
        }
        if(msa->aligned != ALN_STATUS_UNALIGNED){
                RUN(dealign_msa(msa));
        }
        RUN(init_ap(&ap,param,msa->L ));

        RUN(alloc_tasks(&tasks, msa->numseq));
        /* msa =  alloc_msa(); */
        /* test for dna / proten  */
        /* align  */

        RUN(build_tree_kmeans(msa,ap->nthreads,1,&tasks));

        /* by default all protein sequences are converted into a reduced alphabet
           when read from file. Here we turn them back into the default representation. */
        if(msa->L == ALPHA_redPROTEIN){
                //RUN(convert_msa_to_internal(msa, defPROTEIN));
                RUN(convert_msa_to_internal(msa, ALPHA_ambigiousPROTEIN));
        }
        /* allocate aln parameters  */
        RUN(init_ap(&ap,param,msa->L ));
        /* align */
        RUN(create_msa_tree(msa, ap, tasks,param->nthreads));
        msa->aligned = ALN_STATUS_ALIGNED;


        /* echo_aln(msa); */

        int aln_len = 0;
        for (int j = 0; j <= msa->sequences[0]->len;j++){
                aln_len += msa->sequences[0]->gaps[j];
        }
        aln_len += msa->sequences[0]->len;
        aln_len += 1;

        MMALLOC(out, sizeof(char*) * numseq);
        for(int i = 0; i < numseq; i++){
                out[i] = 0;
                MMALLOC(out[i], sizeof(char) * (uint64_t)aln_len);
                int pos = 0;
                for(int j = 0;j < msa->sequences[i]->len;j++){
                        for(int c = 0;c < msa->sequences[i]->gaps[j];c++){
                                out[i][pos] = '-';
                                pos++;
                        }
                        out[i][pos] = msa->sequences[i]->seq[j];
                        pos++;
                }
                for(int c = 0;c < msa->sequences[i]->gaps[ msa->sequences[i]->len];c++){
                        out[i][pos] = '-';
                        pos++;
                }
                out[i][pos] = 0;
                /* LOG_MSG("%d %d %s ", aln_len,pos, out[i]); */
        }

        *aligned = out;
        *out_aln_len = aln_len;
        /* clean up map */
        free_tasks(tasks);
        free_ap(ap);
        /* copy alignment  */
        free_parameters(param);
        free_msa(msa);
        return OK;
ERROR:
        free_parameters(param);
        free_msa(msa);
        free_ap(ap);
        return FAIL;
}

int echo_aln(struct msa* msa)
{
        FILE* f_ptr = stdout;
        ASSERT(msa != NULL,"No alignment");
        for(int i = 0; i < msa->numseq;i++){
                fprintf(f_ptr,">%s\n", msa->sequences[i]->name);
                int f = 0;
                for(int j = 0;j < msa->sequences[i]->len;j++){
                        for(int c = 0;c < msa->sequences[i]->gaps[j];c++){
                                fprintf(f_ptr,"-");
                                f++;
                                if(f == 60){
                                        fprintf(f_ptr, "\n");
                                        f = 0;
                                }
                        }
                        fprintf(f_ptr,"%c", msa->sequences[i]->seq[j]);
                        f++;
                        if(f == 60){
                                fprintf(f_ptr, "\n");
                                f = 0;
                        }
                }
                for(int c = 0;c < msa->sequences[i]->gaps[ msa->sequences[i]->len];c++){
                        fprintf(f_ptr,"-");
                        f++;
                        if(f == 60){
                                fprintf(f_ptr, "\n");
                                f = 0;
                        }
                }
                if(f){
                        fprintf(f_ptr,"\n");
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int to_msa(char** input_sequences, int* len, int numseq,struct msa** multiple_aln)
{
        struct msa* msa = NULL;

        MMALLOC(msa, sizeof(struct msa));
        msa->sequences = NULL;
        msa->alloc_numseq = numseq;
        msa->numseq = numseq;
        msa->num_profiles = 0;
        msa->L = ALPHA_UNDEFINED;
        msa->aligned = 0;
        msa->plen = NULL;
        msa->sip = NULL;
        msa->nsip = NULL;
        msa->quiet = 1;
        MMALLOC(msa->sequences, sizeof(struct msa_seq*) * msa->alloc_numseq);

        for(int i = 0; i < 128; i++){
                msa->letter_freq[i] = 0;
        }

        for(int i = 0; i < msa->alloc_numseq;i++){
                msa->sequences[i] = NULL;
                struct msa_seq* seq = NULL;

                MMALLOC(seq, sizeof(struct msa_seq));
                seq->name = NULL;
                seq->seq = NULL;
                seq->s = NULL;
                seq->gaps = NULL;
                seq->len = len[i];
                seq->alloc_len = len[i]+1;

                MMALLOC(seq->name, sizeof(char)* MSA_NAME_LEN);

                MMALLOC(seq->seq, sizeof(char) * seq->alloc_len);
                MMALLOC(seq->s, sizeof(uint8_t) * seq->alloc_len);
                MMALLOC(seq->gaps, sizeof(int) * (seq->alloc_len+1));
                for(int j = 0;j < seq->alloc_len+1;j++){

                        seq->gaps[j] = 0;

                }
                for(int j = 0; j < len[i];j++){
                        msa->letter_freq[(int)input_sequences[i][j]]++;
                        seq->seq[j] = input_sequences[i][j];
                }
                seq->seq[len[i]] = 0;
                msa->sequences[i] = seq;
                /* LOG_MSG("%s",msa->sequences[i]->seq); */
        }

        *multiple_aln = msa;
        return OK;

ERROR:
        free_msa(msa);
        return FAIL;
}


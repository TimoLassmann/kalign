#include "tldevel.h"
#include "esl_stopwatch.h"
#include "task.h"
#include "msa_struct.h"
#include "msa_op.h"
#include "msa_alloc.h"
#include "alphabet.h"
#include "bisectingKmeans.h"

#include "aln_param.h"
#include "aln_run.h"


#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#define ALN_WRAP_IMPORT
#include "aln_wrap.h"


int kalign(char **seq, int *len, int numseq,int n_threads, int type, float gpo, float gpe, float tgpe, char ***aligned, int *out_aln_len)
{
        struct msa *msa;
        RUN(kalign_arr_to_msa(seq, len,numseq, &msa));
        msa->quiet = 1;
        if(n_threads < 1){
                n_threads = 1;
        }
        RUN(kalign_run(msa,n_threads, type,  gpo, gpe, tgpe));
        RUN(kalign_msa_to_arr(msa, aligned, out_aln_len));
        kalign_free_msa(msa);

        return OK;
ERROR:
        if(msa){
                kalign_free_msa(msa);
        }
        return FAIL;
}

int kalign_run(struct msa *msa, int n_threads, int type, float gpo, float gpe, float tgpe)
{
        struct aln_tasks* tasks = NULL;
        struct aln_param* ap = NULL;

        /* If already aligned unalign ! */
        if(msa->aligned != ALN_STATUS_UNALIGNED){
                RUN(dealign_msa(msa));
        }

        /* Convert into internal representation  */
        if(msa->biotype == ALN_BIOTYPE_DNA){
                msa->L = ALPHA_defDNA;
                RUN(convert_msa_to_internal(msa, ALPHA_defDNA));
        }else if(msa->biotype == ALN_BIOTYPE_PROTEIN){
                msa->L = ALPHA_redPROTEIN;
                RUN(convert_msa_to_internal(msa, ALPHA_redPROTEIN));
        }else{
                ERROR_MSG("Unable to determine what alphabet to use.");
        }
        /* -LOG_MSG("L: %d  threads: %d",msa->L, n_threads); */
        /* Start the heavy lifting  */
        RUN(alloc_tasks(&tasks, msa->numseq));

#ifdef HAVE_OPENMP
        omp_set_num_threads(n_threads);
#endif
        /* Build guide tree */
        RUN(build_tree_kmeans(msa,n_threads,&tasks));


        /* Convert to full alphabet after having converted to reduced alphabet for tree building above  */
        if(msa->biotype == ALN_BIOTYPE_PROTEIN){
                RUN(convert_msa_to_internal(msa, ALPHA_ambigiousPROTEIN));
        }
        /* LOG_MSG("L: %d",msa->L); */
        /* align  */
        /* if gap penalties are not negative they are set below */
        RUN(aln_param_init(&ap,
                           msa->biotype,
                           n_threads,
                           type,
                           gpo,
                           gpe,
                           tgpe));

        DECLARE_TIMER(t1);
        if(!msa->quiet){
                LOG_MSG("Aligning");
        }
        START_TIMER(t1);

        RUN(create_msa_tree(msa, ap, tasks));
        /* Hurrah we have aligned sequences  */
        msa->aligned = ALN_STATUS_ALIGNED;

        RUN(finalise_alignment(msa));

        STOP_TIMER(t1);
        if(!msa->quiet){
                GET_TIMING(t1);
        }
        DESTROY_TIMER(t1);

        aln_param_free(ap);
        free_tasks(tasks);
        return OK;
ERROR:
        aln_param_free(ap);
        free_tasks(tasks);
        return FAIL;
}

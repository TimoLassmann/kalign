#include "tldevel.h"

#include "msa_struct.h"

#define ALN_PARAM_IMPORT
#include "aln_param.h"

static int set_subm_gaps_CorBLOSUM66_13plus(struct aln_param *ap);
static int set_subm_gaps_DNA(struct aln_param *ap);
static int set_subm_gaps_DNA_internal(struct aln_param *ap);
static int set_subm_gaps_RNA(struct aln_param *ap);

int aln_param_init(struct aln_param **aln_param,int biotype , int n_threads, int type, float gpo, float gpe, float tgpe)
{
        struct aln_param* ap = NULL;

        /* Allocate  */
        MMALLOC(ap, sizeof(struct aln_param));
        ap->subm = NULL;

        ap->nthreads = n_threads;
        MMALLOC(ap->subm,sizeof (float*) * 23);

        for (int i = 23;i--;){
                ap->subm[i] = NULL;
                MMALLOC(ap->subm[i],sizeof(float) * 23);
                for (int j = 23;j--;){
                        ap->subm[i][j] = 0.0f;
                }
        }
        if(biotype == ALN_BIOTYPE_DNA){
                /* include/kalign/ */
                if(type == KALIGN_DNA){
                        set_subm_gaps_DNA(ap);
                }else if(type == KALIGN_DNA_INTERNAL){
                        set_subm_gaps_DNA_internal(ap);
                }else if(type == KALIGN_RNA){
                        set_subm_gaps_RNA(ap);
                }else{
                        ERROR_MSG("DNA alignment type (%d) not recognised.");
                }

        }else if(biotype == ALN_BIOTYPE_PROTEIN){
                set_subm_gaps_CorBLOSUM66_13plus(ap);
                /* set_subm_gaps(ap); */
        }else{
                ERROR_MSG("Unable to determine what alphabet to use.");
        }

        if(gpo >= 0.0){
                ap->gpo = gpo;
        }
        if(gpe >= 0.0){
                ap->gpe = gpe;
        }

        if(gpe >= 0.0){
                ap->tgpe = tgpe;
        }
        /* LOG_MSG("%f %f %f", ap->gpo, ap->gpe, ap->tgpe); */
        *aln_param = ap;
        return OK;
ERROR:
        aln_param_free(ap);
        return FAIL;
}


int set_subm_gaps_CorBLOSUM66_13plus(struct aln_param* ap)
{
        int i;
        int j;
        //char aacode[20] = "ACDEFGHIKLMNPQRSTVWY";
        //char aa_order[23] = "ARNDCQEGHILKMFPSTWYVBZX";
        //int num_aa = 23;

        int CorBLOSUM66_13plus[23][23] = {
                /*A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X*/
                {5,-1,-1,-2,-2,-1,-1,0,-2,-1,-1,-1,0,-2,-1,1,0,-2,-2,0,-2,-1,0},
                {-1,6,0,-1,-3,1,1,-2,0,-2,-2,3,-1,-3,-1,-1,-1,-1,-1,-2,0,1,-1},
                {-1,0,6,2,-3,1,0,0,0,-3,-3,0,-2,-2,-1,1,0,-2,-1,-2,4,0,-1},
                {-2,-1,2,7,-3,1,2,-1,-1,-3,-3,0,-3,-3,-1,0,-1,-3,-2,-3,5,2,-1},
                {-2,-3,-3,-3,12,-3,-4,-3,-2,-2,-3,-3,-2,-1,-3,-2,-2,-3,-2,-2,-3,-3,-2},
                {-1,1,1,1,-3,5,2,-2,0,-2,-2,1,0,-2,-1,0,0,-1,-1,-2,1,3,0},
                {-1,1,0,2,-4,2,6,-2,-1,-3,-3,1,-2,-3,0,0,-1,-2,-2,-2,1,4,-1},
                {0,-2,0,-1,-3,-2,-2,7,-2,-4,-4,-2,-3,-3,-2,0,-2,-3,-3,-3,-1,-2,-1},
                {-2,0,0,-1,-2,0,-1,-2,10,-3,-3,0,-2,-2,-2,-1,-1,-2,1,-3,0,0,-1},
                {-1,-2,-3,-3,-2,-2,-3,-4,-3,5,2,-2,2,0,-3,-2,-1,-1,-1,3,-3,-2,-1},
                {-1,-2,-3,-3,-3,-2,-3,-4,-3,2,5,-2,3,1,-3,-3,-2,0,-1,1,-3,-2,-1},
                {-1,3,0,0,-3,1,1,-2,0,-2,-2,5,-1,-3,-1,0,0,-2,-2,-2,0,1,-1},
                {0,-1,-2,-3,-2,0,-2,-3,-2,2,3,-1,6,1,-2,-1,-1,0,-1,1,-2,-1,0},
                {-2,-3,-2,-3,-1,-2,-3,-3,-2,0,1,-3,1,7,-3,-2,-2,2,3,0,-3,-3,-1},
                {-1,-1,-1,-1,-3,-1,0,-2,-2,-3,-3,-1,-2,-3,9,0,-1,-2,-2,-2,-1,-1,-1},
                {1,-1,1,0,-2,0,0,0,-1,-2,-3,0,-1,-2,0,4,2,-2,-2,-1,0,0,0},
                {0,-1,0,-1,-2,0,-1,-2,-1,-1,-2,0,-1,-2,-1,2,5,-1,-1,0,0,0,0},
                {-2,-1,-2,-3,-3,-1,-2,-3,-2,-1,0,-2,0,2,-2,-2,-1,13,3,-2,-2,-2,-1},
                {-2,-1,-1,-2,-2,-1,-2,-3,1,-1,-1,-2,-1,3,-2,-2,-1,3,9,-1,-2,-2,-1},
                {0,-2,-2,-3,-2,-2,-2,-3,-3,3,1,-2,1,0,-2,-1,0,-2,-1,4,-3,-2,-1},
                {-2,0,4,5,-3,1,1,-1,0,-3,-3,0,-2,-3,-1,0,0,-2,-2,-3,4,1,-1},
                {-1,1,0,2,-3,3,4,-2,0,-2,-2,1,-1,-3,-1,0,0,-2,-2,-2,1,4,-1},
                {0,-1,-1,-1,-2,0,-1,-1,-1,-1,-1,-1,0,-1,-1,0,0,-1,-1,-1,-1,-1,-1},
        };

        for(i = 0; i < 23;i++){
                for(j = 0; j < 23;j++){
                        ap->subm[i][j] = (float)(CorBLOSUM66_13plus[i][j]);// *2.0F;
                }
        }
        ap->gpo = 5.5F;// * 2.0F;
        ap->gpe = 2.0F;// * 2.0F;
        ap->tgpe = 1.0F;// * 2.0F;
        return OK;
}

int set_subm_gaps_DNA(struct aln_param *ap)
{
        int i,j;
        for(i = 0; i < 5; i++){
                for(j =0; j < 5;j++){
                        ap->subm[i][j] = -4;
                        if(i == j){
                                ap->subm[i][j] = 5;
                        }
                }
        }
        ap->gpo = 8;
        ap->gpe = 6;
        ap->tgpe = 0;
        return OK;
}

int set_subm_gaps_DNA_internal(struct aln_param *ap)
{
        int i,j;
        for(i = 0; i < 5; i++){
                for(j =0; j < 5;j++){
                        ap->subm[i][j] = -4;
                        if(i == j){
                                ap->subm[i][j] = 5;
                        }
                }
        }
        ap->gpo = 8;
        ap->gpe = 6;
        ap->tgpe = 8;
        return OK;
}

int set_subm_gaps_RNA(struct aln_param* ap)
{
        int i,j;
        for(i = 0; i < 5; i++){
                for(j =0; j < 5;j++){
                        ap->subm[i][j] = 283.0;
                }
        }
//	A   91 -114  -31 -123    0  -43
        ap->subm[0][0] += 91.0;
        ap->subm[0][1] += -114.0;
        ap->subm[0][2] += -31.0;
        ap->subm[0][3] += -123.0;

//	C -114  100 -125  -31    0  -43
        ap->subm[1][0] += -114.0;
        ap->subm[1][1] += 100.0;
        ap->subm[1][2] += -125.0;
        ap->subm[1][3] += -31.0;

//	G  -31 -125  100 -114    0  -43
        ap->subm[2][0] += -31.0;
        ap->subm[2][1] += -125.0;
        ap->subm[2][2] += 100.0;
        ap->subm[2][3] += -114.0;

//	T -123  -31 -114   91    0  -43
        ap->subm[3][0] += -123.0;
        ap->subm[3][1] += -31.0;
        ap->subm[3][2] += -114.0;
        ap->subm[3][3] += 91.0;

        ap->gpo = 217.0;
        ap->gpe = 39.4;
        ap->tgpe = 292.6;
        return OK;
}

void aln_param_free(struct aln_param *ap)
{
        if(ap){
                if(ap->subm){
                        for (int i = 23;i--;){
                                MFREE(ap->subm[i]);
                        }
                        MFREE(ap->subm);
                }
                MFREE(ap);
        }
}

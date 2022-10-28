#include "tldevel.h"

#include "kalign/kalign.h"
#include "msa_struct.h"

#define ALN_PARAM_IMPORT
#include "aln_param.h"

static int set_subm_gaps_CorBLOSUM66_13plus(struct aln_param *ap);
static int set_subm_gaps_gon250(struct aln_param* ap);
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
                switch (type) {
                case KALIGN_TYPE_DNA:
                        set_subm_gaps_DNA(ap);
                        break;
                case KALIGN_TYPE_DNA_INTERNAL:
                        set_subm_gaps_DNA_internal(ap);
                        break;
                case KALIGN_TYPE_RNA:
                        set_subm_gaps_RNA(ap);
                        break;
                case KALIGN_TYPE_PROTEIN:
                        ERROR_MSG("Detected DNA sequences but --type protein option was selected.");
                        break;
                default:
                        set_subm_gaps_RNA(ap);
                        break;
                }
        }else if(biotype == ALN_BIOTYPE_PROTEIN){
                switch (type) {
                case KALIGN_TYPE_PROTEIN:
                        set_subm_gaps_CorBLOSUM66_13plus(ap);
                        break;
                case KALIGN_TYPE_PROTEIN_DIVERGENT:
                         set_subm_gaps_gon250(ap);
                        break;
                case KALIGN_TYPE_DNA:
                        ERROR_MSG("Detected protein sequences but --type dna option was selected.");
                        break;
                case KALIGN_TYPE_DNA_INTERNAL:
                        ERROR_MSG("Detected protein sequences but --type internal  option was selected.");
                        break;
                case KALIGN_TYPE_RNA:
                        ERROR_MSG("Detected protein sequences but --type rna option was selected.");
                        break;
                default:
                        set_subm_gaps_CorBLOSUM66_13plus(ap);
                        /* set_subm_gaps_gon250(ap); */
                        break;
                }
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



int set_subm_gaps_gon250(struct aln_param* ap)
{
        //A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,B,Z,X,
        int gon250mt[23][23] = {
                //  A,   R,   N,   D,   C,   Q,   E,   G,   H,   I,   L,   K,   M,   F,   P,   S,   T,   W,   Y,   V,   B,   Z,   X,
                {  24,  -6,  -3,  -3,   5,  -2,   0,   5,  -8,  -8, -12,  -4,  -7, -23,   3,  11,   6, -36, -22,   1,   0,   0,   0},// A
                {  -6,  47,   3,  -3, -22,  15,   4, -10,   6, -24, -22,  27, -17, -32,  -9,  -2,  -2, -16, -18, -20,   0,   0,   0},// R
                {  -3,   3,  38,  22, -18,   7,   9,   4,  12, -28, -30,   8, -22, -31,  -9,   9,   5, -36, -14, -22,   0,   0,   0},// N
                {  -3,  -3,  22,  47, -32,   9,  27,   1,   4, -38, -40,   5, -30, -45,  -7,   5,   0, -52, -28, -29,   0,   0,   0},// D
                {   5, -22, -18, -32, 115, -24, -30, -20, -13, -11, -15, -28,  -9,  -8, -31,   1,  -5, -10,  -5,   0,   0,   0,   0},// C
                {  -2,  15,   7,   9, -24,  27,  17, -10,  12, -19, -16,  15, -10, -26,  -2,   2,   0, -27, -17, -15,   0,   0,   0},// Q
                {   0,   4,   9,  27, -30,  17,  36,  -8,   4, -27, -28,  12, -20, -39,  -5,   2,  -1, -43, -27, -19,   0,   0,   0},// E
                {   5, -10,   4,   1, -20, -10,  -8,  66, -14, -45, -44, -11, -35, -52, -16,   4, -11, -40, -40, -33,   0,   0,   0},// G
                {  -8,   6,  12,   4, -13,  12,   4, -14,  60, -22, -19,   6, -13,  -1, -11,  -2,  -3,  -8,  22, -20,   0,   0,   0},// H
                {  -8, -24, -28, -38, -11, -19, -27, -45, -22,  40,  28, -21,  25,  10, -26, -18,  -6, -18,  -7,  31,   0,   0,   0},// I
                { -12, -22, -30, -40, -15, -16, -28, -44, -19,  28,  40, -21,  28,  20, -23, -21, -13,  -7,   0,  18,   0,   0,   0},// L
                {  -4,  27,   8,   5, -28,  15,  12, -11,   6, -21, -21,  32, -14, -33,  -6,   1,   1, -35, -21, -17,   0,   0,   0},// K
                {  -7, -17, -22, -30,  -9, -10, -20, -35, -13,  25,  28, -14,  43,  16, -24, -14,  -6, -10,  -2,  16,   0,   0,   0},// M
                { -23, -32, -31, -45,  -8, -26, -39, -52,  -1,  10,  20, -33,  16,  70, -38, -28, -22,  36,  51,   1,   0,   0,   0},// F
                {   3,  -9,  -9,  -7, -31,  -2,  -5, -16, -11, -26, -23,  -6, -24, -38,  76,   4,   1, -50, -31, -18,   0,   0,   0},// P
                {  11,  -2,   9,   5,   1,   2,   2,   4,  -2, -18, -21,   1, -14, -28,   4,  22,  15, -33, -19, -10,   0,   0,   0},// S
                {   6,  -2,   5,   0,  -5,   0,  -1, -11,  -3,  -6, -13,   1,  -6, -22,   1,  15,  25, -35, -19,   0,   0,   0,   0},// T
                { -36, -16, -36, -52, -10, -27, -43, -40,  -8, -18,  -7, -35, -10,  36, -50, -33, -35, 142,  41, -26,   0,   0,   0},// W
                { -22, -18, -14, -28,  -5, -17, -27, -40,  22,  -7,   0, -21,  -2,  51, -31, -19, -19,  41,  78, -11,   0,   0,   0},// Y
                {   1, -20, -22, -29,   0, -15, -19, -33, -20,  31,  18, -17,  16,   1, -18, -10,   0, -26, -11,  34,   0,   0,   0},// V
                {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},// B
                {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},// Z
                {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},// X
        };

        for(int i = 0; i < 23;i++){
                for(int j = 0; j < 23;j++){
                        ap->subm[i][j] = (float)(gon250mt[i][j]);// *2.0F;
                }
        }

        ap->gpo = 55;
        ap->gpe = 8;
        ap->tgpe = 4;
        return OK;
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

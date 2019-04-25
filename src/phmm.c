

#include "phmm.h"


#define INDEXMM 0
#define INDEXGPO 1
#define INDEXGPE 2
#define INDEXTM 3
/* need forward backward  */
/* or be adventurous and use beam sampling with u =1 i.e. no threshold ... worth a try?  */

struct phmm{

        float** fM ;
        float** fX ;
        float** fY ;

        float** bM ;
        float** bX ;
        float** bY ;



        float emit_M[26][26];
        float emit_background[26];

        float emit_M_e[26][26];
        float emit_background_e[26];

        float transition[4];
        float transition_e[4];

        float tau;
        float tau_e;
        float f_score;
        float b_score;

        int alloc_x;
        int alloc_y;
        int L;                  /* alphabet len */
};

int forward_phmm(struct phmm* phmm,int* seq_a,int* seq_b, int len_a,int len_b);
int backward_phmm(struct phmm* phmm,int* seq_a,int* seq_b, int len_a,int len_b);
int collect_phmm(struct phmm* phmm,int* seq_a,int* seq_b, int len_a,int len_b);
int re_estimate(struct phmm* phmm);

int phmm_transitions(struct phmm* phmm);
int print_phmm(struct phmm* phmm,int len_a,int len_b);

int clear_phmm_e(struct phmm* phmm);

int simple_init(struct phmm*phmm);
/* memory functions */
struct phmm* alloc_phmm(int size);

void free_phmm(struct phmm* phmm);

int main(int argc, char *argv[])
{
        struct phmm* phmm = NULL;
        int* seqa = NULL;
        int* seqb = NULL;

        int len_a = 5;
        int len_b = 5;

        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};


        /* example from biological sequence analysis */

        len_a = 7;
        len_b = 10;
        MMALLOC(seqa, sizeof(int) * len_a);
        MMALLOC(seqb, sizeof(int) * len_b);

        seqa[0] = aacode['P'-65];
        seqa[1] = aacode['A'-65];
        seqa[2] = aacode['W'-65];
        seqa[3] = aacode['H'-65];
        seqa[4] = aacode['E'-65];
        seqa[5] = aacode['A'-65];
        seqa[6] = aacode['E'-65];

        seqb[0] = aacode['H'-65];
        seqb[1] = aacode['E'-65];
        seqb[2] = aacode['A'-65];
        seqb[3] = aacode['G'-65];
        seqb[4] = aacode['A'-65];
        seqb[5] = aacode['W'-65];
        seqb[6] = aacode['G'-65];
        seqb[7] = aacode['H'-65];
        seqb[8] = aacode['E'-65];
        seqb[9] = aacode['E'-65];



        RUNP(phmm = alloc_phmm(MACRO_MAX(len_a,len_b)));
        RUN(simple_init(phmm));

        RUN(clear_phmm_e(phmm));

        //RUN(phmm_transitions(phmm));
        for(int iter = 0;iter < 10;iter++){
                RUN(forward_phmm(phmm, seqa, seqb, len_a, len_b));
                RUN(backward_phmm(phmm, seqa, seqb, len_a, len_b));


                RUN(collect_phmm(phmm, seqa, seqb, len_a, len_b));
                RUN(re_estimate(phmm));

                fprintf(stdout,"%f\tforward\n%f\tbackward\n",phmm->f_score, phmm->b_score);
          
        }
        RUN(phmm_transitions(phmm));
        //print_phmm(phmm, len_a, len_b);
        free_phmm(phmm);

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int forward_phmm(struct phmm* phmm,int* seq_a,int* seq_b, int len_a,int len_b)
{

        int i,j;
        int* sa = seq_a -1;
        int* sb = seq_b -1;
        const float MM = phmm->transition[INDEXMM];
        const float MX = phmm->transition[INDEXGPO];
        const float MY = phmm->transition[INDEXGPO];
        const float XX = phmm->transition[INDEXGPE];
        const float XM = phmm->transition[INDEXTM];
        const float YY = phmm->transition[INDEXGPE];
        const float YM = phmm->transition[INDEXTM];





        float** M = phmm->fM;
        float** X = phmm->fX;
        float** Y = phmm->fY;
/*const tAA = prob2scaledprob(1.0f - scaledprob2prob(phmm->eta));
        const tBB = prob2scaledprob(1.0f - scaledprob2prob(phmm->eta));
        const tAM = phmm->eta;
        const tBM = phmm->eta;*/

        /* init first element  */
        M[0][0] = prob2scaledprob(1.0f);
        X[0][0] = prob2scaledprob(0.0f);
        Y[0][0] = prob2scaledprob(0.0f);

        for(j = 1; j <= len_b;j++){
                M[0][j] = prob2scaledprob(0.0f);

                X[0][j] = prob2scaledprob(0.0f);

                Y[0][j] =                 M[0][j-1] + MY;
                Y[0][j] = logsum(Y[0][j], Y[0][j-1] + YY);
                Y[0][j] += phmm->emit_background[sb[j]];
        }

        for(i = 1; i <= len_a;i++){
                M[i][0] = prob2scaledprob(0.0f);

                X[i][0] =                 M[i-1][0] + MX;
                X[i][0] = logsum(X[i][0], X[i-1][0] + XX);
                X[i][0] += phmm->emit_background[sa[i]];

                Y[i][0] = prob2scaledprob(0.0f);

                for(j = 1; j<= len_b;j++){
                        M[i][j] = M[i-1][j-1] + MM;
                        M[i][j] = logsum(M[i][j], X[i-1][j-1] + XM);
                        M[i][j] = logsum(M[i][j], Y[i-1][j-1] + YM);
                        M[i][j] += phmm->emit_M[sa[i]][sb[j]];

                        X[i][j] =                 M[i-1][j] + MX;
                        X[i][j] = logsum(X[i][j], X[i-1][j] + XX);
                        X[i][j] += phmm->emit_background[sa[i]];

                        Y[i][j] =                 M[i][j-1] + MY;
                        Y[i][j] = logsum(Y[i][j], Y[i][j-1] + YY);
                        Y[i][j] += phmm->emit_background[sb[j]];
                }

        }
        //fprintf(stdout,"%f %f %f\n", M[len_a][len_b]+phmm->theta, X[len_a][len_b] + phmm->theta, Y[len_a][len_b] + phmm->theta);
        phmm->f_score =                      M[len_a][len_b] + phmm->tau;
        phmm->f_score = logsum(phmm->f_score,X[len_a][len_b] + phmm->tau);
        phmm->f_score = logsum(phmm->f_score,Y[len_a][len_b] + phmm->tau);

//        LOG_MSG("Forward: %f\n",phmm->f_score);
        /* zero length alignment */
        //phmm->f_score = logsum(phmm->f_score, phmm->tau);
        //LOG_MSG("Forward: %f\n",phmm->f_score);
        return OK;
}


int backward_phmm(struct phmm* phmm,int* seq_a,int* seq_b, int len_a,int len_b)
{

        int i,j;
        int* sa = seq_a -1;
        int* sb = seq_b -1;
        const float MM = phmm->transition[INDEXMM];
        const float MX = phmm->transition[INDEXGPO];
        const float MY = phmm->transition[INDEXGPO];
        const float XX = phmm->transition[INDEXGPE];
        const float XM = phmm->transition[INDEXTM];
        const float YY = phmm->transition[INDEXGPE];
        const float YM = phmm->transition[INDEXTM];


        float** M = phmm->bM;
        float** X = phmm->bX;
        float** Y = phmm->bY;
/*const tAA = prob2scaledprob(1.0f - scaledprob2prob(phmm->eta));
        const tBB = prob2scaledprob(1.0f - scaledprob2prob(phmm->eta));
        const tAM = phmm->eta;
        const
tBM = phmm->eta;*/

        /* init first element  */
        M[len_a][len_b] = phmm->tau;
        X[len_a][len_b] = phmm->tau;
        Y[len_a][len_b] = phmm->tau;


        for(j = len_b-1; j > 0;j--){
                M[len_a][j] = Y[len_a][j+1] + MY +  phmm->emit_background[sb[j+1]];

                X[len_a][j] = prob2scaledprob(0.0f);

                //Y[len_a][j] =                     M[len_a][j+1] + YM + phmm->emit_M[sa[len_a]][sb[j+1]];
                Y[len_a][j] = Y[len_a][j+1] + YY + phmm->emit_background[sb[j+1]];

        }
        M[len_a][0] = prob2scaledprob(0.0f);
        X[len_a][0] = prob2scaledprob(0.0f);
        Y[len_a][0] = prob2scaledprob(0.0f);
        for(i = len_a-1; i >= 1;i--){

                M[i][len_b] = X[i+1][len_b] + MX + phmm->emit_background[sa[i+1]];


                X[i][len_b] = X[i+1][len_b] + XX + phmm->emit_background[sa[i+1]];


                Y[i][len_b] = prob2scaledprob(0.0f);

                for(j = len_b-1; j >=0;j--){

                        M[i][j] = M[i+1][j+1] + MM +phmm->emit_M[sa[i+1]][sb[j+1]];

                        M[i][j] = logsum(M[i][j], X[i+1][j] + MX + phmm->emit_background[sa[i+1]]);
                        M[i][j] = logsum(M[i][j], Y[i][j+1] + MY + phmm->emit_background[sb[j+1]]);


                        X[i][j] =                 M[i+1][j+1] + XM +phmm->emit_M[sa[i+1]][sb[j+1]];
                        X[i][j] = logsum(X[i][j], X[i+1][j] + XX + phmm->emit_background[sa[i+1]]);


                        Y[i][j] =                 M[i+1][j+1] + YM+ phmm->emit_M[sa[i+1]][sb[j+1]];
                        Y[i][j] = logsum(Y[i][j], Y[i][j+1] + YY+ phmm->emit_background[sb[j+1]]);
                }

        }
        i = 0;
        M[i][len_b] = prob2scaledprob(0.0f);
        X[i][len_b] = prob2scaledprob(0.0f);
        Y[i][len_b] = prob2scaledprob(0.0f);
        for(j = len_b-1; j >=1;j--){

                M[i][j] = prob2scaledprob(0.0f);



                X[i][j] = prob2scaledprob(0.0f);

                //X[i][j] =                 M[i+1][j+1] + XM +phmm->emit_M[sa[i+1]][sb[j+1]];
                //      X[i][j] = logsum(X[i][j], X[i+1][j] + XX + phmm->emit_background[sa[i+1]]);


                Y[i][j] =                 M[i+1][j+1] + YM+ phmm->emit_M[sa[i+1]][sb[j+1]];
                Y[i][j] = logsum(Y[i][j], Y[i][j+1] + YY+ phmm->emit_background[sb[j+1]]);
        }

        i = 0;
        j = 0;
        M[i][j] = M[i+1][j+1] + MM +phmm->emit_M[sa[i+1]][sb[j+1]];
        M[i][j] = logsum(M[i][j], X[i+1][j] + MX + phmm->emit_background[sa[i+1]]);
        M[i][j] = logsum(M[i][j], Y[i][j+1] + MY + phmm->emit_background[sb[j+1]]);

        phmm->b_score = M[0][0];

        /* zero length alignment */
//phmm->b_score = logsum(phmm->b_score, phmm->tau);
        return OK;
}


int collect_phmm(struct phmm* phmm,int* seq_a,int* seq_b, int len_a,int len_b)
{

        int i,j;
        int* sa = seq_a -1;
        int* sb = seq_b -1;


        const float MM = phmm->transition[INDEXMM];
        const float MX = phmm->transition[INDEXGPO];
        const float MY = phmm->transition[INDEXGPO];
        const float XX = phmm->transition[INDEXGPE];
        const float XM = phmm->transition[INDEXTM];
        const float YY = phmm->transition[INDEXGPE];
        const float YM = phmm->transition[INDEXTM];




        float eMM = phmm->transition_e[INDEXMM];
        float eGPO = phmm->transition_e[INDEXGPO];
        float eGPE = phmm->transition_e[INDEXGPE];
        float eTM = phmm->transition_e[INDEXTM];

        float** fM = phmm->fM;
        float** fX = phmm->fX;
        float** fY = phmm->fY;

        float** bM = phmm->bM;
        float** bX = phmm->bX;
        float** bY = phmm->bY;


        float score = phmm->f_score;

        /* init first element  */
        //M[0][0] = prob2scaledprob(1.0f);
        //X[0][0] = prob2scaledprob(0.0f);
        //Y[0][0] = prob2scaledprob(0.0f);

        for(j = 1; j<= len_b;j++){

                //Y[0][j] =                 M[0][j-1] + MY;
                eGPO = logsum(eGPO, fM[0][j-1] + MY + bY[0][j] + phmm->emit_background[sb[j]] - score);
                //Y[0][j] = logsum(Y[0][j], Y[0][j-1] + YY);
                eGPE = logsum(eGPE, fY[0][j-1] + YY + bY[0][j] + phmm->emit_background[sb[j]]- score);
                //Y[0][j] += phmm->emit_background[sb[j]];
                phmm->emit_background_e[sb[j]] = logsum(phmm->emit_background_e[sb[j]], fY[0][j] + bY[0][j]-score);
        }

        for(i = 1; i <= len_a;i++){
                //X[i][0] =                 M[i-1][0] + MX;
                eGPO = logsum(eGPO, fM[i-1][0] + MX + bX[i][0] + phmm->emit_background[sa[i]] - score);
                //X[i][0] = logsum(X[i][0], X[i-1][0] + XX);
                eGPE = logsum(eGPE, fX[i-1][0] + XX + bX[i][0] + phmm->emit_background[sa[i]] - score);



                //X[i][0] += phmm->emit_background[sa[i]];
                phmm->emit_background_e[sa[i]] = logsum(phmm->emit_background_e[sa[i]], fX[i][0]+ bX[i][0]-score);


                for(j = 1; j<= len_b;j++){
                        //M[i][j] = M[i-1][j-1] + MM;
                        eMM = logsum(eMM,fM[i-1][j-1] + MM + bM[i][j] + phmm->emit_M[sa[i]][sb[j]] - score);


                        //M[i][j] = logsum(M[i][j], X[i-1][j-1] + XM);
                        eTM = logsum(eTM,fX[i-1][j-1] + XM + bM[i][j] + phmm->emit_M[sa[i]][sb[j]] - score);
                        //M[i][j] = logsum(M[i][j], Y[i-1][j-1] + YM);
                        eTM = logsum(eTM,fY[i-1][j-1] + YM + bM[i][j] + phmm->emit_M[sa[i]][sb[j]] - score);
                        //M[i][j] += phmm->emit_M[sa[i]][sb[j]];
                        phmm->emit_M_e[sa[i]][sb[j]] = logsum(phmm->emit_M_e[sa[i]][sb[j]], fM[i][j] + bM[i][i] - score);


                        //X[i][j] =                 M[i-1][j] + MX;
                        eGPO = logsum(eGPO, fM[i-1][j] + MX + bX[i][j] + phmm->emit_background[sa[i]] - score);
                        //X[i][j] = logsum(X[i][j], X[i-1][j] + XX);
                        eGPE = logsum(eGPE, fX[i-1][j] + XX + bX[i][j] + phmm->emit_background[sa[i]] - score);

                        //X[i][j] += phmm->emit_background[sa[i]];
                        phmm->emit_background_e[sa[i]] = logsum(phmm->emit_background_e[sa[i]], fX[i][j]+ bX[i][j]-score);

                        //Y[i][j] =                 M[i][j-1] + MY;
                        eGPO = logsum(eGPO, fM[i][j-1] + MY + bY[i][j] + phmm->emit_background[sb[j]] - score);
                        //Y[i][j] = logsum(Y[i][j], Y[i][j-1] + YY);
                        eGPE = logsum(eGPE, fY[i][j-1] + YY + bY[i][j] + phmm->emit_background[sb[j]] - score);
                        //Y[i][j] += phmm->emit_background[sb[j]];
                        phmm->emit_background_e[sb[j]] = logsum(phmm->emit_background_e[sb[j]], fY[i][j] + bY[i][j]-score);
                }

        }
        //fprintf(stdout,"%f %f %f\n", M[len_a][len_b]+phmm->theta, X[len_a][len_b] + phmm->theta, Y[len_a][len_b] + phmm->theta);

        //phmm->tau_e = logsum(phmm->tau_e, fM[len_a][len_b] + phmm->tau - score);
        //phmm->tau_e = logsum(phmm->tau_e, fX[len_a][len_b] + phmm->tau - score);
        //phmm->tau_e = logsum(phmm->tau_e, fY[len_a][len_b] + phmm->tau - score);
//        phmm->tau_e = logsum(phmm->tau_e, phmm->tau - score);
        phmm->transition_e[INDEXMM] = eMM;
        phmm->transition_e[INDEXGPO] = eGPO;
        phmm->transition_e[INDEXGPE] = eGPE;
        phmm->transition_e[INDEXTM] = eTM;

        return OK;
}

int re_estimate(struct phmm* phmm)
{
        int i,j;
        int L;
        float sum;
        /* first for match state  */
        /*  */
        L = phmm->L;
        /*fprintf(stdout,"Estimated:\n");
        fprintf(stdout,"%f\tMM\n", scaledprob2prob(phmm->transition_e[INDEXMM]));
        fprintf(stdout,"%f\tGPO\n",scaledprob2prob(phmm->transition_e[INDEXGPO]));
        fprintf(stdout,"%f\tGPE\n",scaledprob2prob(phmm->transition_e[INDEXGPE]));
        fprintf(stdout,"%f\tTM\n",scaledprob2prob(phmm->transition_e[INDEXTM]));
        */


        sum = phmm->transition_e[INDEXMM];
        sum = logsum(sum, phmm->transition_e[INDEXGPO]);
        sum = logsum(sum, phmm->transition_e[INDEXGPO]);
        //sum = logsum(sum, phmm->tau_e);
        phmm->transition[INDEXMM] = phmm->transition_e[INDEXMM] - sum;
        phmm->transition[INDEXGPO] = phmm->transition_e[INDEXGPO] - sum;
        /* add in tau */
        sum = phmm->transition[INDEXMM];
        sum = logsum(sum, phmm->transition[INDEXGPO]);
        sum = logsum(sum, phmm->transition[INDEXGPO]);
        sum = logsum(sum, phmm->tau);

        phmm->transition[INDEXMM] = phmm->transition[INDEXMM] - sum;
        phmm->transition[INDEXGPO] = phmm->transition[INDEXGPO] - sum;


        /* sanity check */
        sum = phmm->transition[INDEXMM];
        sum = logsum(sum, phmm->transition[INDEXGPO]);
        sum = logsum(sum, phmm->transition[INDEXGPO]);
        sum = logsum(sum, phmm->tau);
        //fprintf(stdout,"sanity: M: %f\n", scaledprob2prob(sum));

        /* next delete state */
        sum = phmm->transition_e[INDEXGPE];
        sum = logsum(sum, phmm->transition_e[INDEXTM]);
        //sum = logsum(sum, phmm->tau);
        phmm->transition[INDEXGPE] =  phmm->transition_e[INDEXGPE] - sum;
        phmm->transition[INDEXTM] =  phmm->transition_e[INDEXTM] - sum;
        /* add in tau */
        sum = phmm->transition[INDEXGPE];
        sum = logsum(sum, phmm->transition[INDEXTM]);
        sum = logsum(sum, phmm->tau);

        phmm->transition[INDEXGPE] =  phmm->transition[INDEXGPE] - sum;
        phmm->transition[INDEXTM] =  phmm->transition[INDEXTM] - sum;
        

        sum = phmm->transition[INDEXGPE];
        sum = logsum(sum, phmm->transition[INDEXTM]);
        sum = logsum(sum, phmm->tau);

        //fprintf(stdout,"sanity: X/Y: %f\n", scaledprob2prob(sum));

        sum = prob2scaledprob(0.0f);
        for(i = 0; i < L;i++){

                sum = logsum(sum, phmm->emit_background_e[i]);

        }

        for(i = 0; i < L;i++){
                phmm->emit_background[i] = phmm->emit_background_e[i] -sum;
                fprintf(stdout,"%d\t%f\n",i,scaledprob2prob(phmm->emit_background[i]));
        }
        sum = prob2scaledprob(0.0f);
        for(i = 0; i < L;i++){
                sum = logsum(sum, phmm->emit_background[i]);
        }
        //fprintf(stdout,"sanity:background: %f\n", scaledprob2prob(sum));

        sum = prob2scaledprob(0.0f);
        for(i = 0; i < L;i++){
                for(j = 0; j < L;j++){
                        sum = logsum(sum, phmm->emit_M_e[i][j]);
                }
        }

        for(i = 0; i < L;i++){
                for(j = 0; j < L;j++){
                        phmm->emit_M[i][j] = phmm->emit_M_e[i][j] - sum;
                        fprintf(stdout,"%0.2f ",scaledprob2prob(phmm->emit_M[i][j]));
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        sum = prob2scaledprob(0.0f);
        for(i = 0; i < L;i++){
                for(j = 0; j < L;j++){
                        sum = logsum(sum, phmm->emit_M[i][j]);
                }
        }
        //fprintf(stdout,"sanity:emit: %f\n", scaledprob2prob(sum));
        clear_phmm_e(phmm);
        return OK;
}

int phmm_transitions(struct phmm* phmm)
{
        fprintf(stdout,"%f\tMM\n", scaledprob2prob(phmm->transition[INDEXMM]));
        fprintf(stdout,"%f\tGPO\n",scaledprob2prob(phmm->transition[INDEXGPO]));
        fprintf(stdout,"%f\tGPE\n",scaledprob2prob(phmm->transition[INDEXGPE]));
        fprintf(stdout,"%f\tTM\n",scaledprob2prob(phmm->transition[INDEXTM]));
        return OK;
}

int print_phmm(struct phmm* phmm,int len_a,int len_b)
{
        int i,j;
        float** M = phmm->fM;
        float** X = phmm->fX;
        float** Y = phmm->fY;
        for(i = 0; i <= len_a;i++){
                for(j = 0; j<= len_b;j++){
                        fprintf(stdout,"%2.2f,%2.2f,%2.2f ",M[i][j],X[i][j],Y[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        fprintf(stdout,"\n");

        M = phmm->bM;
        X = phmm->bX;
        Y = phmm->bY;
        for(i = 0; i <= len_a;i++){
                for(j = 0; j<= len_b;j++){
                        fprintf(stdout,"%2.2f,%2.2f,%2.2f ",M[i][j],X[i][j],Y[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        return OK;
}

int simple_init(struct phmm*phmm)
{
        float mismatch = 0.01f;
        float alphabet = 20;
        float sum = 0.0;


        int i,j;

        phmm->L = alphabet;
        for(i =0; i < alphabet;i++){
                phmm->emit_background[i] = prob2scaledprob(1.0 / (float) alphabet);
                phmm->emit_background_e[i] = prob2scaledprob(0.0f);
                sum = prob2scaledprob(0.0f);
                for (j=0; j< alphabet; j++){
                        if(i == j){
                                phmm->emit_M[i][j] = prob2scaledprob(1.0 - (alphabet-1)* mismatch);
                        }else{
                                phmm->emit_M[i][j] = prob2scaledprob(mismatch);
                        }
                        phmm->emit_M_e[i][j] = prob2scaledprob(0.0f);
                        sum = logsum(sum,phmm->emit_M[i][j]);
                }
//                fprintf(stdout,"sanity check %d: %f\n",i, scaledprob2prob(sum));
        }


        int m_pos = 0;
        short *matrix_pointer = 0;
        short blosum50mt[]={
                5,
                -2,  5,
                -1, -3, 13,
                -2,  5, -4,  8,
                -1,  1, -3,  2,  6,
                -3, -4, -2, -5, -3,  8,
                0, -1, -3, -1, -3, -4,  8,
                -2,  0, -3, -1,  0, -1, -2, 10,
                -1, -4, -2, -4, -4,  0, -4, -4,  5,
                -1,  0, -3, -1,  1, -4, -2,  0, -3,  6,
                -2, -4, -2, -4, -3,  1, -4, -3,  2, -3,  5,
                -1, -3, -2, -4, -2,  0, -3, -1,  2, -2,  3,  7,
                -1,  4, -2,  2,  0, -4,  0,  1, -3,  0, -4, -2,  7,
                -1, -2, -4, -1, -1, -4, -2, -2, -3, -1, -4, -3, -2, 10,
                -1,  0, -3,  0,  2, -4, -2,  1, -3,  2, -2,  0,  0, -1,  7,
                -2, -1, -4, -2,  0, -3, -3,  0, -4,  3, -3, -2, -1, -3,  1,  7,
                1,  0, -1,  0, -1, -3,  0, -1, -3,  0, -3, -2,  1, -1,  0, -1,  5,
                0,  0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  2,  5,
                0, -4, -1, -4, -3, -1, -4, -4,  4, -3,  1,  1, -3, -3, -3, -3, -2,  0,  5,
                -3, -5, -5, -5, -3,  1, -3, -3, -3, -3, -2, -1, -4, -4, -1, -3, -4, -3, -3, 15,
                -1, -1, -2, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1, -2, -1, -1, -1,  0, -1, -3, -1,
                -2, -3, -3, -3, -2,  4, -3,  2, -1, -2, -1,  0, -2, -3, -1, -1, -2, -2, -1,  2, -1,  8,
                -1,  2, -3,  1,  5, -4, -2,  0, -3,  1, -3, -1,  0, -1,  4,  0,  0, -1, -3, -2, -1, -2,  5};
                                matrix_pointer = blosum50mt;

        for (i = 0;i < 23;i++){
                for (j = 0;j <= i;j++){
                        if (i == j){
                                //	subm[i][j] += blosum62mt[m_pos]*10;
                                phmm->emit_M[i][j] = exp((float)matrix_pointer[m_pos]) * 1.0 / (float) alphabet *  1.0 / (float) alphabet ;
                        }else{
                                //	subm[i][j] += blosum62mt[m_pos]*10;
                                //	subm[j][i] += blosum62mt[m_pos]*10;
                                phmm->emit_M[i][j] = exp((float)matrix_pointer[m_pos]) * 1.0 / (float) alphabet *  1.0 / (float) alphabet ;
                                phmm->emit_M[j][i] = exp((float)matrix_pointer[m_pos]) * 1.0 / (float) alphabet *  1.0 / (float) alphabet ;
                        }
                        m_pos++;
                }
        }
        sum = 0.0;
        for(i = 0 ; i < alphabet;i++){
                for(j = 0; j < alphabet;j++){
                        sum += phmm->emit_M[i][j];// = prob2scaledprob(phmm->emit_M[i][j] / sum);

                }

        }

        for(i = 0 ; i < alphabet;i++){
                for(j = 0; j < alphabet;j++){
                        phmm->emit_M[i][j] = prob2scaledprob(phmm->emit_M[i][j] / sum);

                }

        }

        sum = prob2scaledprob(0.0);
        for(i = 0 ; i < alphabet;i++){
                for(j = 0; j < alphabet;j++){

                        sum = logsum(sum, phmm->emit_M[i][j]);
                }

        }
        //fprintf(stdout,"emit sanity %f\n", scaledprob2prob(sum));

        phmm->tau = prob2scaledprob(0.01f);
        phmm->transition[INDEXMM] = prob2scaledprob(0.89f);
        phmm->transition[INDEXGPO] = prob2scaledprob(0.05f);



        sum = phmm->transition[INDEXMM];
        sum = logsum(sum, phmm->transition[INDEXGPO]);
        sum = logsum(sum, phmm->transition[INDEXGPO]);
        sum = logsum(sum, phmm->tau);

        //fprintf(stdout,"sanity: M: %f\n", scaledprob2prob(sum));


        phmm->transition[INDEXGPE] = prob2scaledprob(0.09f);
        phmm->transition[INDEXTM] = prob2scaledprob(0.90f);

        sum = phmm->transition[INDEXGPE];
        sum = logsum(sum, phmm->transition[INDEXTM]);
        sum = logsum(sum, phmm->tau);

        //fprintf(stdout,"sanity: X/Y: %f\n", scaledprob2prob(sum));
        //exit(0);
        return OK;
}


int clear_phmm_e(struct phmm* phmm)
{
        int i,j;

        phmm->transition_e[INDEXMM] = prob2scaledprob(1.0f);
        phmm->transition_e[INDEXGPO] = prob2scaledprob(1.0f);
        phmm->transition_e[INDEXGPE] = prob2scaledprob(1.0f);
        phmm->transition_e[INDEXTM] = prob2scaledprob(1.0f);

        for(i = 0;i < phmm->L;i++){
                phmm->emit_background_e[i] = prob2scaledprob(1.0f);
                for(j = 0; j < phmm->L;j++){
                        if(i == j){
                                phmm->emit_M_e[i][j] = prob2scaledprob(100.0f);
                        }else{
                                phmm->emit_M_e[i][j] = prob2scaledprob(1.0f);
                        }
                }
        }

        return OK;
}




struct phmm* alloc_phmm(int size)
{
        struct phmm* phmm = NULL;

        init_logsum();

        MMALLOC(phmm, sizeof(struct phmm));



        phmm->fM = NULL;
        phmm->fX = NULL;
        phmm->fY = NULL;

        phmm->bM = NULL;
        phmm->bX = NULL;
        phmm->bY = NULL;

        phmm->fM = galloc(phmm->fM,size+2, size+2,0.0f);
        phmm->fX = galloc(phmm->fX,size+2, size+2,0.0f);
        phmm->fY = galloc(phmm->fY,size+2, size+2,0.0f);

        phmm->bM = galloc(phmm->bM,size+2, size+2,0.0f);
        phmm->bX = galloc(phmm->bX,size+2, size+2,0.0f);
        phmm->bY = galloc(phmm->bY,size+2, size+2,0.0f);


        return phmm;
ERROR:
        free_phmm(phmm);
        return NULL;
}

void free_phmm(struct phmm* phmm)
{
        if(phmm){

                gfree(phmm->fM);
                gfree(phmm->fX);
                gfree(phmm->fY);

                gfree(phmm->bM);
                gfree(phmm->bX);
                gfree(phmm->bY);

                MFREE(phmm);
        }
}

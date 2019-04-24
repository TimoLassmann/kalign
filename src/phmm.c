

#include "phmm.h"





/* need forward backward  */
/* or be adventurous and use beam sampling with u =1 i.e. no threshold ... worth a try?  */

struct phmm{

        float** fM ;
        float** fX ;
        float** fY ;

        float** bM ;
        float** bX ;
        float** bY ;


        int alloc_x;
        int alloc_y;

        float emit_M[26][26];
        float emit_background[26];

        float emit_M_e[26][26];
        float emit_background_e[26];

        float delta;           /* gap open */
        float epsilon;         /* gap extension */
        float theta;          /* exit model */
        float eta;             /* exit random  */

        float delta_e;           /* gap open */
        float epsilon_e;         /* gap extension */
        float theta_e;          /* exit model */
        float eta_e;             /* exit random  */

        float f_score;
        float b_score;
};

int forward_phmm(struct phmm* phmm,int* seq_a,int* seq_b, int len_a,int len_b);
int backward_phmm(struct phmm* phmm,int* seq_a,int* seq_b, int len_a,int len_b);


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

        RUN(forward_phmm(phmm, seqa, seqb, len_a, len_b));
        RUN(backward_phmm(phmm, seqa, seqb, len_a, len_b));
        fprintf(stdout,"%f\tforward\n%f\tbackward\n",phmm->f_score, phmm->b_score);
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
        const float MM = prob2scaledprob(1.0f - 2.0f * scaledprob2prob(phmm->delta) - scaledprob2prob(phmm->theta));
        const float MX = phmm->delta;
        const float MY = phmm->delta;
        const float XM = prob2scaledprob(1.0 - scaledprob2prob(phmm->epsilon) - scaledprob2prob(phmm->theta));
        const float YM = XM;
        const float XX = phmm->epsilon;
        const float YY = phmm->epsilon;


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

        for(j = 1; j<= len_b;j++){
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

        phmm->f_score =                      M[len_a][len_b] + phmm->theta;
        phmm->f_score = logsum(phmm->f_score,X[len_a][len_b] + phmm->theta);
        phmm->f_score = logsum(phmm->f_score,Y[len_a][len_b] + phmm->theta);

        return OK;
}


int backward_phmm(struct phmm* phmm,int* seq_a,int* seq_b, int len_a,int len_b)
{

        int i,j;
        int* sa = seq_a -1;
        int* sb = seq_b -1;
        const float MM = prob2scaledprob(1.0f - 2.0f * scaledprob2prob(phmm->delta) - scaledprob2prob(phmm->theta));
        const float MX = phmm->delta;
        const float MY = phmm->delta;
        const float XM = prob2scaledprob(1.0 - scaledprob2prob(phmm->epsilon) - scaledprob2prob(phmm->theta));
        const float YM = XM;
        const float XX = phmm->epsilon;
        const float YY = phmm->epsilon;


        float** M = phmm->bM;
        float** X = phmm->bX;
        float** Y = phmm->bY;
/*const tAA = prob2scaledprob(1.0f - scaledprob2prob(phmm->eta));
        const tBB = prob2scaledprob(1.0f - scaledprob2prob(phmm->eta));
        const tAM = phmm->eta;
        const tBM = phmm->eta;*/

        /* init first element  */
        M[len_a][len_b] = phmm->theta + phmm->emit_M[sa[len_a]][sb[len_b]];
        X[len_a][len_b] = phmm->theta + phmm->emit_background[sa[len_a]];
        Y[len_a][len_b] = phmm->theta + phmm->emit_background[sb[len_b]];


        for(j = len_b-1; j >= 1;j--){
                M[len_a][j] = prob2scaledprob(0.0f);

                X[len_a][j] = prob2scaledprob(0.0f);

                Y[len_a][j] =                     M[len_a][j+1] + YM;
                Y[len_a][j] = logsum(Y[len_a][j], Y[len_a][j+1] + YY);
                Y[len_a][j] += phmm->emit_background[sb[j]];
        }
        /* top right corner */
        M[len_a][0] = prob2scaledprob(0.0f);
        X[len_a][0] = prob2scaledprob(0.0f);
        Y[len_a][0] = prob2scaledprob(0.0f);


        for(i = len_a-1; i >= 1;i--){

                M[i][len_b] = X[i+1][len_b] + MX;
                M[i][len_b] += phmm->emit_M[sa[i]][sb[len_b]];

                X[i][len_b] = X[i+1][len_b] + XX;
                X[i][len_b] += phmm->emit_background[sa[i]];

                Y[i][len_a] =prob2scaledprob(0.0f);

                for(j = len_b-1; j >=1;j--){

                        M[i][j] = M[i+1][j+1] + MM;
                        M[i][j] = logsum(M[i][j], X[i+1][j] + MX);
                        M[i][j] = logsum(M[i][j], Y[i][j+1] + MY);
                        M[i][j] += phmm->emit_M[sa[i]][sb[j]];


                        X[i][j] =                 M[i+1][j+1] + XM;
                        X[i][j] = logsum(X[i][j], X[i+1][j] + XX);
                        X[i][j] += phmm->emit_background[sa[i]];

                        Y[i][j] =                 M[i+1][j+1] + YM;
                        Y[i][j] = logsum(Y[i][j], Y[i][j+1] + YY);
                        Y[i][j] += phmm->emit_background[seq_b[j]];
                }
                M[i][0] = prob2scaledprob(0.0f);

                X[i][0] =               M[i+1][j+1] + XM;
                X[i][0] = logsum(X[i][j], X[i+1][j] + XX);
                X[i][0] += phmm->emit_background[sa[i]];

                Y[i][0] = prob2scaledprob(0.0f); /* dead end */
        }

        /* first column */
        /* bottom right corner */

        M[0][len_b] = prob2scaledprob(0.0f);
        X[0][len_b] = prob2scaledprob(0.0f);
        Y[0][len_b] = prob2scaledprob(0.0f);

        for(j = len_b-1; j >=1;j--){


                M[0][j] = prob2scaledprob(0.0f);

                X[0][j] = prob2scaledprob(0.0f); /* dead end */

                Y[0][j] =                 M[1][j+1] + YM;
                Y[0][j] = logsum(Y[0][j], Y[0][j+1] + YY);
                Y[0][j] += phmm->emit_background[seq_b[j]];
        }

        M[0][0] = M[1][1] + MM;
        M[0][0] = logsum(M[0][0], X[1][0] + MX);
        M[0][0] = logsum(M[0][0], Y[0][1] + MY);
        phmm->b_score = M[0][0];

        return OK;
}



int simple_init(struct phmm*phmm)
{
        float mismatch = 0.05f;
        float alphabet = 4;
        float sum = 0.0;


        int i,j;

        for(i =0; i < alphabet;i++){
                phmm->emit_background[i] = prob2scaledprob(0.25f);
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
                fprintf(stdout,"sanity check %d: %f\n",i, scaledprob2prob(sum));

        }
        phmm->delta = prob2scaledprob(0.05f);
        phmm->epsilon = prob2scaledprob(0.5f);
        phmm->theta = prob2scaledprob(0.1f);
        phmm->eta = prob2scaledprob(0.1f);

        phmm->delta_e = prob2scaledprob(0.0f);
        phmm->epsilon_e = prob2scaledprob(0.0f);
        phmm->theta_e = prob2scaledprob(0.0f);
        phmm->eta_e = prob2scaledprob(0.0f);
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

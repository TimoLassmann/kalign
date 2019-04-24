

#include "phmm.h"





/* need forward backward  */
/* or be adventurous and use beam sampling with u =1 i.e. no threshold ... worth a try?  */

struct point{
        float M;
        float X;
        float Y;
        float null;

};

struct phmm{
        struct point* m;
        int alloc_m;

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
};

int forward_phmm(struct phmm* phmm,int* seq_a,int* seq_b, int len_a,int len_b);

int simple_init(struct phmm*phmm);
/* memory functions */
struct phmm* alloc_phmm(int x, int y);
int realloc_phmm(struct phmm* phmm, int x, int y);
void free_phmm(struct phmm* phmm);

int main(int argc, char *argv[])
{
        struct phmm* phmm = NULL;
        int* seqa = NULL;
        int* seqb = NULL;

        int len_a = 5;
        int len_b = 5;

        RUNP(phmm = alloc_phmm(len_a,len_b));
        RUN(simple_init(phmm));
        free_phmm(phmm);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


int forward_phmm(struct phmm* phmm,int* seq_a,int* seq_b, int len_a,int len_b)
{
        struct point* prev;
        struct point* cur;
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

        /*const tAA = prob2scaledprob(1.0f - scaledprob2prob(phmm->eta));
        const tBB = prob2scaledprob(1.0f - scaledprob2prob(phmm->eta));
        const tAM = phmm->eta;
        const tBM = phmm->eta;*/

        cur = phmm->m;


        prev = phmm->m;



        /* init first element  */
        cur[0].M  = prob2scaledprob(1.0f);
        cur[0].X = prob2scaledprob(0.0f);
        cur[0].Y = prob2scaledprob(0.0f);

        for(j = 1;j <= len_b;j++){
                cur[j].M = prob2scaledprob(0.0f);
                cur[j].X =                  cur[j-1].M + MX;
                cur[j].X = logsum(cur[j].X, cur[j-1].X + XX);
                cur[j].X += phmm->emit_background[sb[j]];
                cur[j].Y = prob2scaledprob(0.0f);
                //phmm->m[j].M

        }

        cur += len_b+1;
        /* init top row
 */

        for(i = 1; i <= len_a;i++){
                cur[0].M = prob2scaledprob(0.0f);
                cur[0].X = prob2scaledprob(0.0f);

                cur[0].Y = prev[0].M + MY;
                cur[0].Y = logsum(cur[0].Y,  prev[0].Y + YY);
                cur[0].Y += phmm->emit_background[sa[i]];
                for(j = 1;j <= len_b;j++){

                        cur[j].M =                 prev[j-1].M + MM;
                        cur[j].M = logsum(cur[j].M,prev[j-1].X + XM);
                        cur[j].M = logsum(cur[j].M,prev[j-1].Y + YM);
                        cur[j].M += phmm->emit_M[sa[i]][sb[j]];

                        cur[j].X =                  cur[j-1].M + MX;
                        cur[j].X = logsum(cur[j].X, cur[j-1].X + XX);
                        cur[j].X += phmm->emit_background[sb[j]];

                        cur[j].Y =                   prev[j].M + MY;
                        cur[j].Y = logsum(cur[j].Y,  prev[j].Y + YY);
                        cur[j].Y += phmm->emit_background[sa[i]];


                        //phmm->m[j].M

                }
                prev =cur;
                cur += len_b+1;
        }



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







struct phmm* alloc_phmm(int x, int y)
{
        struct phmm* phmm = NULL;

        init_logsum();

        MMALLOC(phmm, sizeof(struct phmm));

        phmm->m = NULL;
        phmm->alloc_m = (x+2) * (y+2);

        MMALLOC(phmm->m, sizeof(struct point) * phmm->alloc_m);


        return phmm;
ERROR:
        free_phmm(phmm);
        return NULL;
}

int realloc_phmm(struct phmm* phmm, int x, int y)
{
        int size = (x+2)*(y+2);
        if(size > phmm->alloc_m){
                phmm->alloc_m = size;
                MREALLOC(phmm->m, sizeof(struct point) * phmm->alloc_m);
        }
        return OK;
ERROR:

        return FAIL;
}

void free_phmm(struct phmm* phmm)
{
        if(phmm){
                if(phmm->m){
                        MFREE(phmm->m);
                }
                MFREE(phmm);
        }
}



#include "phmm.h"


/* need forward backward  */
/* or be adventurous and use beam sampling with u =1 i.e. no threshold ... worth a try?  */


#ifdef ITEST
int main(int argc, char *argv[])
{
        struct phmm* phmm = NULL;
        uint8_t* seqa = NULL;
        uint8_t* seqb = NULL;

        int len_a = 5;
        int len_b = 5;

        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};


        /* example from biological sequence analysis */

        len_a = 7;
        len_b = 10;
        MMALLOC(seqa, sizeof(uint8_t) * len_a);
        MMALLOC(seqb, sizeof(uint8_t) * len_b);

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
        //RUN(add_pseudocounts(phmm, 20));

        //RUN(phmm_transitions(phmm));
        for(int iter = 0;iter < 50;iter++){
                RUN(forward_phmm(phmm, seqa, seqb, len_a, len_b));
                RUN(backward_phmm(phmm, seqa, seqb, len_a, len_b));


                RUN(collect_phmm(phmm, seqa, seqb, len_a, len_b));
                RUN(re_estimate(phmm));
                RUN(add_pseudocounts(phmm, 20));
                fprintf(stdout,"%f\tforward\n%f\tbackward\n",phmm->f_score, phmm->b_score);
        }
        RUN(phmm_transitions(phmm));
        //print_phmm(phmm, len_a, len_b);
        free_phmm(phmm);

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

#endif
int forward_phmm(struct phmm* phmm, uint8_t* seq_a, uint8_t* seq_b, int len_a,int len_b)
{

        int i,j;
        uint8_t* sa = seq_a -1;
        uint8_t* sb = seq_b -1;
        const float MM = phmm->transition[INDEXMM];
        const float MX = phmm->transition[INDEXGPO];
        const float MY = phmm->transition[INDEXGPO];
        const float XX = phmm->transition[INDEXGPE];
        const float XM = phmm->transition[INDEXTM];
        const float YY = phmm->transition[INDEXGPE];
        const float YM = phmm->transition[INDEXTM];
        /* NEEDS to be changed later */
        const float TGPE = phmm->eta;// phmm->transition[INDEXGPE];



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

        for(j = 1; j < len_b;j++){
                M[0][j] = prob2scaledprob(0.0f);

                X[0][j] = prob2scaledprob(0.0f);

                Y[0][j] =                 M[0][j-1] + MY;
                Y[0][j] = logsum(Y[0][j], Y[0][j-1] + TGPE);
                Y[0][j] += phmm->emit_background[sb[j]];
        }
        M[0][len_b] = prob2scaledprob(0.0f);
        X[0][len_b] = prob2scaledprob(0.0f);
        Y[0][len_b] = prob2scaledprob(0.0f);

        for(i = 1; i <= len_a;i++){
                M[i][0] = prob2scaledprob(0.0f);

                X[i][0] =                 M[i-1][0] + MX;
                X[i][0] = logsum(X[i][0], X[i-1][0] + TGPE);
                X[i][0] += phmm->emit_background[sa[i]];

                Y[i][0] = prob2scaledprob(0.0f);

                for(j = 1; j< len_b;j++){
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

                j = len_b;
                M[i][j] = M[i-1][j-1] + MM;
                M[i][j] = logsum(M[i][j], X[i-1][j-1] + XM);
                M[i][j] = logsum(M[i][j], Y[i-1][j-1] + YM);
                M[i][j] += phmm->emit_M[sa[i]][sb[j]];

                X[i][j] =                 M[i-1][j] + MX;
                X[i][j] = logsum(X[i][j], X[i-1][j] + TGPE);
                X[i][j] += phmm->emit_background[sa[i]];

                Y[i][j] = prob2scaledprob(0.0f);//                M[i][j-1] + MY;
                //Y[i][j] = logsum(Y[i][j], Y[i][j-1] + YY);
                //Y[i][j] += phmm->emit_background[sb[j]];
        }
        i = len_a;
        M[i][0] = prob2scaledprob(0.0f);
        X[i][0] = prob2scaledprob(0.0f);
        Y[i][0] = prob2scaledprob(0.0f);

        for(j = 1; j < len_b;j++){
                M[i][j] = M[i-1][j-1] + MM;
                M[i][j] = logsum(M[i][j], X[i-1][j-1] + XM);
                M[i][j] = logsum(M[i][j], Y[i-1][j-1] + YM);
                M[i][j] += phmm->emit_M[sa[i]][sb[j]];

                X[i][j] = prob2scaledprob(0.0f);//                 M[i-1][j] + MX;
                //X[i][j] = logsum(X[i][j], X[i-1][j] + XX);
                //X[i][j] += phmm->emit_background[sa[i]];

                Y[i][j] =                 M[i][j-1] + MY;
                Y[i][j] = logsum(Y[i][j], Y[i][j-1] + TGPE);
                Y[i][j] += phmm->emit_background[sb[j]];
        }
        /* last cell */
        j = len_b;
        M[i][j] = M[i-1][j-1] + MM;
        M[i][j] = logsum(M[i][j], X[i-1][j-1] + XM);
        M[i][j] = logsum(M[i][j], Y[i-1][j-1] + YM);
        M[i][j] += phmm->emit_M[sa[i]][sb[j]];

        X[i][j] =                 M[i-1][j] + MX;
        X[i][j] = logsum(X[i][j], X[i-1][j] + TGPE);
        X[i][j] += phmm->emit_background[sa[i]];

        Y[i][j] =                 M[i][j-1] + MY;
        Y[i][j] = logsum(Y[i][j], Y[i][j-1] + TGPE);
        Y[i][j] += phmm->emit_background[sb[j]];


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


int backward_phmm(struct phmm* phmm,uint8_t* seq_a, uint8_t* seq_b, int len_a,int len_b)
{

        int i,j;
        uint8_t* sa = seq_a -1;
        uint8_t* sb = seq_b -1;
        const float MM = phmm->transition[INDEXMM];
        const float MX = phmm->transition[INDEXGPO];
        const float MY = phmm->transition[INDEXGPO];
        const float XX = phmm->transition[INDEXGPE];
        const float XM = phmm->transition[INDEXTM];
        const float YY = phmm->transition[INDEXGPE];
        const float YM = phmm->transition[INDEXTM];

        const float TGPE = phmm->eta;//
        //const float TGPE = phmm->transition[INDEXGPE];

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
                Y[len_a][j] = Y[len_a][j+1] + TGPE + phmm->emit_background[sb[j+1]];

        }
        M[len_a][0] = prob2scaledprob(0.0f);
        X[len_a][0] = prob2scaledprob(0.0f);
        Y[len_a][0] = prob2scaledprob(0.0f);

        for(i = len_a-1; i >= 1;i--){

                M[i][len_b] = X[i+1][len_b] + MX + phmm->emit_background[sa[i+1]];


                X[i][len_b] = X[i+1][len_b] + TGPE + phmm->emit_background[sa[i+1]];


                Y[i][len_b] = prob2scaledprob(0.0f);

                for(j = len_b-1; j >=1;j--){

                        M[i][j] = M[i+1][j+1] + MM +phmm->emit_M[sa[i+1]][sb[j+1]];

                        M[i][j] = logsum(M[i][j], X[i+1][j] + MX + phmm->emit_background[sa[i+1]]);
                        M[i][j] = logsum(M[i][j], Y[i][j+1] + MY + phmm->emit_background[sb[j+1]]);


                        X[i][j] =                 M[i+1][j+1] + XM +phmm->emit_M[sa[i+1]][sb[j+1]];
                        X[i][j] = logsum(X[i][j], X[i+1][j] + XX + phmm->emit_background[sa[i+1]]);


                        Y[i][j] =                 M[i+1][j+1] + YM+ phmm->emit_M[sa[i+1]][sb[j+1]];
                        Y[i][j] = logsum(Y[i][j], Y[i][j+1] + YY+ phmm->emit_background[sb[j+1]]);
                }
                j = 0;
                M[i][j] = prob2scaledprob(0.0f);//  M[i+1][j+1] + MM +phmm->emit_M[sa[i+1]][sb[j+1]];

                //M[i][j] = logsum(M[i][j], X[i+1][j] + MX + phmm->emit_background[sa[i+1]]);
                //M[i][j] = logsum(M[i][j], Y[i][j+1] + MY + phmm->emit_background[sb[j+1]]);


                X[i][j] =                 M[i+1][j+1] + XM +phmm->emit_M[sa[i+1]][sb[j+1]];
                X[i][j] = logsum(X[i][j], X[i+1][j] + TGPE + phmm->emit_background[sa[i+1]]);


                Y[i][j] = prob2scaledprob(0.0f);//                 M[i+1][j+1] + YM+ phmm->emit_M[sa[i+1]][sb[j+1]];
                //Y[i][j] = logsum(Y[i][j], Y[i][j+1] + YY+ phmm->emit_background[sb[j+1]]);

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


                Y[i][j] =                 M[i+1][j+1] + YM + phmm->emit_M[sa[i+1]][sb[j+1]];
                Y[i][j] = logsum(Y[i][j], Y[i][j+1] + TGPE + phmm->emit_background[sb[j+1]]);
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


int collect_phmm(struct phmm* phmm,uint8_t* seq_a,uint8_t* seq_b, int len_a,int len_b)
{

        int i,j;
        uint8_t* sa = seq_a -1;
        uint8_t* sb = seq_b -1;


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
//                eGPE = logsum(eGPE, fY[0][j-1] + YY + bY[0][j] + phmm->emit_background[sb[j]]- score);
                //Y[0][j] += phmm->emit_background[sb[j]];
                phmm->emit_background_e[sb[j]] = logsum(phmm->emit_background_e[sb[j]], fY[0][j] + bY[0][j]-score);
        }

        for(i = 1; i <= len_a;i++){
                //X[i][0] =                 M[i-1][0] + MX;
                eGPO = logsum(eGPO, fM[i-1][0] + MX + bX[i][0] + phmm->emit_background[sa[i]] - score);
                //X[i][0] = logsum(X[i][0], X[i-1][0] + XX);
                //eGPE = logsum(eGPE, fX[i-1][0] + XX + bX[i][0] + phmm->emit_background[sa[i]] - score);



                //X[i][0] += phmm->emit_background[sa[i]];
                phmm->emit_background_e[sa[i]] = logsum(phmm->emit_background_e[sa[i]], fX[i][0]+ bX[i][0]-score);


                for(j = 1; j < len_b;j++){
                        //M[i][j] = M[i-1][j-1] + MM;
                        eMM = logsum(eMM,fM[i-1][j-1] + MM + bM[i][j] + phmm->emit_M[sa[i]][sb[j]] - score);


                        //M[i][j] = logsum(M[i][j], X[i-1][j-1] + XM);
                        eTM = logsum(eTM,fX[i-1][j-1] + XM + bM[i][j] + phmm->emit_M[sa[i]][sb[j]] - score);
                        //M[i][j] = logsum(M[i][j], Y[i-1][j-1] + YM);
                        eTM = logsum(eTM,fY[i-1][j-1] + YM + bM[i][j] + phmm->emit_M[sa[i]][sb[j]] - score);
                        //M[i][j] += phmm->emit_M[sa[i]][sb[j]];
                        if(sa[i] == sb[j]){
                                //phmm->emit_M_e[sa[i]][sb[j]] = logsum(phmm->emit_M_e[sa[i]][sb[j]], fM[i][j] + bM[i][j] - score);
                                phmm->emit_M_e[sa[i]][sb[j]] = logsum(phmm->emit_M_e[sa[i]][sb[j]], fM[i][j] + bM[i][j] - score);
                        }else if(sa[i] > sb[j]){
                                phmm->emit_M_e[sa[i]][sb[j]] = logsum(phmm->emit_M_e[sa[i]][sb[j]], fM[i][j] + bM[i][j] - score);
                        }else{
                                phmm->emit_M_e[sb[j]][sa[i]] = logsum(phmm->emit_M_e[sb[j]][sa[i]], fM[i][j] + bM[i][j] - score);
                        }
                        //fprintf(stdout," %d -%d : %f\n", sa[i],sb[j],fM[i][j] + bM[i][j] - score);

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
                j = len_b;
                            //M[i][j] = M[i-1][j-1] + MM;
                eMM = logsum(eMM,fM[i-1][j-1] + MM + bM[i][j] + phmm->emit_M[sa[i]][sb[j]] - score);


                //M[i][j] = logsum(M[i][j], X[i-1][j-1] + XM);
                eTM = logsum(eTM,fX[i-1][j-1] + XM + bM[i][j] + phmm->emit_M[sa[i]][sb[j]] - score);
                //M[i][j] = logsum(M[i][j], Y[i-1][j-1] + YM);
                eTM = logsum(eTM,fY[i-1][j-1] + YM + bM[i][j] + phmm->emit_M[sa[i]][sb[j]] - score);
                //M[i][j] += phmm->emit_M[sa[i]][sb[j]];
                if(sa[i] == sb[j]){
                        //                    phmm->emit_M_e[sa[i]][sb[j]] = logsum(phmm->emit_M_e[sa[i]][sb[j]], fM[i][j] + bM[i][j] - score);
                        phmm->emit_M_e[sa[i]][sb[j]] = logsum(phmm->emit_M_e[sa[i]][sb[j]], fM[i][j] + bM[i][j] - score);
                }else if(sa[i] > sb[j]){
                        phmm->emit_M_e[sa[i]][sb[j]] = logsum(phmm->emit_M_e[sa[i]][sb[j]], fM[i][j] + bM[i][j] - score);
                }else{
                        phmm->emit_M_e[sb[j]][sa[i]] = logsum(phmm->emit_M_e[sb[j]][sa[i]], fM[i][j] + bM[i][j] - score);
                }


                //X[i][j] =                 M[i-1][j] + MX;
                eGPO = logsum(eGPO, fM[i-1][j] + MX + bX[i][j] + phmm->emit_background[sa[i]] - score);
                //X[i][j] = logsum(X[i][j], X[i-1][j] + XX);
                //eGPE = logsum(eGPE, fX[i-1][j] + XX + bX[i][j] + phmm->emit_background[sa[i]] - score);

                //X[i][j] += phmm->emit_background[sa[i]];
                phmm->emit_background_e[sa[i]] = logsum(phmm->emit_background_e[sa[i]], fX[i][j]+ bX[i][j]-score);

                //Y[i][j] =                 M[i][j-1] + MY;
                eGPO = logsum(eGPO, fM[i][j-1] + MY + bY[i][j] + phmm->emit_background[sb[j]] - score);
                //Y[i][j] = logsum(Y[i][j], Y[i][j-1] + YY);
                //eGPE = logsum(eGPE, fY[i][j-1] + YY + bY[i][j] + phmm->emit_background[sb[j]] - score);
                //Y[i][j] += phmm->emit_background[sb[j]];
                phmm->emit_background_e[sb[j]] = logsum(phmm->emit_background_e[sb[j]], fY[i][j] + bY[i][j]-score);

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

        /*fprintf(stdout,"1st round:\n");
        fprintf(stdout,"%f\tMM\n", scaledprob2prob(phmm->transition[INDEXMM]));
        fprintf(stdout,"%f\tGPO\n",scaledprob2prob(phmm->transition[INDEXGPO]));
        fprintf(stdout,"%f\tGPE\n",scaledprob2prob(phmm->transition[INDEXGPE]));
        fprintf(stdout,"%f\tTM\n",scaledprob2prob(phmm->transition[INDEXTM]));
        */
        /* add in tau */
        sum = phmm->transition[INDEXMM];
        sum = logsum(sum, phmm->transition[INDEXGPO]);
        sum = logsum(sum, phmm->transition[INDEXGPO]);
        sum = logsum(sum, phmm->tau);

        //exit(0);

        phmm->transition[INDEXMM] = phmm->transition[INDEXMM] - sum;
        phmm->transition[INDEXGPO] = phmm->transition[INDEXGPO] - sum;


        /*fprintf(stdout,"2nd round:\n");
        fprintf(stdout,"%f\tMM\n", scaledprob2prob(phmm->transition[INDEXMM]));
        fprintf(stdout,"%f\tGPO\n",scaledprob2prob(phmm->transition[INDEXGPO]));
        fprintf(stdout,"%f\tGPE\n",scaledprob2prob(phmm->transition[INDEXGPE]));
        fprintf(stdout,"%f\tTM\n",scaledprob2prob(phmm->transition[INDEXTM]));
        exit(0);*/
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
                //fprintf(stdout,"%d\t%f\n",i,scaledprob2prob(phmm->emit_background[i]));
        }
        sum = prob2scaledprob(0.0f);
        for(i = 0; i < L;i++){
                sum = logsum(sum, phmm->emit_background[i]);
        }
        //fprintf(stdout,"sanity:background: %f\n", scaledprob2prob(sum));

        sum = prob2scaledprob(0.0f);
        for(i = 0; i < L;i++){
                for(j = 0; j <= i;j++){
//                        fprintf(stdout," %0.2f",phmm->emit_M_e[i][j]);
                        sum = logsum(sum, phmm->emit_M_e[i][j]);
                }
                //fprintf(stdout,"\n");
        }

        //fprintf(stdout,"\n");fprintf(stdout,"SUM:%f\n",sum);

        for(i = 0; i < L;i++){
                for(j = 0; j <= i;j++){
                        phmm->emit_M[i][j] = phmm->emit_M_e[i][j] - sum;
                        phmm->emit_M[j][i] = phmm->emit_M[i][j];
                        //              fprintf(stdout," %0.2f",scaledprob2prob(phmm->emit_M[i][j])*100.0);
                        //fprintf(stdout,"%0.2f ",scaledprob2prob(phmm->emit_M[i][j]));

                }
//                fprintf(stdout,"\n");

        }
        //      fprintf(stdout,"\n");
        //    exit(0);
//fprintf(stdout,"\n");

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
        float alphabet = 23;
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
                        //      fprintf(stdout," %0.1f", phmm->emit_M[i][j]);
                }
                //      fprintf(stdout,"\n");

        }
        //fprintf(stdout,"\n");

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

        phmm->eta = prob2scaledprob(0.997738);
        //fprintf(stdout,"sanity: X/Y: %f\n", scaledprob2prob(sum));
        //exit(0);
        return OK;
}


int clear_phmm_e(struct phmm* phmm)
{
        int i,j;

        phmm->transition_e[INDEXMM] = prob2scaledprob(0.0f);
        phmm->transition_e[INDEXGPO] = prob2scaledprob(0.0f);
        phmm->transition_e[INDEXGPE] = prob2scaledprob(0.0f);
        phmm->transition_e[INDEXTM] = prob2scaledprob(0.0f);

        for(i = 0;i < phmm->L;i++){
                phmm->emit_background_e[i] = prob2scaledprob(0.0f);

                for(j = 0; j < phmm->L;j++){
                        phmm->emit_M_e[i][j] = prob2scaledprob(0.0f);
                }
        }

        return OK;
}


int add_pseudocounts(struct phmm* phmm, float w)
{
        int i,j;

        float prior_back[21] = {
-2.462039,
-4.290951,
-2.861703,
-2.670287,
-3.242155,
-2.542811,
-3.743423,
-2.809374,
-2.829935,
-2.397597,
-3.819814,
-3.188450,
-3.071127,
-3.347402,
-2.974614,
-2.837042,
-2.945077,
-2.657298,
-4.388594,
-3.402006,
-8.457309,
0.000000};
float prior_m[21][21] = {
{-3.177299,-6.191718,-5.419304,-4.957334,-5.935829,-4.777886,-6.339544,-5.268996,-5.102014,-4.902290,-6.165040,-5.657795,-5.499347,-5.647688,-5.481493,-4.709463,-5.124226,-4.784890,-7.212927,-6.058960,-10.571393},
{-6.191718,-4.954772,-7.795046,-7.669646,-7.500984,-7.225438,-8.290624,-6.950052,-7.753595,-6.616895,-7.820372,-7.657623,-7.697898,-8.023582,-7.629232,-6.917443,-6.908076,-6.411755,-9.012855,-7.828272,-12.526193},
{-5.419304,-7.795046,-3.489306,-4.657373,-6.812158,-5.450961,-6.417383,-6.554618,-5.414551,-6.006003,-7.269661,-5.261111,-5.965141,-5.845256,-5.910008,-5.426257,-5.781932,-6.281671,-7.790838,-6.722013,-10.435695},
{-4.957334,-7.669646,-4.657373,-3.418172,-6.538280,-5.477660,-6.250037,-6.091206,-4.851527,-5.571180,-6.770990,-5.559409,-5.684911,-5.237435,-5.387835,-5.317072,-5.523394,-5.787795,-7.587077,-6.502786,-10.271364},
{-5.935829,-7.500984,-6.812158,-6.538280,-3.972199,-6.495368,-6.768459,-5.554402,-6.665441,-4.862917,-6.374017,-6.843914,-6.674985,-6.922089,-6.711258,-6.457205,-6.383111,-5.609489,-6.553840,-5.220626,-11.429184},
{-4.777886,-7.225438,-5.450961,-5.477660,-6.495368,-2.764877,-6.624538,-6.330297,-5.525306,-5.821952,-7.006146,-5.458260,-5.848382,-6.147118,-5.940452,-5.205796,-5.740575,-6.053668,-7.598521,-6.578559,-11.280321},
{-6.339544,-8.290624,-6.417383,-6.250037,-6.768459,-6.624538,-4.554317,-6.998195,-6.319529,-6.320387,-7.654377,-6.328598,-7.012725,-6.578886,-6.306406,-6.527787,-6.666181,-6.706566,-7.735720,-6.364377,-12.616707},
{-5.268996,-6.950052,-6.554618,-6.091206,-5.554402,-6.330297,-6.998195,-3.582603,-6.108110,-4.185634,-5.885242,-6.636443,-6.346604,-6.528688,-6.200615,-6.139034,-5.708130,-4.131721,-7.280789,-6.212640,-10.753851},
{-5.102014,-7.753595,-5.414551,-4.851527,-6.665441,-5.525306,-6.319529,-6.108110,-3.773652,-5.540370,-6.757301,-5.643952,-5.906023,-5.433663,-4.837247,-5.491502,-5.600276,-5.781867,-7.720065,-6.526969,-10.771988},
{-4.902290,-6.616895,-6.006003,-5.571180,-4.862917,-5.821952,-6.320387,-4.185634,-5.540370,-2.978001,-5.135827,-6.000395,-5.884718,-5.860362,-5.641406,-5.685140,-5.445278,-4.384289,-6.661456,-5.560313,-10.473549},
{-6.165040,-7.820372,-7.269661,-6.770990,-6.374017,-7.006146,-7.654377,-5.885242,-6.757301,-5.135827,-5.060910,-7.227395,-7.310930,-6.965243,-6.916307,-6.797553,-6.628583,-5.983461,-7.952424,-6.960117,-11.444551},
{-5.657795,-7.657623,-5.261111,-5.559409,-6.843914,-5.458260,-6.328598,-6.636443,-5.643952,-6.000395,-7.227395,-4.134179,-6.332960,-6.149282,-5.957467,-5.540346,-5.865796,-6.339316,-7.938090,-6.818196,-11.043715},
{-5.499347,-7.697898,-5.965141,-5.684911,-6.674985,-5.848382,-7.012725,-6.346604,-5.906023,-5.884718,-7.310930,-6.332960,-3.552121,-6.486493,-6.250751,-5.781425,-6.103059,-6.007130,-8.033943,-6.921557,-12.243330},
{-5.647688,-8.023582,-5.845256,-5.237435,-6.922089,-6.147118,-6.578886,-6.528688,-5.433663,-5.860362,-6.965243,-6.149282,-6.486493,-4.579259,-5.716352,-5.968196,-6.113530,-6.251965,-7.978458,-6.835013,-10.446450},
{-5.481493,-7.629232,-5.910008,-5.387835,-6.711258,-5.940452,-6.306406,-6.200615,-4.837247,-5.641406,-6.916307,-5.957467,-6.250751,-5.716352,-3.633339,-5.783788,-5.878532,-5.973818,-7.650940,-6.450018,-11.836658},
{-4.709463,-6.917443,-5.426257,-5.317072,-6.457205,-5.205796,-6.527787,-6.139034,-5.491502,-5.685140,-6.797553,-5.540346,-5.781425,-5.968196,-5.783788,-3.884545,-4.970873,-5.766789,-7.519829,-6.485903,-10.788344},
{-5.124226,-6.908076,-5.781932,-5.523394,-6.383111,-5.740575,-6.666181,-5.708130,-5.600276,-5.445278,-6.628583,-5.865796,-6.103059,-6.113530,-5.878532,-4.970873,-3.849776,-5.225908,-7.636603,-6.530021,-11.117906},
{-4.784890,-6.411755,-6.281671,-5.787795,-5.609489,-6.053668,-6.706566,-4.131721,-5.781867,-4.384289,-5.983461,-6.339316,-6.007130,-6.251965,-5.973818,-5.766789,-5.225908,-3.417968,-7.217449,-6.098798,-10.813608},
{-7.212927,-9.012855,-7.790838,-7.587077,-6.553840,-7.598521,-7.735720,-7.280789,-7.720065,-6.661456,-7.952424,-7.938090,-8.033943,-7.978458,-7.650940,-7.519829,-7.636603,-7.217449,-4.919579,-6.558340,-13.507022},
{-6.058960,-7.828272,-6.722013,-6.502786,-5.220626,-6.578559,-6.364377,-6.212640,-6.526969,-5.560313,-6.960117,-6.818196,-6.921557,-6.835013,-6.450018,-6.485903,-6.530021,-6.098798,-6.558340,-4.059473,-12.064638},
{-10.571393,-12.526193,-10.435695,-10.271364,-11.429184,-11.280321,-12.616707,-10.753851,-10.771988,-10.473549,-11.444551,-11.043715,-12.243330,-10.446450,-11.836658,-10.788344,-11.117906,-10.813608,-13.507022,-12.064638,-12.004370},
};
float prior_MM = -0.069129;
float prior_GPO = -3.438407;
float prior_GPE = -0.099126;
float prior_TM = -2.388040;

        phmm->tau = prob2scaledprob(0.002569);
        phmm->eta = prob2scaledprob(0.998702);


        w = prob2scaledprob(w);

        phmm->transition_e[INDEXMM] = w + prior_MM;
        phmm->transition_e[INDEXGPO] = w + prob2scaledprob(scaledprob2prob(prior_GPO) /2.0 );
        phmm->transition_e[INDEXGPE] = w+ prior_GPE;
        phmm->transition_e[INDEXTM] = w+ prior_TM;

        /*fprintf(stdout,"%f\tMM\n", scaledprob2prob(phmm->transition_e[INDEXMM]));
        fprintf(stdout,"%f\tGPO\n",scaledprob2prob(phmm->transition_e[INDEXGPO]));
        fprintf(stdout,"%f\tGPE\n",scaledprob2prob(phmm->transition_e[INDEXGPE]));
        fprintf(stdout,"%f\tTM\n",scaledprob2prob(phmm->transition_e[INDEXTM]));
        */
        for(i = 0;i < phmm->L;i++){
                phmm->emit_background_e[i] = w + prior_back[i];
                for(j = 0; j < phmm->L;j++){
                        phmm->emit_M_e[i][j] = w + prior_m[i][j];

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

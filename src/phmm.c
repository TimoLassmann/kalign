

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
                        if(sa[i] > sb[j]){
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
                if(sa[i] > sb[j]){
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
/*                for(j = 0; j < L;j++){
                        fprintf(stdout,"%0.2f ",scaledprob2prob(phmm->emit_M[i][j]));
                }
                fprintf(stdout,"\n");*/
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

        phmm->transition_e[INDEXMM] = prob2scaledprob(1.0f);
        phmm->transition_e[INDEXGPO] = prob2scaledprob(1.0f);
        phmm->transition_e[INDEXGPE] = prob2scaledprob(1.0f);
        phmm->transition_e[INDEXTM] = prob2scaledprob(1.0f);

        for(i = 0;i < phmm->L;i++){
                phmm->emit_background_e[i] = prob2scaledprob(1.0f);

                for(j = 0; j < phmm->L;j++){
                        phmm->emit_M_e[i][j] = prob2scaledprob(1.0f);
                }
        }

        return OK;
}


int add_pseudocounts(struct phmm* phmm, float w)
{
        int i,j;
        float prior_back[23] = {
                -2.457247,
                -10.054303,
                -4.311800,
                -2.851731,
                -2.667678,
                -3.232455,
                -2.529456,
                -3.743798,
                -2.786486,
                -2.834563,
                -2.392979,
                -3.833805,
                -3.192592,
                -3.099623,
                -3.368073,
                -2.978477,
                -2.873273,
                -2.952438,
                -2.645252,
                -4.387171,
                -9.101513,
                -3.388789,
                -9.740719};
        float prior_m[23][23] = {
                {-3.586943,-12.187064,-5.908188,-5.136782,-4.674923,-5.655868,-4.494801,-6.058004,-4.986310,-4.822314,-4.620924,-5.881766,-5.376194,-5.217359,-5.365639,-5.200409,-4.426410,-4.840773,-4.501764,-6.933658,-10.877686,-5.776690,-11.500774},
                {-12.187064,-13.902116,-16.662125,-10.763325,-11.580721,-14.320320,-12.475506,-14.037457,-13.927758,-11.941842,-13.817216,-14.128428,-11.338140,-13.627172,-12.811978,-13.580215,-11.845075,-12.323528,-14.037457,-15.438350,-18.271563,-13.287957,-13.805655},
                {-5.908188,-16.662125,-5.362725,-7.511832,-7.387384,-7.219388,-6.943739,-8.008061,-6.666524,-7.475955,-6.333887,-7.536167,-7.374676,-7.415203,-7.741576,-7.347912,-6.634334,-6.624457,-6.127356,-8.729117,-12.371666,-7.546711,-14.421416},
                {-5.136782,-10.763325,-7.511832,-3.897418,-4.373279,-6.532986,-5.167504,-6.133904,-6.271998,-5.132123,-5.724230,-6.986028,-4.976907,-5.682599,-5.561400,-5.626865,-5.142759,-5.498054,-5.998865,-7.508567,-11.767275,-6.438745,-11.500774},
                {-4.674923,-11.580721,-7.387384,-4.373279,-3.826822,-6.261415,-5.197017,-5.967257,-5.809501,-4.569603,-5.290121,-6.487093,-5.277190,-5.403094,-4.953707,-5.104881,-5.034592,-5.240110,-5.505589,-7.306318,-11.421437,-6.220421,-10.568105},
                {-5.655868,-14.320320,-7.219388,-6.532986,-6.261415,-4.381587,-6.214397,-6.486064,-5.271128,-6.387944,-4.579815,-6.092168,-6.563507,-6.395155,-6.641015,-6.430450,-6.179578,-6.101191,-5.326118,-6.272281,-11.222308,-4.936403,-14.487373},
                {-4.494801,-12.475506,-6.943739,-5.167504,-5.197017,-6.214397,-3.173316,-6.343441,-6.052318,-5.246617,-5.543826,-6.722393,-5.175848,-5.569360,-5.866119,-5.662332,-4.922862,-5.457891,-5.772978,-7.318024,-11.684013,-6.297417,-12.297754},
                {-6.058004,-14.037457,-8.008061,-6.133904,-5.967257,-6.486064,-6.343441,-4.962453,-6.716211,-6.038401,-6.039263,-7.370851,-6.045585,-6.731568,-6.295930,-6.023299,-6.245401,-6.383391,-6.423738,-7.451925,-12.943687,-6.080437,-13.608125},
                {-4.986310,-13.927758,-6.666524,-6.271998,-5.809501,-5.271128,-6.052318,-6.716211,-3.990818,-5.828819,-3.901964,-5.601018,-6.356064,-6.063987,-6.246509,-5.919146,-5.857242,-5.424184,-3.847097,-6.997848,-10.530030,-5.929418,-14.008883},
                {-4.822314,-11.941842,-7.475955,-5.132123,-4.569603,-6.387944,-5.246617,-6.038401,-5.828819,-4.183959,-5.262429,-6.474452,-5.363887,-5.628259,-5.151264,-4.554500,-5.212001,-5.319683,-5.501311,-7.440450,-11.439610,-6.245677,-11.460320},
                {-4.620924,-13.817216,-6.333887,-5.724230,-5.290121,-4.579815,-5.543826,-6.039263,-3.901964,-5.262429,-3.386889,-4.851864,-5.720275,-5.603582,-5.579082,-5.360677,-5.404397,-5.162668,-4.100563,-6.379024,-10.270208,-5.276715,-13.165618},
                {-5.881766,-14.128428,-7.536167,-6.986028,-6.487093,-6.092168,-6.722393,-7.370851,-5.601018,-6.474452,-4.851864,-5.469102,-6.943968,-7.029148,-6.681353,-6.634069,-6.514391,-6.344557,-5.699755,-7.668524,-11.261251,-6.676229,-14.211121},
                {-5.376194,-11.338140,-7.374676,-4.976907,-5.277190,-6.563507,-5.175848,-6.045585,-6.356064,-5.363887,-5.720275,-6.943968,-4.542977,-6.051548,-5.866882,-5.675854,-5.257514,-5.582679,-6.058252,-7.656034,-12.006262,-6.535695,-12.633209},
                {-5.217359,-13.627172,-7.415203,-5.682599,-5.403094,-6.395155,-5.569360,-6.731568,-6.063987,-5.628259,-5.603582,-7.029148,-6.051548,-3.961577,-6.205253,-5.972292,-5.500381,-5.821604,-5.725087,-7.754730,-12.914977,-6.641433,-12.799293},
                {-5.365639,-12.811978,-7.741576,-5.561400,-4.953707,-6.641015,-5.866119,-6.295930,-6.246509,-5.151264,-5.579082,-6.681353,-5.866882,-6.205253,-4.987673,-5.433273,-5.685860,-5.830549,-5.969690,-7.698121,-12.011982,-6.552836,-10.419514},
                {-5.200409,-13.580215,-7.347912,-5.626865,-5.104881,-6.430450,-5.662332,-6.023299,-5.919146,-4.554500,-5.360677,-6.634069,-5.675854,-5.972292,-5.433273,-4.042877,-5.502720,-5.595900,-5.693975,-7.369413,-11.922424,-6.167623,-13.267617},
                {-4.426410,-11.845075,-6.634334,-5.142759,-5.034592,-6.179578,-4.922862,-6.245401,-5.857242,-5.212001,-5.404397,-6.514391,-5.257514,-5.500381,-5.685860,-5.502720,-4.293262,-4.686790,-5.484488,-7.238175,-11.165777,-6.203707,-12.002467},
                {-4.840773,-12.323528,-6.624457,-5.498054,-5.240110,-6.101191,-5.457891,-6.383391,-5.424184,-5.319683,-5.162668,-6.344557,-5.582679,-5.821604,-5.830549,-5.595900,-4.686790,-4.257891,-4.941537,-7.354350,-11.499628,-6.247732,-12.169004},
                {-4.501764,-14.037457,-6.127356,-5.998865,-5.505589,-5.326118,-5.772978,-6.423738,-3.847097,-5.501311,-4.100563,-5.699755,-6.058252,-5.725087,-5.969690,-5.693975,-5.484488,-4.941537,-3.826304,-6.935471,-10.674167,-5.815472,-12.766232},
                {-6.933658,-15.438350,-8.729117,-7.508567,-7.306318,-6.272281,-7.318024,-7.451925,-6.997848,-7.440450,-6.379024,-7.668524,-7.656034,-7.754730,-7.698121,-7.369413,-7.238175,-7.354350,-6.935471,-5.328221,-13.337090,-6.273943,-17.578417},
                {-10.877686,-18.271563,-12.371666,-11.767275,-11.421437,-11.222308,-11.684013,-12.943687,-10.530030,-11.439610,-10.270208,-11.261251,-12.006262,-12.914977,-12.011982,-11.922424,-11.165777,-11.499628,-10.674167,-13.337090,-16.885269,-12.108249,-18.271563},
                {-5.776690,-13.287957,-7.546711,-6.438745,-6.220421,-4.936403,-6.297417,-6.080437,-5.929418,-6.245677,-5.276715,-6.676229,-6.535695,-6.641433,-6.552836,-6.167623,-6.203707,-6.247732,-5.815472,-6.273943,-12.108249,-4.467856,-14.582684},
                {-11.500774,-13.805655,-14.421416,-11.500774,-10.568105,-14.487373,-12.297754,-13.608125,-14.008883,-11.460320,-13.165618,-14.211121,-12.633209,-12.799293,-10.419514,-13.267617,-12.002467,-12.169004,-12.766232,-17.578417,-18.271563,-14.582684,-13.056627}
        };
        float prior_MM = -0.062269;
        float prior_GPO = -3.552166;
        float prior_GPE = -0.256485;
        float prior_TM = -1.499741;

        w = prob2scaledprob(w);

        phmm->transition_e[INDEXMM] = w + prior_MM;
        phmm->transition_e[INDEXGPO] = w + prior_GPO ;
        phmm->transition_e[INDEXGPE] = w+ prior_GPE;
        phmm->transition_e[INDEXTM] = w+ prior_TM;

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

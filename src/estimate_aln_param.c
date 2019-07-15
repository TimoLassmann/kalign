/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019 Timo Lassmann

    This file is part of kalign.

    Kalign is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#include "estimate_aln_param.h"

#include "phmm.h"

#include "pick_anchor.h"

int get_background(struct alignment*aln, struct phmm* phmm);
int create_subm(struct aln_param* ap,struct phmm* phmm);


int estimate_aln_param(struct alignment*aln, struct aln_param*ap)
{
        struct phmm* phmm = NULL;
        uint8_t* seqa = NULL;
        uint8_t* seqb = NULL;

        int* anchors = NULL;
        int num_anchor = 0;
        int len_a,len_b;

        int i,j;
        ASSERT(aln != NULL, "No alignment.");
        ASSERT(ap != NULL, "No parameters.");



        RUNP(anchors = pick_anchor(aln, &num_anchor));
        i = aln->sl[anchors[0]] + 2;
        RUNP(phmm = alloc_phmm(i));
        phmm->L = aln->L;

        RUN(get_background(aln,phmm));
        //RUN(simple_init(phmm));

        RUN(add_pseudocounts(phmm, 1));
        RUN(re_estimate(phmm));
        RUN(phmm_transitions(phmm));

        RUN(clear_phmm_e(phmm));
        RUN(add_pseudocounts(phmm, 1));

        /* set  eta  */
        LOG_MSG("%f eta", scaledprob2prob(phmm->eta));
        phmm->eta = 0.0;
        for(i = 0; i < aln->numseq;i++){
                phmm->eta += aln->sl[i];

        }

        phmm->eta = 1.0  - 1.0 / (phmm->eta / aln->numseq);
        phmm->eta = prob2scaledprob(phmm->eta);
        LOG_MSG("%f eta", scaledprob2prob(phmm->eta));
//exit(0);

        //DECLARE_TIMER(t1);
        for(int iter = 0;iter < 10;iter++){
                //START_TIMER(t1);
                for(i = 0; i < aln->numseq;i++){
                        seqa = aln->s[i];
                        len_a = aln->sl[i];
                        for(j = 0; j < num_anchor;j++){
                                seqb = aln->s[anchors[j]];
                                len_b = aln->sl[anchors[j]];
                                //fprintf(stdout,"%d vs %d\n", i, anchors[j]);
                                RUN(forward_phmm(phmm, seqa, seqb, len_a, len_b));
                                RUN(backward_phmm(phmm, seqa, seqb, len_a, len_b));
                                RUN(collect_phmm(phmm, seqa, seqb, len_a, len_b));

                        }
                }
                RUN(re_estimate(phmm));
                RUN(add_pseudocounts(phmm, 1));// * aln->numseq * num_anchor));
                //STOP_TIMER(t1);
                //fprintf(stdout,"%d %f sec\n", iter, GET_TIMING(t1));
                RUN(phmm_transitions(phmm));

        }

        RUN(phmm_transitions(phmm));
        RUN(create_subm(ap, phmm));
        //print_phmm(phmm, len_a, len_b);
        free_phmm(phmm);




        MFREE(anchors);

        return OK;
ERROR:
        return FAIL;
}



int get_background(struct alignment*aln, struct phmm* phmm)
{

        uint8_t* seq = NULL;
        int len;

        int i,j;
        float sum;
        for(i = 0; i < aln->L;i++){
                phmm->emit_background[i] = 0.0f;
        }
        /* get background */
        for(i =0; i <aln->numseq;i++){
                seq = aln->s[i];
                len = aln->sl[i];
                for(j = 0; j < len;j++){
                        phmm->emit_background[seq[j]]++;
                }
        }
        sum = 0.0;
        for(i = 0; i < aln->L;i++){
                sum+= phmm->emit_background[i];
        }

        for(i = 0; i < aln->L;i++){
                phmm->emit_background[i] = prob2scaledprob(phmm->emit_background[i] / sum);
        }
        return OK;
}

int create_subm(struct aln_param* ap,struct phmm* phmm)
{
        float tau = 0.003055;
        float eta = 0.998702;
        int i,j;
        double sum;

        ASSERT(phmm!= NULL, "No phmm");
        ASSERT(ap != NULL, "No param");
        sum = prob2scaledprob(0.0);

        for(i = 0; i < phmm->L;i++){

                sum = logsum(sum, phmm->emit_background[i]);
        }
        for(i = 0; i < phmm->L;i++){
                phmm->emit_background[i] = scaledprob2prob(phmm->emit_background[i] - sum);

        }
        sum = prob2scaledprob(0.0f);
        for(i = 0; i < phmm->L;i++){
                for(j = 0; j <= i;j++){
                        sum = logsum(sum, phmm->emit_M[i][j]);
                }
        }
        for(i = 0; i < phmm->L;i++){
                for(j = 0; j <= i;j++){
                        phmm->emit_M[i][j] = scaledprob2prob(phmm->emit_M[i][j] -   sum);
                        phmm->emit_M[j][i] = phmm->emit_M[i][j];
                }
        }


        sum = phmm->transition[INDEXMM];
        sum =  logsum(sum,phmm->transition[INDEXGPO]);
        sum =  logsum(sum,phmm->transition[INDEXGPO]);

        phmm->transition[INDEXMM] = scaledprob2prob(phmm->transition[INDEXMM] - sum);
        phmm->transition[INDEXGPO] = scaledprob2prob(phmm->transition[INDEXGPO] - sum);


        /* add in tau */
        sum = phmm->transition[INDEXMM];
        sum += phmm->transition[INDEXGPO];
        sum += phmm->transition[INDEXGPO];
        sum += tau;

        phmm->transition[INDEXMM] = phmm->transition[INDEXMM]  / sum;
        phmm->transition[INDEXGPO] = phmm->transition[INDEXGPO] / sum;

        sum = phmm->transition[INDEXGPE];
        sum = logsum(sum, phmm->transition[INDEXTM]);


        phmm->transition[INDEXGPE] = scaledprob2prob(phmm->transition[INDEXGPE] - sum);
        phmm->transition[INDEXTM] = scaledprob2prob(phmm->transition[INDEXTM] - sum);
        //ap->GPE = ap->GPE / sum;
        //ap->TM = ap->TM / sum;

        sum = phmm->transition[INDEXGPE];
        sum += phmm->transition[INDEXTM];
        sum += tau;
        phmm->transition[INDEXGPE] = phmm->transition[INDEXGPE] /sum;

        phmm->transition[INDEXTM] = phmm->transition[INDEXTM] /sum;



        /* set subm and penalties  */

        for(i = 0; i < phmm->L;i++){
                //fprintf(stdout,"%d",i);
                for(j = 0; j <= i;j++){
                        sum = log2(phmm->emit_M[i][j] / ( phmm->emit_background[i] * phmm->emit_background[j])) + log2(phmm->transition[INDEXMM]/((eta)*(eta)));
                        ap->subm[i][j] = sum*2.0f;
                        ap->subm[j][i] = sum*2.0f;

                        //fprintf(stdout," %f,", sum);
                }
                //fprintf(stdout,"\n");
        }
        //fprintf(stdout,"\n");
        /*fprintf(stdout,"%f\tMM\n", ap->MM);
        fprintf(stdout,"%f\tGPO\n", ap->GPO);
        fprintf(stdout,"%f\tGPE\n", ap->GPE);
        fprintf(stdout,"%f\tTM\n", ap->TM);*/
        /* taushould be 1/ average length */
        //fprintf(stdout,"%f\ttau\n", tau);
        //fprintf(stdout,"%f\teta\n", eta);
        sum = 0.0;

        sum = -1.0 * log2( (phmm->transition[INDEXGPO] * phmm->transition[INDEXTM]) / ((eta) * phmm->transition[INDEXMM]));
        //fprintf(stdout,"%f;\n", sum);
        ap->gpo = sum / 2.0f;

        sum = -1.0 *log2(phmm->transition[INDEXGPE]  /(1.0 - tau));
        //fprintf(stdout,"ap->gpe =  %f;\n", sum);
        ap->gpe = sum;
        //fprintf(stdout,"ap->tgpe =  %f;\n", 0.0);
        //

        return OK;
ERROR:
        return FAIL;
}


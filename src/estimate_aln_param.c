
#include "estimate_aln_param.h"

#include "phmm.h"

#include "pick_anchor.h"


int create_subm(struct aln_param* ap,struct phmm* phmm);

int estimate_aln_param(struct alignment*aln, struct aln_param*ap)
{
        struct phmm* phmm = NULL;
        int* anchors = NULL;
        int num_anchor = 0;

        int i,j;
        ASSERT(aln != NULL, "No alignment.");
        ASSERT(ap != NULL, "No parameters.");



        RUNP(anchors = pick_anchor(aln, &num_anchor));
        i = aln->sl[anchors[0]] + 2;
        RUNP(phmm = alloc_phmm(i));

        RUN(simple_init(phmm));

        RUN(clear_phmm_e(phmm));
        RUN(add_pseudocounts(phmm, 1 * aln->numseq * num_anchor));
        RUN(phmm_transitions(phmm));

        uint8_t* seqa = NULL;
        uint8_t* seqb = NULL;
        int len_a,len_b;
        DECLARE_TIMER(t1);
        for(int iter = 0;iter < 5;iter++){
                START_TIMER(t1);
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
                RUN(add_pseudocounts(phmm, 1 * aln->numseq * num_anchor));
                STOP_TIMER(t1);
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


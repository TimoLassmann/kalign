#include "tldevel.h"
#include "pick_anchor.h"
#include "sequence_distance.h"

#define IDATA_IMPORT
#include "idata.h"

int get_internal_data(struct msa* msa, struct aln_param* ap,double** out,int* n_out)
{
    //struct drand48_data rand Buffer;
        struct node* root = NULL;
        double* stats = NULL;

        int* samples = NULL;
        int* anchors = NULL;
        int num_anchors;
        int numseq;

        float d;

        int i,j;

        ASSERT(msa != NULL, "No alignment.");
        //ASSERT(param != NULL, "No input parameters.");
        ASSERT(ap != NULL, "No alignment parameters.");


        MMALLOC(stats, sizeof(double) * 13);


        for(i= 0; i < 13;i++){
                stats[i] = 0.0;
        }
        numseq = msa->numseq;



        RUNP(anchors = pick_anchor(msa, &num_anchors));

        for(i = 0; i < numseq;i++){
                stats[0]++;
                stats[1] += msa->sequences[i]->len;
                stats[2] += msa->sequences[i]->len * msa->sequences[i]->len;
                for(j = 0;j < num_anchors;j++){

                        d = calc_distance(msa->sequences[i]->s,
                                          msa->sequences[anchors[j]]->s,
                                          msa->sequences[i]->len,
                                          msa->sequences[anchors[j]]->len,
                                          msa->L);
                        //LOG_MSG("%d %d %f", i,j,d);
                                //dm[i][j] = dist;
                        if(d > 200){
                                d = 199.9F;
                        }
                        stats[3+  (int) ( d / 20.0)]++;

                }
        }

        d  =  sqrt ( (stats[0] * stats[2] -  pow(stats[1], 2.0)) /  (stats[0] * ( stats[0] - 1.0)));
        stats[2] = d;
        stats[1] = stats[1] /stats[0];

        d = 0.0;
        for(i = 3; i < 13;i++){
                d+= stats[i];
        }
        for(i = 3; i < 13;i++){
                stats[i] = stats[i] / d * 10.0;
        }

        *out = stats;
        *n_out = 13;
        return OK;
ERROR:
        return FAIL;
}

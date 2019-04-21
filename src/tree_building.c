#include "tree_building.h"

int build_tree(struct alignment* aln,struct parameters* param, struct aln_param* ap)
{
        float** dm = NULL;
        int* tree = NULL;
        int numseq;
        int a,b,c;
        int i;

        ASSERT(aln != NULL, "No alignment.");
        ASSERT(param != NULL, "No input parameters.");
        ASSERT(ap != NULL, ?"No alignment parameters.");

        tree = ap->tree;
        numseq = aln->numseq;

        /* if we don't build a tree fill and return  */
        if(param->ntree <= 1){
                tree[0] = 0;
                tree[1] = 1;

                c = numseq;
                tree[2] = c;
                a = 2;
                for ( i = 3; i < (numseq-1)*3;i+=3){
                        tree[i] = c;
                        tree[i+1] = a;
                        c++;
                        tree[i+2] = c;
                        a++;
                }
                return OK;
        }

        /* calculate distances  */

        if(param->ntree > 1){
                if(byg_start(param->distance,"pairclustalPAIRCLUSTAL") != -1){
                        if(byg_start(param->tree,"njNJ") != -1){
                                dm = protein_pairwise_alignment_distance(aln,dm,param,submatrix,1);
                        }else{
                                dm = protein_pairwise_alignment_distance(aln,dm,param,submatrix,0);
                        }
                }else if(byg_start("wu",param->alignment_type) != -1){
                        dm =  protein_wu_distance2(aln,dm,param);
                }else if(param->dna == 1){
                        if(byg_start(param->tree,"njNJ") != -1){
                                dm =  dna_distance(aln,dm,param,1);
                        }else{
                                dm =  dna_distance(aln,dm,param,0);
                        }
                }else{
                        if(byg_start(param->tree,"njNJ") != -1){
                                dm =  protein_wu_distance(aln,dm,param,1);
                        }else{
                                dm =  protein_wu_distance(aln,dm,param,0);
                        }
                }
                /*int j;
                  for (i = 0; i< numseq;i++){
                  for (j = 0; j< numseq;j++){
                  fprintf(stderr,"%f	",dm[i][j]);
                  }
                  fprintf(stderr,"\n");
                  }*/

                if(byg_start(param->tree,"njNJ") != -1){
                        tree2 = real_nj(dm,param->ntree);
                }else{
                        tree2 = real_upgma(dm,param->ntree);
                }
                if(param->print_tree){
                        print_tree(tree2,aln,param->print_tree);
                }
        }

        tree = malloc(sizeof(int)*(numseq*3+1));
        for ( i = 1; i < (numseq*3)+1;i++){
                tree[i] = 0;
        }
        tree[0] = 1;

        if(param->ntree < 2){
                tree[0] = 0;
                tree[1] = 1;

                c = numseq;
                tree[2] = c;
                a = 2;
                for ( i = 3; i < (numseq-1)*3;i+=3){
                        tree[i] = c;
                        tree[i+1] = a;
                        c++;
                        tree[i+2] = c;
                        a++;
                }
        }else if(param->ntree > 2){
                ntreeify(tree2,param->ntree);
        }else{
                tree = readtree(tree2,tree);
                for (i = 0; i < (numseq*3);i++){
                        tree[i] = tree[i+1];
                }
                free(tree2->links);
                free(tree2->internal_lables);
                free(tree2);
        }



        return OK;
}

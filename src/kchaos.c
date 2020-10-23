#include "tldevel.h"

#include "queue.h"


int treeify_samples(int*** tree,int* samples,int numseq, int chaos);

int main(int argc, char *argv[])
{

        int** tree = NULL;

        int* samples = NULL;
        int numseq = 7;
        MMALLOC(samples, sizeof(int) * numseq);



        treeify_samples(&tree, samples, numseq, 3);

        gfree(tree);
        MFREE(samples);
        return OK;
ERROR:
        return FAIL;
}




int treeify_samples(int*** tree,int* samples,int numseq, int chaos)
{
        int** t = NULL;
        int i,j,c;
        int work;
        int cur_node;
        int bin_index;
        ASSERT(samples != NULL, "No samples");

        if(*tree){
                t = *tree;
        }else{
                RUN(galloc(&t,numseq*2,chaos*2));
        }

        for(i = 0; i < numseq*2;i++){
                for(j = 0; j <  chaos*2;j++){
                        t[i][j] = 0;
                }
        }

        queue q = q_new();
        for(i = 0; i < numseq;i++){
                enqueue(q, i+1);
        }

        print_queue(q);

        work = 1;
        j= 0;
        cur_node = numseq;
        bin_index = numseq;
        while(work){
                print_queue(q);
                if (!dequeue(q, &work)) {
                        break;
                }
                t[cur_node][j] = work-1;
                if(j){
                        t[cur_node][j-1+chaos] = bin_index;
                        bin_index++;
                }
                j++;
                if(j == chaos){
                        enqueue(q, bin_index cur_node+1);
                        cur_node = bin_index;
                        j = 0;

                }
        }
        for(i = 0; i < numseq*2;i++){
                for(j = 0; j <  chaos*2;j++){
                        fprintf(stdout,"%d ",t[i][j]);
                }
                fprintf(stdout,"\n");
        }

        exit(0);

        /* c = numseq;             /\* Keep track of implied  *\/ */

        /* cur_node = 0; */
        /* j = 0; */
        /* for(i = 0; i < numseq;i++){ */

        /*         t[cur_node][j] = i; */

        /* } */

        return OK;
ERROR:
        return FAIL;
}

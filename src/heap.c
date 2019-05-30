#include "heap.h"

#include "rng.h"

int min (heap_t *h, int i, int j, int k);

int print_heap(heap_t *h, int num);
#ifdef ITEST

int main(int argc, char *argv[])
{
        struct rng_state* rng;

        int i,j;

        int heap_size = 100;


        RUNP(rng = init_rng(0));
        //lrand48_r(randBuffer, &r);
        //r = r % num_samples;
        heap_t* h = NULL;
        RUNP(h = create_heap(heap_size));
        for(i = 0; i < 88;i++){
                j = tl_random_int(rng,100);
                push_heap(h, i,j);

        }

        print_heap(h,heap_size);

        for(i = 0; i < h->len;i++){
                j = pop_heap(h);
                fprintf(stdout,"%d\t%d\n", i,j);
        }
        free_heap(h);
        MFREE(rng);
        return EXIT_SUCCESS;

ERROR:
        return EXIT_FAILURE;
}

#endif


heap_t *create_heap (int n)
{
        heap_t *h = NULL;
        int i;
        MMALLOC(h, sizeof(heap_t));
        h->len = 0;
        h->size = 0;
        h->data = NULL;
        h->index = NULL;
        h->prio = NULL;

        MMALLOC(h->data, sizeof(int) * (n+1));
        MMALLOC(h->prio, sizeof(int) * (n+1));
        MMALLOC(h->index, sizeof(int) * n);
        for(i = 0; i < n;i++){
                h->data[i] = 0;
                h->prio[i] = 0;
                h->index[i] = 0;
        }
        h->data[n] = 0;
        h->prio[n] = 0;

        return h;
ERROR:
        return NULL;
}

void free_heap(heap_t* h)
{
        if(h){
                MFREE(h->data);
                MFREE(h->prio);
                MFREE(h->index);
                MFREE(h);
        }
}


void push_heap(heap_t *h, int v, int p)
{
        int i = h->index[v] == 0 ? ++h->len : h->index[v];
        int j = i / 2;
        while (i > 1) {
                if (h->prio[j] < p){
                        break;
                }
                h->data[i] = h->data[j];
                h->prio[i] = h->prio[j];
                h->index[h->data[i]] = i;
                i = j;
                j = j / 2;
        }
        h->data[i] = v;
        h->prio[i] = p;
        h->index[v] = i;
}


int min(heap_t *h, int i, int j, int k)
{
        int m = i;
        if (j <= h->len && h->prio[j] < h->prio[m])
                m = j;
        if (k <= h->len && h->prio[k] < h->prio[m])
                m = k;
        return m;
}

int print_heap(heap_t *h, int num)
{

        int i;
        for(i = 0; i < num;i++){

                fprintf(stdout,"%d %d\n",h->prio[i], h->data[i]);
        }
        return OK;
}

int pop_heap(heap_t *h)
{
        int v = h->data[1];
        int i = 1;

        fprintf(stdout,"data: %d pri:%d\n",h->data[1],h->prio[1]);
        while (1) {
                int j = min(h, h->len, 2 * i, 2 * i + 1);
                if (j == h->len){
                        break;
                }
                h->data[i] = h->data[j];
                h->prio[i] = h->prio[j];
                h->index[h->data[i]] = i;
                i = j;
        }
        h->data[i] = h->data[h->len];
        h->prio[i] = h->prio[h->len];
        h->index[h->data[i]] = i;
        h->len--;
        return v;
}

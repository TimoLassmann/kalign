#include "tldevel.h"

#define  CORETRALIGN_IMPORT
#include "coretralign.h"


typedef struct aln_scheduler {
        aln_elem_queue* task_alloc;
        aln_elem_queue* task_queue;
        pthread_t* threads;
        int64_t* thread_id_map;
        int thread_id_idx;
        pthread_mutex_t lock;
        int n_threads;
} aln_scheduler;



int aln_scheduler_alloc(aln_scheduler **scheduler)
{
        aln_scheduler* n = NULL;
        return OK;
ERROR:
        return FAIL;
}



static int aln_elem_alloc(aln_elem **node);
static int aln_elem_free(aln_elem *n);

int queue_push(aln_elem_queue *q, aln_elem *n)
{
        if(q->head == NULL){
                q->head = n;
                q->tail = q->head;
                q->head->p = NULL;
        }else{
                n->p = q->head;
                q->head = n;
        }
        q->n++;
        return OK;
}

int queue_pop(aln_elem_queue *q, aln_elem **node)
{
        aln_elem * tmp = NULL;
        if(q->head == NULL){
                if(q->is_allocator){
                        /* RUN(eq_add_mem(q)); */
                }else{
                        *node = NULL;
                        return FAIL;
                }
        }

        *node = q->head;
        tmp = q->head;
        q->head = tmp->p;

        if(q->head == NULL){
                q->tail = NULL;
        }
        q->n--;

        return OK;
}

int queue_append(aln_elem_queue *q, aln_elem* n)
{

        if(q->head == NULL){
                q->head = n;
                q->tail = q->head;

        }else{
                q->tail->p = n;
                q->tail = q->tail->p;
        }
        q->tail->p = NULL;
        q->n++;
        return OK;
}


int queue_add_mem(aln_elem_queue *q)
{
        ASSERT(q->alloc_block_size != 0,"This queue is not an allocator!");

        for(int i = 0; i < q->alloc_block_size;i++){
                /* Alloc e_node  */
                aln_elem* n = NULL;
                aln_elem_alloc(&n);
                /* Add to queue  */
                queue_append(q, n);
        }
        return OK;
ERROR:

        return FAIL;
}

int queue_alloc(aln_elem_queue **queue, uint8_t is_allocator)
{
        aln_elem_queue* l = NULL;
        MMALLOC(l, sizeof(aln_elem_queue));
        l->head = NULL;
        l->tail = 0;
        l->n = 0;
        l->alloc_block_size = 0;

        l->is_allocator = is_allocator;
        if(is_allocator){
                l->alloc_block_size =  4096;
                RUN(queue_add_mem(l));
                /* RUN(eq_add_mem(l)); */
        }
        *queue = l;
        return OK;
ERROR:
        queue_free(l);
        return FAIL;
}

void queue_free(aln_elem_queue *l)
{
        if(l){
                aln_elem* n = NULL;
                while(l->n){
                        n = NULL;
                        queue_pop(l, &n);
                        aln_elem_free(n);
                }
                MFREE(l);
        }
}

int aln_elem_alloc(aln_elem **node)
{
        aln_elem* n = NULL;

        MMALLOC(n, sizeof(aln_elem));
        n->l = NULL;
        n->r = NULL;
        n->p = NULL;
        n->wait = 0;
        n->type = AE_TYPE_UNDEF;
        n->aln_mem = NULL;
        *node = n;
        return OK;
ERROR:
        aln_elem_free(n);
        return FAIL;
}

int aln_elem_free(aln_elem *n)
{
        if(n){
                MFREE(n);
        }
}

#include "tldevel.h"

#define  CORETRALIGN_IMPORT
#include "coretralign.h"

int aln_scheduler_get_tid(aln_scheduler *s)
{
        int i;
        int64_t tid = (int64_t) pthread_self();
        int ID = -1;
        aln_scheduler_lock(s);
        for(i = 0; i < s->thread_id_idx;i++){
                if(s->thread_id_map[i] == tid){
                        ID = i;
                        break;
                }
        }
        if(ID == -1){
                ID = s->thread_id_idx;
                s->thread_id_map[s->thread_id_idx] = tid;
                s->thread_id_idx++;
        }
        aln_scheduler_unlock(s);
        return ID;
}

int  aln_scheduler_lock(aln_scheduler* s)
{
        ASSERT(s != NULL, "No thread controll");
        /* pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL); */
        if(pthread_mutex_lock(&s->lock) != 0){
                ERROR_MSG("Can't get lock");
        }
        return OK;
ERROR:
        return FAIL;
}

int  aln_scheduler_trylock(aln_scheduler* s)
{
        return pthread_mutex_trylock(&s->lock);
}

int  aln_scheduler_unlock(aln_scheduler* s)
{
        ASSERT(s != NULL, "No thread controll");
        if(pthread_mutex_unlock(&s->lock) != 0){
                ERROR_MSG("Can't get lock");
        }
        /* pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL); */
        /* pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL); */
        /* pthread_testcancel(); */
        return OK;
ERROR:
        return FAIL;
}


int aln_scheduler_alloc(aln_scheduler **scheduler, int n_threads)
{
        aln_scheduler* n = NULL;

        MMALLOC(n, sizeof(aln_scheduler));
        n->task_alloc = NULL;
        n->task_queue = NULL;
        n->root = NULL;
        n->threads = NULL;
        n->thread_id_map = NULL;
        n->thread_id_idx = 0;
        pthread_mutex_init(&n->lock, NULL);
        n->n_threads = n_threads;

        MMALLOC(n->threads, sizeof(pthread_t) * n->n_threads);

        MMALLOC(n->thread_id_map, sizeof(int64_t) * n->n_threads);


        queue_alloc(&n->task_alloc, 1);
        queue_alloc(&n->task_queue, 0);



        *scheduler = n;
        return OK;
ERROR:
        return FAIL;
}

void aln_scheduler_free(aln_scheduler *n)
{
        if(n){
                if(n->task_alloc){
                        queue_free(n->task_alloc);
                }
                if(n->task_queue){
                        queue_free(n->task_queue);
                }
                if(n->threads){
                        MFREE(n->threads);
                }
                if(n->thread_id_map){
                        MFREE(n->thread_id_map);
                }
                MFREE(n);
        }

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

#ifndef CORETRALIGN_H
#define CORETRALIGN_H

#include <pthread.h>
#include <stdint.h>
#include "aln_mem.h"

#ifdef CORETRALIGN_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

/* Idea:
   construct hirschberg alignment via trees
   -> only because this would make threading easier to control
*/

typedef enum {
        AE_TYPE_UNDEF,
        AE_TYPE_ALN_BACK,
        AE_TYPE_ALN_FORWARD,
        AE_TYPE_SPLIT
        /* tree stuff here  */
}aln_elem_type;

typedef struct aln_elem aln_elem;
typedef struct aln_elem {
        aln_elem* l;           /* left OR backward in dyn*/
        aln_elem* r;           /* right OR backward in dyn*/
        aln_elem* p;           /* parent */
        uint8_t wait;
        aln_elem_type type;
        struct aln_mem* aln_mem;
} aln_elem;

typedef struct aln_elem_queue {
        aln_elem* head;
        aln_elem* tail;
        int alloc_block_size;
        int n;
        uint8_t is_allocator;
} aln_elem_queue;

typedef struct aln_scheduler {
        aln_elem_queue* task_alloc;
        aln_elem_queue* task_queue;
        aln_elem* root;
        pthread_t* threads;
        int64_t* thread_id_map;
        int thread_id_idx;
        pthread_mutex_t lock;
        int n_threads;
        double** dm;
        struct msa* msa;
} aln_scheduler;

EXTERN  int aln_scheduler_lock(aln_scheduler* s);
EXTERN  int aln_scheduler_trylock(aln_scheduler* s);
EXTERN  int aln_scheduler_unlock(aln_scheduler* s);

EXTERN  int aln_scheduler_alloc(aln_scheduler **scheduler, int n_threads);
EXTERN void aln_scheduler_free(aln_scheduler *n);

EXTERN  int queue_alloc(aln_elem_queue **queue, uint8_t is_allocator);
EXTERN void queue_free(aln_elem_queue *l);

#undef CORETRALIGN_IMPORT
#undef EXTERN

#endif

/* tl_scratch.c -- free thread-local scratch buffers when a thread exits. */
#include <pthread.h>
#include <stdlib.h>

#include "tl_scratch.h"

struct tl_scratch {
        void* p[TL_SCRATCH_SLOTS][TL_SCRATCH_PTRS];
};

static pthread_key_t tl_scratch_key;
static pthread_once_t tl_scratch_once = PTHREAD_ONCE_INIT;

/* runs when a thread that tracked scratch exits; touches only heap memory */
static void tl_scratch_free(void* v)
{
        struct tl_scratch* s = v;
        for(int i = 0; i < TL_SCRATCH_SLOTS; i++){
                for(int j = 0; j < TL_SCRATCH_PTRS; j++){
                        free(s->p[i][j]);
                }
        }
        free(s);
}

static void tl_scratch_init(void)
{
        if(pthread_key_create(&tl_scratch_key, tl_scratch_free)) abort();
}

void tl_scratch_track(int slot, void* const* ptrs, int n)
{
        struct tl_scratch* s;
        if(slot < 0 || slot >= TL_SCRATCH_SLOTS || n < 0 || n > TL_SCRATCH_PTRS) abort();
        pthread_once(&tl_scratch_once, tl_scratch_init);
        s = pthread_getspecific(tl_scratch_key);
        if(!s){
                s = calloc(1, sizeof(*s));
                if(!s || pthread_setspecific(tl_scratch_key, s)) abort();
        }
        for(int j = 0; j < TL_SCRATCH_PTRS; j++){
                s->p[slot][j] = j < n ? ptrs[j] : NULL;
        }
}

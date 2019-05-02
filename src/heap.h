#ifndef HEAP_H
#define HEAP_H

#include "global.h"


typedef struct {
        int *data;
        int *prio;
        int *index;
        int len;
        int size;
} heap_t;



heap_t *create_heap (int n);

void free_heap(heap_t* h);

void push_heap (heap_t *h, int v, int p);
int pop_heap (heap_t *h);


#endif

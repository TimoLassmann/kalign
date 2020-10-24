#ifndef QUEUE_H
#define QUEUE_H


#ifdef QUEUE_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

typedef struct node_t node_t, *node, *queue;

EXTERN queue q_new();
EXTERN void enqueue(queue q, int n);
EXTERN int dequeue(queue q, int *val);
EXTERN void print_queue(queue q);
EXTERN void free_queue(queue q);
#undef QUEUE_IMPORT
#undef EXTERN

#endif

#include "stdio.h"
#include "stdlib.h"

#define QUEUE_IMPORT
#include "queue.h"
struct node_t { int val; node prev, next; };

#define HEAD(q) q->prev
#define TAIL(q) q->next



queue q_new()
{
        node q = malloc(sizeof(node_t));
        q->next = q->prev = 0;
        return q;
}

int empty(queue q)
{
        return !HEAD(q);
}

void enqueue(queue q, int n)
{
        node nd = malloc(sizeof(node_t));
        nd->val = n;
        if (!HEAD(q)){
                HEAD(q) = nd;
        }
        nd->prev = TAIL(q);
        if (nd->prev){
                nd->prev->next = nd;
        }
        TAIL(q) = nd;
        nd->next = 0;
}

int dequeue(queue q, int *val)
{
        node tmp = HEAD(q);
        if (!tmp){
                return 0;
        }
        *val = tmp->val;

        HEAD(q) = tmp->next;
        if (TAIL(q) == tmp){
                TAIL(q) = 0;
        }
        free(tmp);
        return 1;
}


void print_queue(queue q)
{
        node tmp = HEAD(q);
        int n = 0;
        while(tmp){
                fprintf(stdout,"%d ", tmp->val);
                n++;
                tmp = tmp->next;
        }
        fprintf(stdout,"\t(%d)\n",n);
}

void free_queue(queue q)
{
        node tmp = HEAD(q);
        node p = NULL;
        while(tmp){
                /* fprintf(stdout,"%d ", tmp->val); */
                p = tmp;
                tmp = tmp->next;

                free(p);
        }
        free(q);

}

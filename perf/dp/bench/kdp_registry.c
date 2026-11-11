#include <string.h>
#include <stdio.h>

#include "kdp_kernel.h"

#define KDP_MAX_KERNELS 128

void* kdp_kernel_state = NULL;

static const struct kdp_kernel* kernels[KDP_MAX_KERNELS];
static int n_kernels = 0;

void kdp_register_kernel(const struct kdp_kernel* k)
{
        if(n_kernels < KDP_MAX_KERNELS){
                kernels[n_kernels++] = k;
        }else{
                fprintf(stderr, "kdp: too many kernels, ignoring %s\n", k->name);
        }
}

const struct kdp_kernel* kdp_find_kernel(const char* name)
{
        int i;
        for(i = 0; i < n_kernels; i++){
                if(strcmp(kernels[i]->name, name) == 0){
                        return kernels[i];
                }
        }
        return NULL;
}

int kdp_n_kernels(void)
{
        return n_kernels;
}

const struct kdp_kernel* kdp_kernel_at(int i)
{
        return (i >= 0 && i < n_kernels) ? kernels[i] : NULL;
}

int kdp_resolve(const struct kdp_kernel* k, struct kdp_resolved* r)
{
        const struct kdp_kernel* ref = kdp_find_kernel("ref");
        int i;
        if(!ref){
                fprintf(stderr, "kdp: reference kernel not registered\n");
                return 1;
        }
        r->k = k;
        for(i = 0; i < 3; i++){
                r->fwd[i] = k->fwd[i] ? k->fwd[i] : ref->fwd[i];
                r->bwd[i] = k->bwd[i] ? k->bwd[i] : ref->bwd[i];
                r->meet[i] = k->meet[i] ? k->meet[i] : ref->meet[i];
        }
        return 0;
}

#ifndef RNG_H
#define RNG_H


#include "tldevel.h"
#include "stdint.h"


struct rng_state{
        uint64_t s[4];

};
extern struct rng_state* init_rng(uint64_t seed);
extern void free_rng(struct rng_state* rng);

extern double tl_random_double(struct rng_state* rng);
extern int tl_random_int(struct rng_state* rng,int a);

#endif

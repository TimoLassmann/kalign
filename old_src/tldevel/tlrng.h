#ifndef TLRNG_H
#define TLRNG_H

#include "stdint.h"

#ifndef kalign_extern
#ifdef __cplusplus
#define kalign_extern extern "C"
#else
#define kalign_extern extern
#endif
#endif


typedef struct rng_state rng_state;

kalign_extern struct rng_state* init_rng(uint64_t seed);
kalign_extern struct rng_state* init_rng_from_rng(struct rng_state* rng);

kalign_extern void free_rng(struct rng_state* rng);

kalign_extern double tl_random_double(struct rng_state* rng);
kalign_extern double tl_random_gaussian(struct rng_state* rng, double mu, double sigma);
kalign_extern double tl_random_gamma(struct rng_state* rng, double shape, double scale);
kalign_extern int tl_random_int(struct rng_state* rng,int a);


#endif

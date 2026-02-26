#ifndef SEQUENCE_DISTANCE_H
#define SEQUENCE_DISTANCE_H

#include <stdint.h>

#ifdef SEQUENCE_DISTANCE_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif





struct msa;
/* #include "alignment_parameters.h" */

EXTERN float** d_estimation(struct msa* msa, int* samples, int num_samples,int pair);
EXTERN float calc_distance(uint8_t *seq_a, uint8_t *seq_b, int len_a,
                           int len_b);

#undef SEQUENCE_DISTANCE_IMPORT
#undef EXTERN

#endif

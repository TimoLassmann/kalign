#ifndef EUCLIDIAN_DIST_H
#define EUCLIDIAN_DIST_H

#ifdef EUCLIDEAN_DIST_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

EXTERN int edist_256(const float* a,const float* b, const int len, float* ret);
EXTERN int edist_serial(const float* a,const float* b,const int len, float* ret);
EXTERN int edist_serial_d(const double* a,const double* b,const int len, double* ret);


#undef EUCLIDEAN_DIST_IMPORT
#undef EXTERN

#endif

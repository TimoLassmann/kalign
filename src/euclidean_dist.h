#ifndef EUCLIDIAN_DIST_H
#define EUCLIDIAN_DIST_H

#include "global.h"



extern  int edist_256(const float* a,const float* b, const int len, float* ret);
extern int edist_serial(const float* a,const float* b,const int len, float* ret);
extern int edist_serial_d(const double* a,const double* b,const int len, double* ret);

#endif

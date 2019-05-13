#ifndef KMEANS_H
#define KMEANS_H


#include <float.h>
#include "tldevel.h"


#include "misc.h"
#include "euclidean_dist.h"

double** kmeans(double** data, int len_a,int len_b, int k);


#endif

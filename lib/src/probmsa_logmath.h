#ifndef PROBMSA_LOGMATH_H
#define PROBMSA_LOGMATH_H

#include <math.h>
#include <float.h>

/* Log-space arithmetic for pair-HMM.
   All probabilities are stored as natural logarithms.
   -INFINITY represents zero probability. */

static inline double prob2log(double p)
{
        if(p <= 0.0){
                return -INFINITY;
        }
        return log(p);
}

static inline double log2prob(double lp)
{
        if(lp == -INFINITY){
                return 0.0;
        }
        return exp(lp);
}

/* Numerically stable log(exp(a) + exp(b)) */
static inline double logsum(double a, double b)
{
        if(a == -INFINITY){
                return b;
        }
        if(b == -INFINITY){
                return a;
        }
        if(a > b){
                return a + log1p(exp(b - a));
        }
        return b + log1p(exp(a - b));
}

#define LOGSUM3(a,b,c) logsum(logsum(a,b),c)
#define LOGSUM5(a,b,c,d,e) logsum(LOGSUM3(a,b,c), logsum(d,e))

#endif

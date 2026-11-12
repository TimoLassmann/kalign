/* The plain version of edist must round exactly like the AVX2 version, so no
   compiler may fuse a multiply and an add (GCC does so by default on Arm). */
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("fp-contract=off")
#endif
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#endif

#include "kalign_simd.h"

#if KALIGN_AVX2
#include <immintrin.h>
#endif

#include <math.h>
#include "float.h"
#include "tldevel.h"
#include "tlrng.h"
#include "esl_stopwatch.h"
#define EUCLIDEAN_DIST_IMPORT

#include "euclidean_dist.h"







#ifdef UTEST_EDIST
int main(void)
{
        struct rng_state* rng;

        float** mat = NULL;
        double r;
        float d1;
        float d2;
        int i,j,c;
        int max_iter = 10;
        int num_element = 128;
//        mat = galloc(mat,1000,8,0.0);
//void* _mm_malloc (size_t size, size_t align)

        MMALLOC(mat, sizeof(float*)* 100);
        for(i = 0; i < 100;i++){
                mat[i] = NULL;
                mat[i] = kalign_aligned_malloc(sizeof(float)*num_element);
        }

        RUNP( rng =init_rng(0));
        //srand48_r(time(NULL), &randBuffer);

        for(i =0; i < 100;i++){
                for(j = 0; j <num_element;j++){
                        //drand48_r(&randBuffer,&r);
                        r = tl_random_double(rng);
                        mat[i][j] = (float) r;
                }
        }
        LOG_MSG("Check for correctness.");
        for(i = 0; i < 100;i++){
                for(j = 0; j <= i;j++){
                        edist_serial(mat[i], mat[j], num_element, &d1);
                        edist(mat[i], mat[j], num_element, &d2);

                        if(fabsf(d1-d2) > 10e-6){
                                ERROR_MSG("DIFFER: %d\t%d\t%f\t%f  (%e %e)\n", i,j,d1,d2, fabsf(d1-d2), FLT_EPSILON);

                        }
                }
        }
        DECLARE_TIMER(t);

        LOG_MSG("Timing serial");
        START_TIMER(t);
        for(c = 0; c < max_iter;c++){
        for(i = 0; i < 100;i++){
                for(j = 0; j <= i;j++){
                        edist_serial(mat[i], mat[j], num_element, &d1);
                }
        }
        }
        STOP_TIMER(t);
        GET_TIMING(t);
        //LOG_MSG("%f\tsec.",GET_TIMING(t));

        LOG_MSG("Timing edist");
        START_TIMER(t);
        for(c = 0; c < max_iter; c++){
        for(i = 0; i < 100;i++){
                for(j = 0; j <= i;j++){

                        edist(mat[i], mat[j], num_element, &d2);

                }
        }
        }
        STOP_TIMER(t);
        GET_TIMING(t);
        //LOG_MSG("%f\tsec.",GET_TIMING(t));

        for(i = 0; i < 100;i++){
                kalign_aligned_free(mat[i]);
        }
        MFREE(mat);

        MFREE(rng);
        DESTROY_TIMER(t);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

#endif



int edist_serial(const float* a,const float* b,const int len, float* ret)
{
        int i;
        float d = 0.0f;
        float t;

        for(i = 0; i < len;i++){
                t = (a[i] - b[i]);
                d += t *t;
        }
        *ret = sqrtf(d);
        return OK;
}


int edist_serial_d(const double* a,const double* b,const int len, double* ret)
{
        int i;
        double d = 0.0f;
        double t;

        for(i = 0; i < len;i++){
                t = (a[i] - b[i]);
                d += t *t;
        }

        *ret = sqrt(d);
        return OK;
}

#if KALIGN_AVX2

/* AVX2 horizontal sum of 8 lanes: ((s0+s4)+(s1+s5))+((s2+s6)+(s3+s7)). */
KALIGN_TARGET_AVX2
static float hsum_ps_sse3(__m128 v)
{
        __m128 shuf = _mm_movehdup_ps(v);        /* broadcast elements 3,1 to 2,0 */
        __m128 sums = _mm_add_ps(v, shuf);
        shuf        = _mm_movehl_ps(shuf, sums); /* high half -> low half */
        sums        = _mm_add_ss(sums, shuf);
        return        _mm_cvtss_f32(sums);
}

KALIGN_TARGET_AVX2
static float hsum256_ps_avx(__m256 v)
{
        __m128 vlow  = _mm256_castps256_ps128(v);
        __m128 vhigh = _mm256_extractf128_ps(v, 1); /* high 128 */
        vlow  = _mm_add_ps(vlow, vhigh);     /* add the low 128 */
        return hsum_ps_sse3(vlow);
}

/* a and b must be 32-byte aligned and padded to a multiple of 8 floats with
   equal values (zeros) in the padding. */
KALIGN_TARGET_AVX2
int edist_256(const float* a,const float* b, const int len, float* ret)
{
        float d = 0.0f;
        int i;
        __m256 xmm1;
        __m256 xmm2;
        __m256 r = _mm256_set1_ps(0.0f);
        for(i = 0;i < len;i+=8){
                xmm1 = _mm256_load_ps(a);
                xmm2 = _mm256_load_ps(b);

                xmm1 =  _mm256_sub_ps(xmm1, xmm2);

                xmm1 = _mm256_mul_ps(xmm1, xmm1);

                r = _mm256_add_ps(r, xmm1);
                a+=8;
                b+=8;
        }
        d = hsum256_ps_avx(r);

        *ret = sqrtf(d);
        return OK;
}

#endif

/* Plain version of edist_256: same operations in the same order.  Eight lane
   accumulators over blocks of 8 (the padding of the AVX2 version adds exactly
   0 to a lane, so a partial last block just skips the missing lanes), then
   the same horizontal sum.  Reads only a[0..len-1], b[0..len-1]. */
int edist_plain(const float* a,const float* b, const int len, float* ret)
{
        float s[8];
        float t;
        float lo[4];
        float h0,h1;
        int i,k;
        int full = len & ~7;

        for(k = 0; k < 8;k++){
                s[k] = 0.0f;
        }
        for(i = 0; i < full;i+=8){
                for(k = 0; k < 8;k++){
                        t = a[i+k] - b[i+k];
                        t = t * t;
                        s[k] = s[k] + t;
                }
        }
        for(k = 0; k < len - full;k++){
                t = a[full+k] - b[full+k];
                t = t * t;
                s[k] = s[k] + t;
        }
        for(k = 0; k < 4;k++){
                lo[k] = s[k] + s[k+4];
        }
        h0 = lo[0] + lo[1];
        h1 = lo[2] + lo[3];
        *ret = sqrtf(h0 + h1);
        return OK;
}

/* Euclidean distance of two float vectors, identical result on every CPU.
   With AVX2 the arrays must be aligned (kalign_aligned_malloc) and padded
   to a multiple of 8 floats. */
int edist(const float* a,const float* b, const int len, float* ret)
{
#if KALIGN_AVX2
        if(kalign_simd_on()){
                return edist_256(a,b,len,ret);
        }
#endif
        return edist_plain(a,b,len,ret);
}

/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019 Timo Lassmann

    This file is part of kalign.

    Kalign is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#include "euclidean_dist.h"
#include "rng.h"
#include <xmmintrin.h>
#include <immintrin.h>
#include "float.h"
/* These functions were taken from:  */
/* https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86 */
float hsum256_ps_avx(__m256 v);
float hsum_ps_sse3(__m128 v);


#ifdef ITEST
int main(int argc, char *argv[])
{
        struct rng_state* rng;

        float** mat = NULL;
        double r;
        float d1,d2;
        int i,j,c;
        int max_iter = 10;
        int num_element = 128;
//        mat = galloc(mat,1000,8,0.0);
//void* _mm_malloc (size_t size, size_t align)

        MMALLOC(mat, sizeof(float*)* 1000);
        for(i = 0; i < 1000;i++){
                mat[i] = _mm_malloc(sizeof(float)*num_element, 32);
        }

        RUNP( rng =init_rng(0));
        //srand48_r(time(NULL), &randBuffer);

        for(i =0; i < 1000;i++){
                for(j = 0; j <num_element;j++){
                        //drand48_r(&randBuffer,&r);
                        r = tl_random_double(rng);
                        mat[i][j] = (float) r;
                }
        }
        LOG_MSG("Check for correctness.");
#ifdef HAVE_AVX2
        for(i = 0; i < 1000;i++){
                for(j = 0; j <= i;j++){
                        edist_serial(mat[i], mat[j], num_element, &d1);
                        edist_256(mat[i], mat[j], num_element, &d2);
                        if(fabsf(d1-d2) > 10e-6){
                                ERROR_MSG("DIFFER: %d\t%d\t%f\t%f  (%e %e)\n", i,j,d1,d2, fabsf(d1-d2), FLT_EPSILON);

                        }
                }
        }
#endif
        DECLARE_TIMER(t);

        LOG_MSG("Timing serial");
        START_TIMER(t);
        for(c = 0; c < max_iter;c++){
        for(i = 0; i < 1000;i++){
                for(j = 0; j <= i;j++){
                        edist_serial(mat[i], mat[j], num_element, &d1);


                }
        }
        }
        STOP_TIMER(t);
        LOG_MSG("%f\tsec.",GET_TIMING(t));

#ifdef HAVE_AVX2
        LOG_MSG("Timing AVX");
        START_TIMER(t);
        for(c = 0; c < max_iter; c++){
        for(i = 0; i < 1000;i++){
                for(j = 0; j <= i;j++){

                        edist_256(mat[i], mat[j], num_element, &d2);

                }
        }
        }
        STOP_TIMER(t);
        LOG_MSG("%f\tsec.",GET_TIMING(t));

#endif
        for(i = 0; i < 1000;i++){
                _mm_free(mat[i]);
        }
        MFREE(mat);

        MFREE(rng);
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

#ifdef HAVE_AVX2

int edist_256(const float* a,const float* b, const int len, float* ret)
{

        float d = 0.0f;
        register int i;
        __m256 xmm1;// = _mm256_load_ps(a);
        __m256 xmm2;// = _mm256_load_ps(b);
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




float hsum256_ps_avx(__m256 v)
{
        __m128 vlow  = _mm256_castps256_ps128(v);
        __m128 vhigh = _mm256_extractf128_ps(v, 1); // high 128
        vlow  = _mm_add_ps(vlow, vhigh);     // add the low 128
        return hsum_ps_sse3(vlow);         // and inline the sse3 version, which is optimal for AVX
        // (no wasted instructions, and all of them are the 4B minimum)
}

float hsum_ps_sse3(__m128 v)
{
        __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
        __m128 sums = _mm_add_ps(v, shuf);
        shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
        sums        = _mm_add_ss(sums, shuf);
        return        _mm_cvtss_f32(sums);
}

#endif

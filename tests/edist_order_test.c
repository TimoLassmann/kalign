/* edist_order_test.c -- the plain edist must round exactly like edist_256.
 *
 * Compares edist_plain() and edist() bitwise (memcmp of the float) with an
 * independent emulation of the AVX2 order (zero-padded to a multiple of 8,
 * eight accumulators, ((s0+s4)+(s1+s5))+((s2+s6)+(s3+s7)), sqrtf) for lengths
 * 0..200 and random data.  When built on x86-64 with AVX2 code compiled in
 * (KALIGN_AVX2) and run on a CPU with AVX2, it also compares edist_256().
 *
 * Link: tests/edist_order_test.c lib/src/euclidean_dist.c lib/src/kalign_simd.c
 *       lib/src/tldevel.c lib/src/tlrng.c lib/src/esl_stopwatch.c -lm -lpthread
 *       (include path lib/src; tldevel needs the usual project defines)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "kalign_simd.h"
#include "euclidean_dist.h"

/* The reference sum below must not be fused into multiply-adds.  Clang obeys
   the pragma; GCC ignores it, so tests/CMakeLists.txt also compiles this file
   with -ffp-contract=off. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#endif

static uint64_t st = 88172645463325252ULL;
static uint32_t rnd(void){ st ^= st << 13; st ^= st >> 7; st ^= st << 17; return (uint32_t)(st >> 11); }
static float rf(void){ return ((float)rnd() / 2097152.0f) * 3.0f - 1.5f; }

/* emulation of edist_256 on padded input */
static float emu(const float* a, const float* b, int len)
{
        float s0=0,s1=0,s2=0,s3=0,s4=0,s5=0,s6=0,s7=0;
        int i;
        for(i = 0; i < len; i += 8){
                float t0 = a[i]-b[i];     t0 = t0*t0; s0 = s0 + t0;
                float t1 = a[i+1]-b[i+1]; t1 = t1*t1; s1 = s1 + t1;
                float t2 = a[i+2]-b[i+2]; t2 = t2*t2; s2 = s2 + t2;
                float t3 = a[i+3]-b[i+3]; t3 = t3*t3; s3 = s3 + t3;
                float t4 = a[i+4]-b[i+4]; t4 = t4*t4; s4 = s4 + t4;
                float t5 = a[i+5]-b[i+5]; t5 = t5*t5; s5 = s5 + t5;
                float t6 = a[i+6]-b[i+6]; t6 = t6*t6; s6 = s6 + t6;
                float t7 = a[i+7]-b[i+7]; t7 = t7*t7; s7 = s7 + t7;
        }
        float l0 = s0 + s4, l1 = s1 + s5, l2 = s2 + s6, l3 = s3 + s7;
        float h0 = l0 + l1, h1 = l2 + l3;
        float d = h0 + h1;
        return sqrtf(d);
}

int main(void)
{
        int len, rep, bad = 0, n = 0, nsimd = 0;
        float* a = kalign_aligned_malloc(sizeof(float) * 208);
        float* b = kalign_aligned_malloc(sizeof(float) * 208);
        if(!a || !b) return 2;
        for(len = 0; len <= 200; len++){
                for(rep = 0; rep < 200; rep++){
                        int i;
                        float e, p, d;
                        for(i = 0; i < 208; i++){
                                a[i] = (i < len) ? rf() : 0.0f;
                                b[i] = (i < len) ? rf() : 0.0f;
                        }
                        e = emu(a, b, len);
                        edist_plain(a, b, len, &p);
                        edist(a, b, len, &d);
                        n++;
                        if(memcmp(&e, &p, sizeof(float)) || memcmp(&e, &d, sizeof(float))){
                                if(bad++ < 10) printf("MISMATCH plain/edist len %d rep %d %.9g %.9g %.9g\n", len, rep, e, p, d);
                        }
#if KALIGN_AVX2
                        if(kalign_simd_on()){
                                float v;
                                edist_256(a, b, len, &v);
                                nsimd++;
                                if(memcmp(&e, &v, sizeof(float))){
                                        if(bad++ < 10) printf("MISMATCH avx2 len %d rep %d %.9g %.9g\n", len, rep, e, v);
                                }
                        }
#endif
                }
        }
        kalign_aligned_free(a);
        kalign_aligned_free(b);
        printf("edist_order_test: %d cases, %d vs AVX2, %d mismatches (simd=%s)\n", n, nsimd, bad, kalign_simd_name());
        return bad ? 1 : 0;
}

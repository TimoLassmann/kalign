/* Run-time vector switch and aligned allocator (lib/src/kalign_simd.h). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "kalign_simd.h"

static int fails = 0;

#define CHECK(cond, ...) do{ if(!(cond)){ fprintf(stderr, "FAIL line %d: ", __LINE__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); fails++; } }while(0)

int main(void)
{
        size_t sizes[] = {0, 1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 255, 256, 1000, 4095, 4096, 4097, 65536, 1000003};
        size_t i, j;
        int on0;
        const char* name0;

        on0 = kalign_simd_on();
        name0 = kalign_simd_name();
        printf("kalign_simd_name: %s (kalign_simd_on = %d)\n", name0, on0);
        CHECK(strcmp(name0, "avx2") == 0 || strcmp(name0, "neon") == 0 || strcmp(name0, "scalar") == 0, "unknown name '%s'", name0);
        CHECK((on0 != 0) == (strcmp(name0, "scalar") != 0), "name '%s' does not match kalign_simd_on() = %d", name0, on0);
#if !KALIGN_AVX2 && !KALIGN_NEON
        CHECK(on0 == 0, "no vector kernels are compiled in but kalign_simd_on() = %d", on0);
#endif
        {
                const char* e = getenv("KALIGN_SIMD");
                if(e && strcmp(e, "scalar") == 0){
                        CHECK(on0 == 0, "KALIGN_SIMD=scalar but kalign_simd_on() = %d", on0);
                }
        }

        /* aligned allocator */
        for(i = 0; i < sizeof(sizes) / sizeof(sizes[0]); i++){
                unsigned char* p = kalign_aligned_malloc(sizes[i]);
                if(p == NULL){
                        /* size 0 may return NULL or a valid pointer; others must succeed */
                        CHECK(sizes[i] == 0, "kalign_aligned_malloc(%zu) returned NULL", sizes[i]);
                        continue;
                }
                CHECK(((uintptr_t)p % KALIGN_SIMD_ALIGN) == 0, "kalign_aligned_malloc(%zu) = %p is not %d-byte aligned", sizes[i], (void*)p, KALIGN_SIMD_ALIGN);
                for(j = 0; j < sizes[i]; j++){
                        p[j] = (unsigned char)(j * 31u + 7u);
                }
                for(j = 0; j < sizes[i]; j++){
                        if(p[j] != (unsigned char)(j * 31u + 7u)){
                                CHECK(0, "memory of size %zu not writable at %zu", sizes[i], j);
                                break;
                        }
                }
                kalign_aligned_free(p);
        }
        kalign_aligned_free(NULL);

        /* many live blocks at once */
        {
                void* q[64];
                for(i = 0; i < 64; i++){
                        q[i] = kalign_aligned_malloc(1 + i * 37);
                        CHECK(q[i] != NULL && ((uintptr_t)q[i] % KALIGN_SIMD_ALIGN) == 0, "block %zu bad", i);
                }
                for(i = 0; i < 64; i++){
                        kalign_aligned_free(q[i]);
                }
        }

        /* force scalar */
        kalign_simd_force_scalar(1);
        CHECK(kalign_simd_on() == 0, "kalign_simd_on() = %d after force_scalar(1)", kalign_simd_on());
        CHECK(strcmp(kalign_simd_name(), "scalar") == 0, "name '%s' after force_scalar(1)", kalign_simd_name());
        kalign_simd_force_scalar(0);
        if(getenv("KALIGN_SIMD") == NULL){
                CHECK(kalign_simd_on() == on0, "kalign_simd_on() = %d after force_scalar(0), was %d", kalign_simd_on(), on0);
                CHECK(strcmp(kalign_simd_name(), name0) == 0, "name '%s' after force_scalar(0), was '%s'", kalign_simd_name(), name0);
        }

        if(fails){
                fprintf(stderr, "%d check(s) failed\n", fails);
                return 1;
        }
        printf("kalign_simd_test: all checks passed\n");
        return 0;
}

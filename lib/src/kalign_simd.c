/* kalign_simd.c -- see kalign_simd.h */
#include "kalign_simd.h"

#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

#if KALIGN_AVX2
#include <cpuid.h>

/* AVX2 is usable when the CPU has it (CPUID leaf 7) and the operating system
   saves the 256-bit registers on a task switch (OSXSAVE + XGETBV).  POPCNT is
   required as well: GCC switches it on inside every target("avx2") function,
   so a bit count in an AVX2 kernel becomes a popcnt instruction.  Every CPU
   with AVX2 has POPCNT; only a virtual machine with a masked CPUID can
   report one without the other.  The
   compiler's __builtin_cpu_supports is not used: it depends on run-time
   support code that is missing or not initialised in some tool chains
   (static musl builds report 0 on a CPU with AVX2). */
static int cpu_has_avx2(void)
{
        unsigned int a = 0, b = 0, c = 0, d = 0;
        unsigned int lo = 0, hi = 0;

        if(!__get_cpuid(1, &a, &b, &c, &d)){
                return 0;
        }
        if(!(c & (1u << 27)) || !(c & (1u << 28))){     /* OSXSAVE, AVX */
                return 0;
        }
        if(!(c & (1u << 23))){                          /* POPCNT */
                return 0;
        }
        /* xgetbv with ecx = 0 (as bytes: old assemblers do not know the name) */
        __asm__ volatile(".byte 0x0f, 0x01, 0xd0" : "=a"(lo), "=d"(hi) : "c"(0));
        (void)hi;
        if((lo & 6u) != 6u){                            /* XMM and YMM state saved */
                return 0;
        }
        if(__get_cpuid_max(0, NULL) < 7){
                return 0;
        }
        __cpuid_count(7, 0, a, b, c, d);
        return (b >> 5) & 1u;                           /* AVX2 */
}
#endif

static _Atomic int simd_state = -1;     /* -1 not asked yet, 0 plain, 1 vector; every writer stores the same value */
static _Atomic int simd_forced_scalar = 0;

static int simd_detect(void)
{
        const char* e = getenv("KALIGN_SIMD");
        if(e && strcmp(e, "scalar") == 0){
                return 0;
        }
#if KALIGN_AVX2
        return cpu_has_avx2();
#elif KALIGN_NEON
        return 1;
#else
        return 0;
#endif
}

int kalign_simd_on(void)
{
        int s;
        if(atomic_load_explicit(&simd_forced_scalar, memory_order_relaxed)){
                return 0;
        }
        s = atomic_load_explicit(&simd_state, memory_order_relaxed);
        if(s < 0){
                s = simd_detect();
                atomic_store_explicit(&simd_state, s, memory_order_relaxed);
        }
        return s;
}

const char* kalign_simd_name(void)
{
        if(!kalign_simd_on()){
                return "scalar";
        }
#if KALIGN_AVX2
        return "avx2";
#elif KALIGN_NEON
        return "neon";
#else
        return "scalar";
#endif
}

void kalign_simd_force_scalar(int on)
{
        atomic_store_explicit(&simd_forced_scalar, on ? 1 : 0, memory_order_relaxed);
}

#ifdef KALIGN_DEBUG_ALIGNED
/* Checked variant (ASAN builds): the block gets a header of one alignment unit
   and the caller receives the address after it.  Releasing such a block with
   free() or MFREE then frees an address that malloc never returned, and
   kalign_aligned_free on memory from malloc or MMALLOC reads before the
   block; the address sanitizer reports both.  The header holds a mark that
   kalign_aligned_free checks, so a mismatch also stops without a sanitizer. */
#include <stdint.h>
#include <stdio.h>
#define ALIGNED_MARK 0x4b414c49474e3332ULL      /* "KALIGN32" */

void* kalign_aligned_malloc(size_t size)
{
        void* p = NULL;
        uint64_t mark = ALIGNED_MARK;
        if(size == 0){
                size = KALIGN_SIMD_ALIGN;
        }
        if(size > SIZE_MAX - KALIGN_SIMD_ALIGN){
                return NULL;
        }
        if(posix_memalign(&p, KALIGN_SIMD_ALIGN, size + KALIGN_SIMD_ALIGN) != 0){
                return NULL;
        }
        memcpy(p, &mark, sizeof(mark));
        return (char*)p + KALIGN_SIMD_ALIGN;
}

void kalign_aligned_free(void* p)
{
        uint64_t mark = 0;
        char* base;
        if(p == NULL){
                return;
        }
        base = (char*)p - KALIGN_SIMD_ALIGN;
        memcpy(&mark, base, sizeof(mark));
        if(mark != ALIGNED_MARK){
                fprintf(stderr, "kalign_aligned_free: %p does not come from kalign_aligned_malloc\n", p);
                abort();
        }
        mark = 0;
        memcpy(base, &mark, sizeof(mark));      /* a second free of the same block fails the check */
        free(base);
}
#else
void* kalign_aligned_malloc(size_t size)
{
        void* p = NULL;
        if(size == 0){
                size = KALIGN_SIMD_ALIGN;
        }
        if(posix_memalign(&p, KALIGN_SIMD_ALIGN, size) != 0){
                return NULL;
        }
        return p;
}

void kalign_aligned_free(void* p)
{
        free(p);        /* posix_memalign memory is released with free() */
}
#endif

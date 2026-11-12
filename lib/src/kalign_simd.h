/* kalign_simd.h -- the one place that decides which vector code runs, and that
 * hands out the aligned memory the vector code needs.
 *
 * Kalign is compiled for the plain instruction set of the platform (x86-64:
 * the original 64-bit instruction set, no -mavx2 anywhere).  The x86 vector
 * kernels are separate functions marked KALIGN_TARGET_AVX2; they are called
 * only when kalign_simd_on() says the CPU and the operating system support
 * AVX2.  On 64-bit Arm every CPU has NEON, so the NEON kernels need no check,
 * but they go through the same switch so that the plain kernels can be tested
 * on every machine.
 *
 * Rules for code that uses this header:
 *  - A function that contains AVX2 intrinsics is marked KALIGN_TARGET_AVX2 (so
 *    are the always-inline helpers it calls) and is reached only through
 *    "if(kalign_simd_on())".  Nothing else may contain AVX2 code: the compiler
 *    rejects an intrinsic in an unmarked function, because no file is built
 *    with -mavx2.
 *  - The marker is AVX2 only.  Never add FMA or any other feature: the vector
 *    kernels must round exactly like the plain kernels (the alignments are
 *    byte-identical on every CPU), and a fused multiply-add rounds differently.
 *  - Call kalign_simd_on() once per kernel call (a whole DP pass, a whole
 *    distance), never per element.
 *  - Memory that a vector kernel may load with aligned loads comes from
 *    kalign_aligned_malloc() and goes back through kalign_aligned_free(),
 *    whichever kernel runs.  Do not mix with _mm_malloc/_mm_free, free or
 *    MFREE.  In release builds a mismatch goes unnoticed (both end in
 *    free()); builds with -DKALIGN_DEBUG_ALIGNED (the ASAN build type) give
 *    each aligned block a header, so a mismatch is reported.
 *  - On x86-64 the library is compiled with -ffp-contract=off (CMake and
 *    build.zig): a build with FMA switched on for the whole program
 *    (-march=x86-64-v3, -march=native) must not fuse a multiply and an add
 *    in the plain kernels, or they would round differently from the AVX2
 *    kernels and from a build without FMA.
 */
#ifndef KALIGN_SIMD_H
#define KALIGN_SIMD_H

#include <stddef.h>

/* KALIGN_AVX2: the AVX2 kernels are compiled in (x86-64, GCC or Clang, not
   switched off with -DKALIGN_NO_AVX2 = CMake -DENABLE_AVX2=OFF). */
#if defined(__x86_64__) && !defined(KALIGN_NO_AVX2) && (defined(__GNUC__) || defined(__clang__))
#define KALIGN_AVX2 1
#define KALIGN_TARGET_AVX2 __attribute__((target("avx2")))
#else
#define KALIGN_AVX2 0
#define KALIGN_TARGET_AVX2
#endif

/* Plain C code that the compiler vectorises itself (no intrinsics) is made
   faster on AVX2 CPUs by compiling it twice: the body goes in a
   KALIGN_ALWAYS_INLINE function, and two thin wrappers call it, one plain and
   one marked KALIGN_TARGET_AVX2; the caller picks one with kalign_simd_on().
   The AVX2 copy is the same C code, so for integer code (and float code,
   which the compiler does not reorder: no -ffast-math, no FMA) it gives the
   same results. */
#if defined(__GNUC__) || defined(__clang__)
#define KALIGN_ALWAYS_INLINE static inline __attribute__((always_inline))
#else
#define KALIGN_ALWAYS_INLINE static inline
#endif

/* KALIGN_NEON: the NEON kernels are compiled in. */
#if defined(__ARM_NEON)
#define KALIGN_NEON 1
#else
#define KALIGN_NEON 0
#endif

/* 1: use the vector kernels.  0: use the plain kernels (no usable vector
   instructions, or KALIGN_SIMD=scalar in the environment, or forced by
   kalign_simd_force_scalar).  Cheap and thread-safe; the CPU is asked once. */
int kalign_simd_on(void);

/* "avx2", "neon" or "scalar": what kalign_simd_on() selects now. */
const char* kalign_simd_name(void);

/* For tests: on != 0 selects the plain kernels until called again with 0.
   Call it only while no alignment is running. */
void kalign_simd_force_scalar(int on);

/* Aligned memory for arrays the vector kernels read (32 bytes: one AVX2
   register).  kalign_aligned_malloc returns NULL on failure.  The memory is
   released with kalign_aligned_free (NULL is allowed). */
#define KALIGN_SIMD_ALIGN 32
void* kalign_aligned_malloc(size_t size);
void kalign_aligned_free(void* p);

#endif

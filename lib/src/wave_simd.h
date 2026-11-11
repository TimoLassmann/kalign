/* wave_simd.h -- 8-lane float vectors for the wavefront DP passes (aln_wave.c).
 *
 * One type, vf8, and the few operations the kernels need, implemented with
 * NEON (two float32x4_t) or AVX2 (one __m256).  WAVE_SIMD is 1 when one of
 * them is available, 0 otherwise (callers then use the scalar kernels).
 *
 * Exactness rules (the kernels must reproduce the scalar DP bit for bit):
 *  - vf8_max(a, b) returns b when a > b is false, like MAX(a, b) in the DP
 *    code ((a) > (b) ? (a) : (b)); this is the semantics of AVX max_ps.  NEON
 *    fmax differs only for +0/-0 ties, which the DP never produces.
 *  - vf8_muladd(acc, w, x) must round like the scalar statement
 *    acc += w * x as the compiler builds it for the reference kernels: fused
 *    on AArch64 and on x86 with FMA enabled (compilers contract it there),
 *    a separate multiply and add otherwise (x86 -mavx2 without -mfma, the
 *    default Kalign build).  bench_dp --check verifies this per platform.
 */
#ifndef WAVE_SIMD_H
#define WAVE_SIMD_H

/* Force inlining so that unoptimised (-O0, ASAN, Debug) builds do not turn
   every vector operation into a call: without it the wave kernels are a few
   hundred times slower than the scalar code in such builds. */
#if defined(__GNUC__) || defined(__clang__)
#define VF8_INLINE static inline __attribute__((always_inline))
#else
#define VF8_INLINE static inline
#endif

#if defined(_MSC_VER) && !defined(__clang__)
/* MSVC: the kernels use GNU attributes and C11 thread-local storage */
#define WAVE_SIMD 0

#elif defined(__ARM_NEON)
#include <arm_neon.h>
#define WAVE_SIMD 1

typedef struct { float32x4_t lo, hi; } vf8;

VF8_INLINE vf8 vf8_set1(float x) { vf8 r = { vdupq_n_f32(x), vdupq_n_f32(x) }; return r; }
VF8_INLINE vf8 vf8_loadu(const float* p) { vf8 r = { vld1q_f32(p), vld1q_f32(p + 4) }; return r; }
VF8_INLINE void vf8_storeu(float* p, vf8 v) { vst1q_f32(p, v.lo); vst1q_f32(p + 4, v.hi); }
VF8_INLINE vf8 vf8_add(vf8 a, vf8 b) { vf8 r = { vaddq_f32(a.lo, b.lo), vaddq_f32(a.hi, b.hi) }; return r; }
VF8_INLINE vf8 vf8_sub(vf8 a, vf8 b) { vf8 r = { vsubq_f32(a.lo, b.lo), vsubq_f32(a.hi, b.hi) }; return r; }
VF8_INLINE vf8 vf8_max(vf8 a, vf8 b) { vf8 r = { vmaxq_f32(a.lo, b.lo), vmaxq_f32(a.hi, b.hi) }; return r; }
VF8_INLINE vf8 vf8_div(vf8 a, vf8 b) { vf8 r = { vdivq_f32(a.lo, b.lo), vdivq_f32(a.hi, b.hi) }; return r; }
VF8_INLINE vf8 vf8_muladd(vf8 acc, vf8 w, vf8 x)
{
        vf8 r = { vfmaq_f32(acc.lo, w.lo, x.lo), vfmaq_f32(acc.hi, w.hi, x.hi) };
        return r;
}
/* [x, v0, v1, ..., v6] */
VF8_INLINE vf8 vf8_shift_in(vf8 v, float x)
{
        vf8 r = { vextq_f32(vdupq_n_f32(x), v.lo, 3), vextq_f32(v.lo, v.hi, 3) };
        return r;
}
VF8_INLINE float vf8_last(vf8 v) { return vgetq_lane_f32(v.hi, 3); }
VF8_INLINE vf8 vf8_set_lane0(vf8 v, float x) { v.lo = vsetq_lane_f32(x, v.lo, 0); return v; }
/* lanes i: p[i][off] */
VF8_INLINE vf8 vf8_gather(const float* const* p, int off)
{
        vf8 r;
        r.lo = vdupq_n_f32(0.0F);
        r.hi = vdupq_n_f32(0.0F);
        r.lo = vld1q_lane_f32(p[0] + off, r.lo, 0);
        r.lo = vld1q_lane_f32(p[1] + off, r.lo, 1);
        r.lo = vld1q_lane_f32(p[2] + off, r.lo, 2);
        r.lo = vld1q_lane_f32(p[3] + off, r.lo, 3);
        r.hi = vld1q_lane_f32(p[4] + off, r.hi, 0);
        r.hi = vld1q_lane_f32(p[5] + off, r.hi, 1);
        r.hi = vld1q_lane_f32(p[6] + off, r.hi, 2);
        r.hi = vld1q_lane_f32(p[7] + off, r.hi, 3);
        return r;
}

#elif defined(__AVX2__)
#include <immintrin.h>
#define WAVE_SIMD 1

typedef __m256 vf8;

VF8_INLINE vf8 vf8_set1(float x) { return _mm256_set1_ps(x); }
VF8_INLINE vf8 vf8_loadu(const float* p) { return _mm256_loadu_ps(p); }
VF8_INLINE void vf8_storeu(float* p, vf8 v) { _mm256_storeu_ps(p, v); }
VF8_INLINE vf8 vf8_add(vf8 a, vf8 b) { return _mm256_add_ps(a, b); }
VF8_INLINE vf8 vf8_sub(vf8 a, vf8 b) { return _mm256_sub_ps(a, b); }
VF8_INLINE vf8 vf8_max(vf8 a, vf8 b) { return _mm256_max_ps(a, b); }
VF8_INLINE vf8 vf8_div(vf8 a, vf8 b) { return _mm256_div_ps(a, b); }
VF8_INLINE vf8 vf8_muladd(vf8 acc, vf8 w, vf8 x)
{
#if defined(__FMA__)
        return _mm256_fmadd_ps(w, x, acc);
#else
        return _mm256_add_ps(acc, _mm256_mul_ps(w, x));
#endif
}
VF8_INLINE vf8 vf8_shift_in(vf8 v, float x)
{
        const __m256i rot = _mm256_setr_epi32(7, 0, 1, 2, 3, 4, 5, 6);
        return _mm256_blend_ps(_mm256_permutevar8x32_ps(v, rot), _mm256_set1_ps(x), 1);
}
VF8_INLINE float vf8_last(vf8 v)
{
        return _mm_cvtss_f32(_mm_shuffle_ps(_mm256_extractf128_ps(v, 1), _mm256_extractf128_ps(v, 1), 3));
}
VF8_INLINE vf8 vf8_set_lane0(vf8 v, float x) { return _mm256_blend_ps(v, _mm256_set1_ps(x), 1); }
VF8_INLINE vf8 vf8_gather(const float* const* p, int off)
{
        return _mm256_setr_ps(p[0][off], p[1][off], p[2][off], p[3][off],
                              p[4][off], p[5][off], p[6][off], p[7][off]);
}

#else
#define WAVE_SIMD 0
#endif

#if WAVE_SIMD
/* lane i of v (via memory; used off the hot path) */
VF8_INLINE float vf8_lane(vf8 v, int i)
{
        float t[8];
        vf8_storeu(t, v);
        return t[i];
}
/* v with lane i replaced by x */
VF8_INLINE vf8 vf8_with_lane(vf8 v, int i, float x)
{
        float t[8];
        vf8_storeu(t, v);
        t[i] = x;
        return vf8_loadu(t);
}
#endif

#endif

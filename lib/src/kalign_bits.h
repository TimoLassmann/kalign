/* kalign_bits.h -- bit counting that works on every CPU.
 *
 * On x86-64 Kalign is compiled for the original instruction set, which has
 * no bit-count instruction.  GCC then turns __builtin_popcount into a call
 * to a library routine (__popcountdi2), which is slow in the consensus code.
 * There the count is done inline with shifts, adds and one multiply (the
 * standard "SWAR" method).  Where the build allows the CPU's own instruction
 * (POPCNT on x86, CNT on 64-bit Arm), the built-in is used.  Both give the
 * same count.
 */
#ifndef KALIGN_BITS_H
#define KALIGN_BITS_H

#include <stdint.h>

#if (defined(__GNUC__) || defined(__clang__)) && (!defined(__x86_64__) || defined(__POPCNT__))
#define KALIGN_POPCOUNT_BUILTIN 1
#else
#define KALIGN_POPCOUNT_BUILTIN 0
#endif

static inline int kalign_popcount64(uint64_t x)
{
#if KALIGN_POPCOUNT_BUILTIN
        return __builtin_popcountll(x);
#else
        x = x - ((x >> 1) & 0x5555555555555555ULL);
        x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
        x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
        return (int)((x * 0x0101010101010101ULL) >> 56);
#endif
}

static inline int kalign_popcount32(uint32_t x)
{
#if KALIGN_POPCOUNT_BUILTIN
        return __builtin_popcount(x);
#else
        x = x - ((x >> 1) & 0x55555555U);
        x = (x & 0x33333333U) + ((x >> 2) & 0x33333333U);
        x = (x + (x >> 4)) & 0x0f0f0f0fU;
        return (int)((x * 0x01010101U) >> 24);
#endif
}


#endif

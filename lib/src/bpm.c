#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "tldevel.h"
#define BPM_IMPORT
#include "bpm.h"
#include  <stdalign.h>
#include <string.h>
#include "kalign_simd.h"

#define SIGMA 13
#define DIV_CEIL(a,b) (a == 0 ? 1 : a/b+(a%b == 0 ? 0 : 1))


#if KALIGN_AVX2
#include <immintrin.h>

static KALIGN_TARGET_AVX2 void bitShiftLeft256ymm (__m256i *data, int count);

/* taken from Alexander Yee: http://www.numberworld.org/y-cruncher/internals/addition.html#ks_add */
static KALIGN_TARGET_AVX2 __m256i add256(uint32_t carry, __m256i A, __m256i B);
#endif


uint8_t dyn_256(const uint8_t* t,const uint8_t* p,int n,int m)
{
        uint8_t* prev = NULL;
        uint8_t* cur = NULL;

        uint8_t* tmp = NULL;
        int i,j,c;
        if(m > 255){
                m = 255;
        }


        MMALLOC(prev, sizeof(uint8_t)* 257);
        MMALLOC(cur, sizeof(uint8_t)* 257);
        cur[0] = 0;
        //fprintf(stdout,"%d ", cur[0]);
        for(j = 1; j <= m;j++){

                cur[j] = cur[j-1] +1;
                //fprintf(stdout,"%d ", cur[j]);
        }
        //fprintf(stdout,"\n");
        tmp  = cur;
        cur = prev;
        prev = tmp;

        for(i = 1; i <= n;i++){
                cur[0] = prev[0] ;
                //fprintf(stdout,"%d ", cur[0]);
                for(j = 1; j < m;j++){
                        c = 1;
                        if(t[i-1] == p[j-1]){
                                c = 0;

                        }
                        cur[j] = prev[j-1] +c ;
                        cur[j] = MACRO_MIN(cur[j], prev[j]+1);
                        cur[j] = MACRO_MIN(cur[j], cur[j-1]+1);
                }
                j = m;

                c = 1;
                if(t[i-1] == p[j-1]){
                        c = 0;

                }
                cur[j] = prev[j-1] +c ;

                cur[j] = MACRO_MIN(cur[j], prev[j]);
                cur[j] = MACRO_MIN(cur[j], cur[j-1]+1);
                //fprintf(stdout,"\n");
                tmp  = cur;
                cur = prev;
                prev = tmp;

        }
        c = prev[m];
        MFREE(prev);
        MFREE(cur);
        return c;
ERROR:
        return 255;

}

uint8_t bpm(const uint8_t* t,const uint8_t* p,int n,int m)
{
        register uint64_t VP,VN,D0,HN,HP,X;
        register uint64_t i;//,c;
        uint64_t MASK = 0;
        int64_t diff;
        uint64_t B[13];
        int8_t k;

        if(m > 63){
                m = 63;
        }
        diff = m;
        k = m;
        for(i = 0; i < 13;i++){
                B[i] = 0;
        }

        for(i = 0; i < (uint64_t)m;i++){
                B[p[i]] |= (1ul << i);
        }

        VP = (1ul << (m))-1 ;

        VN = 0ul;

        m--;
        MASK = 1ul << (m);

        //fprintf(stdout,"BEGINNING\t%lu %lu\n",VN,VP);
        for(i = 0; i < (uint64_t) n;i++){
                //        fprintf(stdout,"%lu:\t",i);
                X = (B[t[i]] | VN);
                //fprintf(stdout,"%lu ", X);
                D0 = ((VP+(X&VP)) ^ VP) | X ;
                //fprintf(stdout,"%lu ", D0);
                HN = VP & D0;
                //fprintf(stdout,"%lu ", HN);
                HP = VN | (~(VP | D0));
                //fprintf(stdout,"%lu ", HP);
                X = HP << 1ul;
                //fprintf(stdout,"%lu ", X);
                VN = X & D0;
                //fprintf(stdout,"%lu ", VN);

                VP = (HN << 1ul) | (~(X | D0));
                //fprintf(stdout,"%lu ", VP);


                diff += (HP & MASK)? 1 : 0;// >> m;
                diff -= (HN & MASK)? 1 : 0;// >> m;
                if(diff < k){
                        k = diff;

                }
                //fprintf(stdout,"%lu ", diff);
                //fprintf(stdout,"\n");

        }
        return k;
}


#if KALIGN_AVX2
static KALIGN_TARGET_AVX2 uint8_t bpm_256_avx2(const uint8_t* t,const uint8_t* p,int n,int m)
{
        __m256i VP,VN,D0,HN,HP,X,NOTONE;
        __m256i xmm1,xmm2;
        __m256i MASK;
        __m256i B[13];

        int i,j, k,diff;

        alignas(32)  uint32_t f[13][8];
        //int ALIGNED_(64) f[8];
        if(m > 255){
                m = 255;
        }

        for(i = 0; i < 13;i++){
                for(j = 0;j < 8;j++){
                        f[i][j] =0u;
                }
        }

        for(i = 0; i < m;i++){
                f[p[i]][i/32] |= (1 << (i % 32));
        }

        /* for(i = 0; i < 13;i++){ */
        /*         B[i] = _mm256_load_si256((__m256i const*) &f[i]); */
        /* } */

        B[0] = _mm256_load_si256((__m256i const*) &f[0]);
        B[1] = _mm256_load_si256((__m256i const*) &f[1]);
        B[2] = _mm256_load_si256((__m256i const*) &f[2]);
        B[3] = _mm256_load_si256((__m256i const*) &f[3]);
        B[4] = _mm256_load_si256((__m256i const*) &f[4]);
        B[5] = _mm256_load_si256((__m256i const*) &f[5]);
        B[6] = _mm256_load_si256((__m256i const*) &f[6]);
        B[7] = _mm256_load_si256((__m256i const*) &f[7]);
        B[8] = _mm256_load_si256((__m256i const*) &f[8]);
        B[9] = _mm256_load_si256((__m256i const*) &f[9]);
        B[10] = _mm256_load_si256((__m256i const*) &f[10]);
        B[11] = _mm256_load_si256((__m256i const*) &f[11]);
        B[12] = _mm256_load_si256((__m256i const*) &f[12]);

        diff = m;
        k = m;

        VP     = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFFul);
        VN     = _mm256_setzero_si256();
        NOTONE = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFFul);
        MASK   = _mm256_set_epi64x (0ul,0ul,0ul,1);
        m--;

        i = m / 64;
        while(i){
                bitShiftLeft256ymm(&MASK,64);
                i--;
        }
        bitShiftLeft256ymm(&MASK,m%64);

        for(i = 0; i < n ;i++){
                //X = (B[(int) *t] | VN);

                X = _mm256_or_si256(B[t[i]], VN);
                                //fprintf(stdout,"%lu ", X);
                //D0 = ((VP+(X&VP)) ^ VP) | X ;

                //print_256(X);
                xmm1 = _mm256_and_si256(X, VP);

                xmm2 = add256(0, VP, xmm1);
                //xmm2 = _mm256_add_epi64(VP, xmm1);
                xmm1 = _mm256_xor_si256(xmm2, VP);
                D0 = _mm256_or_si256(xmm1, X);
                //print_256(D0);
                //HN = VP & D0;
                HN =_mm256_and_si256(VP, D0);
                //print_256(HN);
                //HP = VN | ~(VP | D0);

                xmm1 = _mm256_or_si256(VP, D0);
                xmm2 = _mm256_andnot_si256(xmm1, NOTONE);
                HP = _mm256_or_si256(VN, xmm2);
                //print_256(HP);

                //X = HP << 1ul;
                X = HP;
                bitShiftLeft256ymm(&X, 1);

                //print_256(X);
                //VN = X & D0;
                VN= _mm256_and_si256(X, D0);
                //print_256(VN);
                //VP = (HN << 1ul) | ~(X | D0);
                xmm1 = HN;
                bitShiftLeft256ymm(&xmm1, 1);

                xmm2 = _mm256_or_si256(X, D0);

                xmm2 = _mm256_andnot_si256(xmm2, NOTONE);
                //xmm2 = _mm_andnot_si128 (xmm2,NOTONE);
                VP = _mm256_or_si256(xmm1, xmm2);
                //print_256(VP);


                //diff += (HP & MASK) >> m;
                diff += 1- _mm256_testz_si256(HP, MASK);

                ///diff -= (HN & MASK) >> m;
                diff -= 1- _mm256_testz_si256(HN,MASK);

                //fprintf(stdout,"%d ",diff);
                //xmm1 = _mm256_cmpgt_epi64(K, diff);
                k = MACRO_MIN(k, diff);
        }
        return k;
}




/* Carry-correction mask for add256: lane j is 0x8000000000000000 plus bit j
   of m.  Built on the stack, so there is no global state to initialise. */
static KALIGN_TARGET_AVX2 __m256i broadcast_mask(uint32_t m)
{
        uint64_t v[4];
        int j;
        for(j = 0; j < 4; j++){
                v[j] = 0x8000000000000000ull | (uint64_t)((m >> j) & 1u);
        }
        return _mm256_loadu_si256((const __m256i*)v);
}


static KALIGN_TARGET_AVX2 __m256i add256(uint32_t carry, __m256i A, __m256i B)
{
        A = _mm256_xor_si256(A, _mm256_set1_epi64x(0x8000000000000000));
        __m256i s = _mm256_add_epi64(A, B);
        __m256i cv = _mm256_cmpgt_epi64(A, s);
        __m256i mv = _mm256_cmpeq_epi64(s, _mm256_set1_epi64x(0x7fffffffffffffff));
        uint32_t c = _mm256_movemask_pd(_mm256_castsi256_pd(cv));
        uint32_t m = _mm256_movemask_pd(_mm256_castsi256_pd(mv));

        {
                c = m + 2*c; //  lea
                carry += c;
                m ^= carry;
                carry >>= 4;
                m &= 0x0f;
        }
        return _mm256_add_epi64(s, broadcast_mask(m));
}


//----------------------------------------------------------------------------
// bit shift left a 256-bit value using ymm registers
//          __m256i *data - data to shift
//          int count     - number of bits to shift
// return:  __m256i       - carry out bit(s)

static KALIGN_TARGET_AVX2 void bitShiftLeft256ymm (__m256i *data, int count)
{
        __m256i innerCarry, rotate;

        innerCarry = _mm256_srli_epi64 (*data, 64 - count);                        // carry outs in bit 0 of each qword
        rotate     = _mm256_permute4x64_epi64 (innerCarry, 0x93);                  // rotate ymm left 64 bits
        innerCarry = _mm256_blend_epi32 (_mm256_setzero_si256 (), rotate, 0xFC);   // clear lower qword
        *data    = _mm256_slli_epi64 (*data, count);                               // shift all qwords left
        *data    = _mm256_or_si256 (*data, innerCarry);                            // propagate carrys from low qwords
        //carryOut   = _mm256_xor_si256 (innerCarry, rotate);                        // clear all except lower qword
        //return carryOut;
}

#endif /* KALIGN_AVX2 */

/* Plain version of bpm_256: the same 256-bit Myers recurrence on four 64-bit
   words, same pattern cap (255) and same result as bpm_256_avx2. */
uint8_t bpm_256_plain(const uint8_t* t,const uint8_t* p,int n,int m)
{
        uint64_t B[13][4];
        uint64_t VP[4],VN[4],D0[4],HN[4],HP[4],MASK[4];
        int i,j,k,diff;

        if(m > 255){
                m = 255;
        }
        if(m <= 0){
                return 0;
        }
        for(i = 0; i < 13;i++){
                for(j = 0;j < 4;j++){
                        B[i][j] = 0;
                }
        }
        for(i = 0; i < m;i++){
                B[p[i]][i / 64] |= (1ull << (i % 64));
        }
        for(j = 0; j < 4; j++){
                VP[j] = ~0ull;
                VN[j] = 0ull;
                MASK[j] = 0ull;
        }
        MASK[(m - 1) / 64] = 1ull << ((m - 1) % 64);
        diff = m;
        k = m;
        for(i = 0; i < n; i++){
                const uint64_t* Bt = B[t[i]];
                uint64_t carry = 0;
                uint64_t hpc = 0;
                uint64_t hnc = 0;
                int hit_p = 0;
                int hit_n = 0;
                for(j = 0; j < 4; j++){
                        const uint64_t X = Bt[j] | VN[j];
                        const uint64_t a = VP[j] + (X & VP[j]);
                        const uint64_t sum = a + carry;
                        carry = (uint64_t)((a < VP[j]) | (sum < a));
                        D0[j] = (sum ^ VP[j]) | X;
                        HN[j] = VP[j] & D0[j];
                        HP[j] = VN[j] | ~(VP[j] | D0[j]);
                }
                for(j = 0; j < 4; j++){
                        if(HP[j] & MASK[j]) hit_p = 1;
                        if(HN[j] & MASK[j]) hit_n = 1;
                }
                for(j = 0; j < 4; j++){
                        const uint64_t xs = (HP[j] << 1) | hpc;
                        const uint64_t hs = (HN[j] << 1) | hnc;
                        hpc = HP[j] >> 63;
                        hnc = HN[j] >> 63;
                        VN[j] = xs & D0[j];
                        VP[j] = hs | ~(xs | D0[j]);
                }
                diff += hit_p;
                diff -= hit_n;
                k = MACRO_MIN(k, diff);
        }
        return (uint8_t)k;
}

/* bpm_256 for every CPU: the AVX2 kernel when the CPU has it, else plain. */
uint8_t bpm_256(const uint8_t* t,const uint8_t* p,int n,int m)
{
#if KALIGN_AVX2
        if(kalign_simd_on()){
                return bpm_256_avx2(t,p,n,m);
        }
#endif
        return bpm_256_plain(t,p,n,m);
}

/* Kept for source compatibility: the AVX2 kernel keeps no global constants
   any more, so there is nothing to set up.  Safe on any CPU, any number of
   times, from any thread. */
void set_broadcast_mask(void)
{
}




#if defined(__GNUC__) || defined(__clang__)
#define U4_INLINE static inline __attribute__((always_inline)) KALIGN_TARGET_AVX2
#else
#define U4_INLINE static inline
#endif
/* ---- SIMD block wavefront for bpm_block -------------------------------
   4 x 64-bit lanes (NEON: 2 x uint64x2_t, AVX2: __m256i).  Lane b holds block
   b and processes character s - b at step s, so the carry from block b-1 is a
   one-lane shift of the previous step's output.  Per block, the arithmetic is
   exactly bpm_block's (branch-free form).  Blocks that have not started are
   masked during the first steps; finished blocks only feed higher blocks that
   have finished too, and the result is read from the last block while it is
   active.  Used from BPM_SIMD_MIN_BLOCKS blocks on (NEON pays off only for
   long patterns: the scalar integer path is already fast on Apple cores). */
#if KALIGN_NEON
#include <arm_neon.h>
#define BPM_SIMD 1
#define BPM_SIMD_MIN_BLOCKS 6
typedef struct { uint64x2_t lo, hi; } u64x4;
U4_INLINE u64x4 u4_set(uint64_t a, uint64_t b, uint64_t c, uint64_t d)
{
        u64x4 r;
        r.lo = vcombine_u64(vcreate_u64(a), vcreate_u64(b));
        r.hi = vcombine_u64(vcreate_u64(c), vcreate_u64(d));
        return r;
}
U4_INLINE u64x4 u4_dup(uint64_t a) { u64x4 r = { vdupq_n_u64(a), vdupq_n_u64(a) }; return r; }
#define U4_BIN(name, op) U4_INLINE u64x4 name(u64x4 a, u64x4 b) { u64x4 r = { op(a.lo, b.lo), op(a.hi, b.hi) }; return r; }
U4_BIN(u4_and, vandq_u64)
U4_BIN(u4_or, vorrq_u64)
U4_BIN(u4_xor, veorq_u64)
U4_BIN(u4_add, vaddq_u64)
U4_BIN(u4_sub, vsubq_u64)
U4_INLINE u64x4 u4_not(u64x4 a) { return u4_xor(a, u4_dup(~(uint64_t)0)); }
U4_INLINE u64x4 u4_shl1(u64x4 a) { u64x4 r = { vshlq_n_u64(a.lo, 1), vshlq_n_u64(a.hi, 1) }; return r; }
U4_INLINE u64x4 u4_shr63(u64x4 a) { u64x4 r = { vshrq_n_u64(a.lo, 63), vshrq_n_u64(a.hi, 63) }; return r; }
/* [prev3, cur0, cur1, cur2] */
U4_INLINE u64x4 u4_shift_in(u64x4 prev, u64x4 cur)
{
        u64x4 r = { vextq_u64(prev.hi, cur.lo, 1), vextq_u64(cur.lo, cur.hi, 1) };
        return r;
}
/* lanes with mask all-ones take a, others b */
U4_INLINE u64x4 u4_sel(u64x4 mask, u64x4 a, u64x4 b)
{
        u64x4 r = { vbslq_u64(mask.lo, a.lo, b.lo), vbslq_u64(mask.hi, a.hi, b.hi) };
        return r;
}
U4_INLINE uint64_t u4_lane(u64x4 a, int l)
{
        switch(l){
        case 0: return vgetq_lane_u64(a.lo, 0);
        case 1: return vgetq_lane_u64(a.lo, 1);
        case 2: return vgetq_lane_u64(a.hi, 0);
        default: return vgetq_lane_u64(a.hi, 1);
        }
}
#elif KALIGN_AVX2
#define BPM_SIMD 1
#define BPM_SIMD_MIN_BLOCKS 2
typedef __m256i u64x4;
U4_INLINE u64x4 u4_set(uint64_t a, uint64_t b, uint64_t c, uint64_t d)
{
        return _mm256_set_epi64x((long long)d, (long long)c, (long long)b, (long long)a);
}
U4_INLINE u64x4 u4_dup(uint64_t a) { return _mm256_set1_epi64x((long long)a); }
U4_INLINE u64x4 u4_and(u64x4 a, u64x4 b) { return _mm256_and_si256(a, b); }
U4_INLINE u64x4 u4_or(u64x4 a, u64x4 b) { return _mm256_or_si256(a, b); }
U4_INLINE u64x4 u4_xor(u64x4 a, u64x4 b) { return _mm256_xor_si256(a, b); }
U4_INLINE u64x4 u4_add(u64x4 a, u64x4 b) { return _mm256_add_epi64(a, b); }
U4_INLINE u64x4 u4_sub(u64x4 a, u64x4 b) { return _mm256_sub_epi64(a, b); }
U4_INLINE u64x4 u4_not(u64x4 a) { return _mm256_xor_si256(a, _mm256_set1_epi64x(-1)); }
U4_INLINE u64x4 u4_shl1(u64x4 a) { return _mm256_slli_epi64(a, 1); }
U4_INLINE u64x4 u4_shr63(u64x4 a) { return _mm256_srli_epi64(a, 63); }
U4_INLINE u64x4 u4_shift_in(u64x4 prev, u64x4 cur)
{
        const u64x4 c = _mm256_permute4x64_epi64(cur, _MM_SHUFFLE(2, 1, 0, 3));
        const u64x4 p = _mm256_permute4x64_epi64(prev, _MM_SHUFFLE(2, 1, 0, 3));
        return _mm256_blend_epi32(c, p, 0x03);
}
U4_INLINE u64x4 u4_sel(u64x4 mask, u64x4 a, u64x4 b) { return _mm256_blendv_epi8(b, a, mask); }
U4_INLINE uint64_t u4_lane(u64x4 a, int l)
{
        uint64_t t[4];
        _mm256_storeu_si256((__m256i*)t, a);
        return t[l];
}
#else
#define BPM_SIMD 0
#endif


#if BPM_SIMD
/* peq: SIGMA x 16 table as built by bpm_block; B blocks; L = n + padding */
static KALIGN_TARGET_AVX2 int bpm_block_simd(const uint8_t* t, int n, const uint64_t peq[SIGMA][16], int B, int L, int m)
{
        const int V = (B + 3) / 4;
        int k = m;
        /* text with 16 zero characters before and after: char(i) = tp[i] */
        uint8_t tbuf[4096 + 64];
        uint8_t* heap = NULL;
        uint8_t* tp;
        if(L + 32 <= (int)sizeof(tbuf)){
                tp = tbuf + 16;
        }else{
                heap = malloc((size_t)L + 32);
                if(!heap) return -1;
                tp = heap + 16;
        }
        memset(tp - 16, 0, 16);
        if(n > 0) memcpy(tp, t, (size_t)n);
        memset(tp + n, 0, (size_t)(L - n) + 16);

        u64x4 P[4], M[4], SC[4], HP[4], HN[4];
        for(int v = 0; v < V; v++){
                P[v] = u4_dup(~(uint64_t)0);
                M[v] = u4_dup(0);
                SC[v] = u4_set(64 * (4 * v + 1), 64 * (4 * v + 2), 64 * (4 * v + 3), 64 * (4 * v + 4));
                HP[v] = u4_dup(0);
                HN[v] = u4_dup(0);
        }
        const int lv = (B - 1) / 4;
        const int ll = (B - 1) % 4;
        const u64x4 zero = u4_dup(0);
        for(int s = 0; s < L + B - 1; s++){
                u64x4 prevp = zero;
                u64x4 prevn = zero;
                for(int v = 0; v < V; v++){
                        const int b0 = 4 * v;
                        const u64x4 hp = u4_shift_in(prevp, HP[v]);
                        const u64x4 hn = u4_shift_in(prevn, HN[v]);
                        prevp = HP[v];
                        prevn = HN[v];
                        u64x4 Eq = u4_set(peq[tp[s - b0]][b0], peq[tp[s - b0 - 1]][b0 + 1],
                                          peq[tp[s - b0 - 2]][b0 + 2], peq[tp[s - b0 - 3]][b0 + 3]);
                        const u64x4 Pv = P[v];
                        const u64x4 Mv = M[v];
                        const u64x4 Xv = u4_or(Eq, Mv);
                        Eq = u4_or(Eq, hn);
                        const u64x4 Xh = u4_or(u4_xor(u4_add(u4_and(Eq, Pv), Pv), Pv), Eq);
                        u64x4 Ph = u4_or(Mv, u4_not(u4_or(Xh, Pv)));
                        u64x4 Mh = u4_and(Pv, Xh);
                        const u64x4 op = u4_shr63(Ph);
                        const u64x4 on = u4_shr63(Mh);
                        Ph = u4_or(u4_shl1(Ph), hp);
                        Mh = u4_or(u4_shl1(Mh), hn);
                        u64x4 nP = u4_or(Mh, u4_not(u4_or(Xv, Ph)));
                        u64x4 nM = u4_and(Ph, Xv);
                        u64x4 nS = u4_sub(u4_add(SC[v], op), on);
                        if(s < b0 + 3){
                                const u64x4 act = u4_set(s >= b0 ? ~(uint64_t)0 : 0, s >= b0 + 1 ? ~(uint64_t)0 : 0,
                                                         s >= b0 + 2 ? ~(uint64_t)0 : 0, s >= b0 + 3 ? ~(uint64_t)0 : 0);
                                nP = u4_sel(act, nP, Pv);
                                nM = u4_sel(act, nM, Mv);
                                nS = u4_sel(act, nS, SC[v]);
                        }
                        P[v] = nP;
                        M[v] = nM;
                        SC[v] = nS;
                        HP[v] = op;
                        HN[v] = on;
                }
                if(s >= B - 1){
                        const int sc = (int)(int64_t)u4_lane(SC[lv], ll);
                        if(sc < k){
                                k = sc;
                        }
                }
        }
        free(heap);
        return k;
}
#endif

static int bpm_block_impl(const uint8_t *t, const uint8_t *p, int n, int m, int allow_simd)
{
        int32_t w;
        int32_t w_bytes;
        int32_t b_max;
        int32_t maxd;
        int32_t k;
        int32_t W;
        uint64_t ONE = 1;

        if(m > 1024){
                m = 1024;
        }

        w_bytes = sizeof(uint64_t);
        w =  w_bytes * 8;
        b_max = DIV_CEIL(m,w);

        W = w * b_max - m;
        

        k = m;
        maxd = m;

        /* Precompute  */

        uint64_t Peq[13][16] = {
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0},
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0},
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0},
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0},
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0},
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0},
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0},
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0},
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0},
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0},
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0},
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0},
                {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0}
        };


        uint64_t P[16] = {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0};
        uint64_t M[16] = {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0};

        int32_t score[16] = {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0};



        /* m = 1024; *- > 16 */
        /* w_bytes = sizeof(uint64_t); */
        /* w =  w_bytes * 8; */

        /* b_max = DIV_CEIL(m,w); */
        /* LOG_MSG("%d -> %d",m,b_max); */
        /* W = w * b_max - m; */
        /* MMALLOC(score,b_max * sizeof(int32_t)); */
        /* MMALLOC(P ,b_max * w_bytes); */
        /* MMALLOC(M ,b_max * w_bytes); */

        /* MMALLOC(Peq,SIGMA * sizeof(uint64_t *)); */

        /* Peq[c][block] bit j: pattern position block*w + j holds c, or lies
           past the end of the pattern (> m matches anything).  One pass over
           the pattern; this builds exactly the table of the former
           SIGMA x (b_max * w) comparison loop. */
        for (int32_t i = 0; i < m; i++) {
                if (p[i] < SIGMA) {
                        Peq[p[i]][i / w] |= ONE << (i % w);
                }
        }
        for (int32_t i = m; i < b_max * w; i++) {
                for (int c = 0; c < SIGMA; c++) {
                        Peq[c][i / w] |= ONE << (i % w);
                }
        }

#if BPM_SIMD
        if(allow_simd && b_max >= BPM_SIMD_MIN_BLOCKS && kalign_simd_on()){
                int ks = bpm_block_simd(t, n, (const uint64_t (*)[16])Peq, b_max, n + W, m);
                if(ks >= 0){
                        return ks;
                }
        }
#else
        (void)allow_simd;
#endif
        int y = DIV_CEIL(maxd, w) - 1; /* 256 is max edit distance  */
        /* init blocks  */
        for (int b = 0; b <= y; b++) {
                /* bpm_init_block(s,b); */
                P[b] = (uint64_t) -1;//Ones
                M[b] = 0;
                score[b] = (uint32_t)(b + 1) * w;
        }

        /* Main loop.  With maxd = m the last block y = b_max - 1 is always
           active: its score is a semi-global edit distance of the (padded)
           pattern, at most m + W < maxd + w, so the block count never shrinks,
           and it starts at b_max - 1, so it never grows.  The block step is
           branch-free: the horizontal delta handed to the next block is two
           bits, hp (+1) and hn (-1); Ph and Mh are disjoint, so at most one is
           set.  The arithmetic per block is unchanged. */
        for (int i = 0; i < n+W; i++) {
                const uint8_t c = (i >= n) ? 0 : (uint8_t)t[i];
                uint64_t hp = 0;
                uint64_t hn = 0;
                for (int b = 0; b <= y; b++) {
                        const uint64_t Pv = P[b];
                        const uint64_t Mv = M[b];
                        uint64_t Eq = Peq[c][b];
                        const uint64_t Xv = Eq | Mv;
                        Eq |= hn;
                        const uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
                        uint64_t Ph = Mv | (~(Xh | Pv));
                        uint64_t Mh = Pv & Xh;
                        const uint64_t op = Ph >> (w - 1);
                        const uint64_t on = Mh >> (w - 1);
                        Ph = (Ph << 1) | hp;
                        Mh = (Mh << 1) | hn;
                        P[b] = Mh | (~(Xv | Ph));
                        M[b] = Ph & Xv;
                        hp = op;
                        hn = on;
                        score[b] += (int32_t)op - (int32_t)on;
                }
                if(score[y] < k){
                        k = score[y];
                }
        }
        return k;
}

/* The entry point used by BPM(): SIMD kernel when the CPU allows it. */
int bpm_block(const uint8_t *t, const uint8_t *p, int n, int m)
{
        return bpm_block_impl(t, p, n, m, 1);
}

/* Plain kernel only (reference for tests). */
int bpm_block_plain(const uint8_t *t, const uint8_t *p, int n, int m)
{
        return bpm_block_impl(t, p, n, m, 0);
}

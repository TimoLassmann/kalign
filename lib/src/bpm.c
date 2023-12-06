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

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "tldevel.h"
#define BPM_IMPORT
#include "bpm.h"
#include  <stdalign.h>
#include <string.h>

#define SIGMA 13
#define DIV_CEIL(a,b) (a == 0 ? 1 : a/b+(a%b == 0 ? 0 : 1))


#ifdef HAVE_AVX2
#include <immintrin.h>

__m256i BROADCAST_MASK[16];

 void bitShiftLeft256ymm (__m256i *data, int count);
__m256i bitShiftRight256ymm (__m256i *data, int count);

/* taken from Alexander Yee: http://www.numberworld.org/y-cruncher/internals/addition.html#ks_add */
 __m256i add256(uint32_t carry, __m256i A, __m256i B);
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


#ifdef HAVE_AVX2
uint8_t bpm_256(const uint8_t* t,const uint8_t* p,int n,int m)
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




/* Must be called before BPM_256 is!!!  */
void set_broadcast_mask(void)
{
        BROADCAST_MASK[0] =  _mm256_set_epi64x(0x8000000000000000, 0x8000000000000000, 0x8000000000000000, 0x8000000000000000);
        BROADCAST_MASK[1] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000000, 0x8000000000000000, 0x8000000000000001);
        BROADCAST_MASK[2] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000000, 0x8000000000000001, 0x8000000000000000);
        BROADCAST_MASK[3] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000000, 0x8000000000000001, 0x8000000000000001);
        BROADCAST_MASK[4] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000001, 0x8000000000000000, 0x8000000000000000);
        BROADCAST_MASK[5] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000001, 0x8000000000000000, 0x8000000000000001);
        BROADCAST_MASK[6] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000001, 0x8000000000000001, 0x8000000000000000);
        BROADCAST_MASK[7] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000001, 0x8000000000000001, 0x8000000000000001);
        BROADCAST_MASK[8] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000000, 0x8000000000000000, 0x8000000000000000);
        BROADCAST_MASK[9] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000000, 0x8000000000000000, 0x8000000000000001);
        BROADCAST_MASK[10] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000000, 0x8000000000000001, 0x8000000000000000);
        BROADCAST_MASK[11] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000000, 0x8000000000000001, 0x8000000000000001);
        BROADCAST_MASK[12] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000001, 0x8000000000000000, 0x8000000000000000);
        BROADCAST_MASK[13] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000001, 0x8000000000000000, 0x8000000000000001);
        BROADCAST_MASK[14] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000001, 0x8000000000000001, 0x8000000000000000);
        BROADCAST_MASK[15] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000001, 0x8000000000000001, 0x8000000000000001);
}


__m256i add256(uint32_t carry, __m256i A, __m256i B)
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
        return _mm256_add_epi64(s, BROADCAST_MASK[m]);
}


//----------------------------------------------------------------------------
// bit shift left a 256-bit value using ymm registers
//          __m256i *data - data to shift
//          int count     - number of bits to shift
// return:  __m256i       - carry out bit(s)

void bitShiftLeft256ymm (__m256i *data, int count)
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

__m256i bitShiftRight256ymm (__m256i *data, int count)
{
        __m256i innerCarry, carryOut, rotate;


        innerCarry = _mm256_slli_epi64(*data, 64 - count);
        rotate =  _mm256_permute4x64_epi64 (innerCarry, 0x39);
        innerCarry = _mm256_blend_epi32 (_mm256_setzero_si256 (), rotate, 0x3F);
        *data = _mm256_srli_epi64(*data, count);
        *data = _mm256_or_si256(*data,  innerCarry);

        carryOut   = _mm256_xor_si256 (innerCarry, rotate);                        //FIXME: not sure if this is correct!!!
        return carryOut;
}
#endif



int bpm_block(const uint8_t *t, const uint8_t *p, int n, int m)
{
        int32_t w;
        int32_t w_bytes;
        int32_t b_max;
        int32_t maxd;
        int32_t k;
        int32_t W;
        uint64_t HIGH_BIT;
        uint64_t ONE = 1;

        if(m > 1024){
                m = 1024;
        }

        w_bytes = sizeof(uint64_t);
        w =  w_bytes * 8;
        b_max = DIV_CEIL(m,w);
        HIGH_BIT = ONE << (w - 1);

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
        HIGH_BIT = ONE << (w - 1);
        /* MMALLOC(score,b_max * sizeof(int32_t)); */
        /* MMALLOC(P ,b_max * w_bytes); */
        /* MMALLOC(M ,b_max * w_bytes); */

        /* MMALLOC(Peq,SIGMA * sizeof(uint64_t *)); */

        for (int c = 0; c < SIGMA; c++) {
                /* Peq[c] = NULL; */
                /* MMALLOC(Peq[c], b_max*w_bytes); */
                /* Peq[c] = (uint64_t *) malloc(b_max*w_bytes); */

                /* memset(Peq[c], 0, b_max*w_bytes); */
                for (int32_t block = 0; block < b_max; block++) {
                        uint64_t bitPos = (uint64_t) 1;
                        for (int32_t i = block * w; i < (block + 1) * w; ++i) {
                                // fill the remainder after the last block with 1 (>m matches anything)
                                if (i >= m || p[i] == c) {
                                        Peq[c][block] |= bitPos;
                                }
                                /* if (i >= m ){ */
                                /*         LOG_MSG("Am filling %d",i); */
                                /* } */
                                bitPos <<= 1;
                        }
                }
        }

        int y = DIV_CEIL(maxd, w) - 1; /* 256 is max edit distance  */
        /* init blocks  */
        for (int b = 0; b <= y; b++) {
                /* bpm_init_block(s,b); */
                P[b] = (uint64_t) -1;//Ones
                M[b] = 0;
                score[b] = (uint32_t)(b + 1) * w;
        }

        for (int i = 0; i < n+W; i++) {
                /* LOG_MSG(" (%d)",i); */
                uint8_t c = 0;
                if(i >= n){
                        /* seq padding  */
                        c = 0;
                }else{
                        c = (uint8_t)t[i];
                }
                /* LOG_MSG(" c: %d  (%d)", c,i); */
                int carry = 0;

                for (int b = 0; b <= y; b++) {
                        /* LOG_MSG("y: %d c: %d b: %d n: %d m:%d",y, c,b,n,m); */
                        uint64_t Pv = P[b];
                        uint64_t Mv = M[b];
                        uint64_t Eq = Peq[c][b];

                        uint64_t Xv, Xh;
                        uint64_t Ph, Mh;
                        int hIn = carry;
                        int h_out = 0;

                        Xv = Eq | Mv;
                        if (hIn < 0) {
                                Eq |= ONE;
                        }
                        Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

                        Ph = Mv | (~(Xh | Pv));
                        Mh = Pv & Xh;

                        if (Ph & HIGH_BIT) {
                                h_out += 1;
                        }
                        if (Mh & HIGH_BIT) {
                                h_out -= 1;
                        }

                        Ph <<= 1;
                        Mh <<= 1;

                        if (hIn < 0) {
                                Mh |= ONE;
                        } else if (hIn > 0) {
                                Ph |= ONE;
                        }

                        Pv = Mh | (~(Xv | Ph));
                        Mv = Ph & Xv;

                        P[b] = Pv;
                        M[b] = Mv;

                        carry = h_out;
                        score[b] += carry;
                }

                if ((score[y] - carry <= maxd) && (y < (b_max - 1)) && ((Peq[c][y + 1] & ONE) || (carry < 0))) {
                        y += 1;
                        /* bpm_init_block(s,y); */
                        P[y] = (uint64_t) -1;//Ones
                        M[y] = 0;

                        int b = y;

                        uint64_t Pv = P[b];
                        uint64_t Mv = M[b];
                        uint64_t Eq = Peq[c][b];

                        uint64_t Xv, Xh;
                        uint64_t Ph, Mh;
                        int hIn = carry;
                        int h_out = 0;

                        Xv = Eq | Mv;
                        if (hIn < 0) {
                                Eq |= ONE;
                        }
                        Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

                        Ph = Mv | (~(Xh | Pv));
                        Mh = Pv & Xh;

                        if (Ph & HIGH_BIT) {
                                h_out += 1;
                        }
                        if (Mh & HIGH_BIT) {
                                h_out -= 1;
                        }

                        Ph <<= 1;
                        Mh <<= 1;

                        if (hIn < 0) {
                                Mh |= ONE;
                        } else if (hIn > 0) {
                                Ph |= ONE;
                        }

                        Pv = Mh | (~(Xv | Ph));
                        Mv = Ph & Xv;

                        P[b] = Pv;
                        M[b] = Mv;

                        /* carry = h_out; */

                        score[y] = score[y - 1] + w - carry + h_out;//advanceBlock(s,y, c, carry);
                } else {
                        while (score[y] >= (maxd + w)) {
                                if (y == 0) break;
                                y -= 1;
                                /* LOG_MSG("y in loop: %d", y); */
                        }
                }
                if(score[y] < k){
                        /* LOG_MSG("%d", s->score[y]); */
                        k = score[y];
                }


                /* if (y == (b_max - 1) && s->score[y] <= maxd) { */
                /*         assert(i - s->W >= 0); */
                /*         /\* result->push_front(SearchResult(score[y], (uint32_t)(i - W))); *\/ */
                /* } */
        }
        return k;
}



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

#include "bpm.h"
#include  <stdalign.h>

#include "rng.h"

#ifdef HAVE_AVX2
#include <immintrin.h>

__m256i BROADCAST_MASK[16];

 void bitShiftLeft256ymm (__m256i *data, int count);
__m256i bitShiftRight256ymm (__m256i *data, int count);

/* taken from Alexander Yee: http://www.numberworld.org/y-cruncher/internals/addition.html#ks_add */
 __m256i add256(uint32_t carry, __m256i A, __m256i B);
#endif

/* Below are test functions  */
#ifdef BPM_UTEST

#include "alphabet.h"

/* Functions needed for the unit test*/
uint8_t dyn_256(const uint8_t* t,const uint8_t* p,int n,int m);
uint8_t dyn_256_print(const uint8_t* t,const uint8_t* p,int n,int m);
int  mutate_seq(uint8_t* s, int len,int k,int L, struct rng_state* rng);

#ifdef HAVE_AVX2
/* For debugging */
void print_256(__m256i X);
void print_256_all(__m256i X);
#endif

/* The actual test.  */
int bpm_test(void);

int main(int argc, char *argv[])
{

        /* Important set_broadcast_mask has to be called before using bpm_256!!! */
#ifdef HAVE_AVX2
        set_broadcast_mask();
#endif
        RUN(bpm_test());
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int bpm_test(void)
{

        /* idea: start with identical sequences add errors run bpm */
        struct alphabet* alphabet = NULL;
        struct rng_state* rng;
        //long int r;
        int len = 0;
        int i,j,c;
        uint8_t* a = NULL;
        uint8_t* b = NULL;

        int test_iter = 100;

        int calc_errors = 0;
        int total_calc = 0;
        int dyn_score,bpm_score;

        RUNP(rng = init_rng(0));

        RUNP(alphabet = create_alphabet(redPROTEIN));

        char seq[] = "MLLRNSFISQDENDDQLSPTQLGRAAIASSLPPEASLAIFEDLNSASRAIALDTELHMLYLVYFYKNSRAQIIQKIFKIYSIFILKKFKNLEPKFKKKISENITVHITNSIRKKQHFWHVTPINVSVWQECDWHHLFSIFSKLPSDHKRIAKLVGVSEKFILDQLQGRRNDKLLQIHIRFFSALALFDLISEMSIYEVSHKYRIPRGCLQTLQSQSATYAAMIVAFCLRLGWTYLKALLDGFATRLLFGVRSELSELVAIEGIDGQRARILHERGVTCLSHLSACDSSKLAHFLTLAVPYSSSNSNDGLGEWLFGEPRMRVDVAARTLKERARKVLIRRVQELGISVELPKFEENEENIQESCDSGLPDSCEGMEDELEEKENIVKMEEMTKSVTEMSLTDNTISFKSEDDLFKKEIKVEEDEVFIKKEIDEDEEEIVEETVIECLETSLLKLKASTDEVFLRRLSQTFSPIGRSRSILNNSLLEDSFDRPVPRSSIPILNFITPKRESPTPYFEDSFDRPIPGSLPISSSRRKSVLTNIANLDSSRRESINSNASDNNSFDVFVTPPTKSAKEEKRRIAVKHPRVGNIIYSPLTSSPVIKHPKLEINHFYLKDVCHDHNAWNLWTKSSTSTSSCSIRVSDDYTGIAIRTDAGNTFIPLLETFGGENSVPQFFKNISNFYVLKMPLKIFWKQPSPASKYFESFSKCIIPLNTRLEFLKTLAVTVEMYISSMEDAFLIFEKFGIKIFRLKVVRIAAYLNNVIDVEQEENSNFLPILMDRYSILDPEIRKTCSSSLHKAAVEVYSLKPIFEKMCCSGASLQLEMESCQTVLNIFYSGIVFDQALCNSFIYKIRKQIENLEENIWRLAYGKFNIHSSNEVANVLFYRLGLIYPETSGCKPKLRHLPTNKLILEQMNTQHPIVGKILEYRQIQHTLTQCLMPLAKFIGRIHCWFEMCTSTGRILTSVPNLQNVPKRISSDGMSARQLFIANSENLLIGADYKQLELRVLAHLSNDSNLVNLITSDRDLFEELSIQWNFPRDAVKQLCYGLIYGMGAKSLSELTRMSIEDAEKMLKAFFAMFPGVRSYINETKEKVCKEEPISTIIGRRTIIKASGIGEERARIERVAVNYTIQGSASEIFKTAIVDIESKIKEFGAQIVLTIHDEVLVECPEIHVAAASESIENCMQNALSHLLRVPMRVSMKTGRSWADLK";
        len = strlen(seq);
        if(len > 256){
                len = 256;
        }

        MMALLOC(a , sizeof(uint8_t) * len) ;
        MMALLOC(b , sizeof(uint8_t) * len) ;

        for(i = 0;i < len;i++){
                a[i] = alphabet->to_internal[(int)seq[i]];
                b[i] =a[i];
                //fprintf(stdout,"%d %d\n", a[i],b[i]);
        }
        //fprintf(stdout,"LEN: %d\n", len);


        //LOG_MSG("Testing correctness of serial bpm.");
//        len = 63;
        calc_errors = 0;
        total_calc = 0;
        len =63; //could fix;
        for(i = 0; i < 10;i++){
                for (j =0 ; j < test_iter; j++){
                        RUN(mutate_seq(b,len,i,alphabet->L,rng));
                        dyn_score = dyn_256(a,b,len,len);
                        bpm_score = bpm(a,b,len,len);

                        if( abs( dyn_score - bpm_score) != 0){
                                fprintf(stdout,"Scores differ: %d (dyn) %d (bpm) (%d out of %d)\n", dyn_score,bpm_score, calc_errors , total_calc);
                                calc_errors++;
                        }
                        /* restore sequence b */
                        for(c = 0;c < len;c++){
                                b[c] = a[c];
                        }

                        total_calc++;
                }

        }
        //fprintf(stdout,"%d errors out of %d\n", calc_errors , total_calc);
        //LOG_MSG("Testing correctness of AVX bpm.");
        len = strlen(seq);
        if(len > 255){
                len = 255;
        }
//        len = 63;
        calc_errors = 0;
        total_calc = 0;
        for(i = 0; i < 10;i++){
                for (j =0 ; j < test_iter; j++){
                        RUN(mutate_seq(b,len,i,alphabet->L,rng));
                        dyn_score = dyn_256(a,b,len,len);
#ifdef HAVE_AVX2
                        bpm_score = bpm_256(a,b,len,len);
#else
                        bpm_score = dyn_score;
#endif
                        if( abs( dyn_score - bpm_score) != 0){
                                fprintf(stdout,"Scores differ: %d (dyn) %d (bpm) (%d out of %d)\n", dyn_score,bpm_score, calc_errors , total_calc);
                                calc_errors++;
                        }
                        /* restore sequence b */
                        for(c = 0;c < len;c++){
                                b[c] = a[c];
                        }

                        total_calc++;
                }

        }
        //fprintf(stdout,"%d errors out of %d\n", calc_errors , total_calc);


        len = strlen(seq);
        if(len > 255){
                len = 255;
        }




        double dyn_timing;
        double bpm_timing;
        int timing_iter = 10000;
        fprintf(stdout,"Dyn\tAVX\tDYN/AVX\n");
        DECLARE_TIMER(t);
        for(i = 0; i < 100;i+=10){
                RUN(mutate_seq(b,len,i,alphabet->L,rng));
                START_TIMER(t);
                for(j = 0; j < timing_iter;j++){
                        dyn_score = dyn_256(a,b,len,len);
                }
                STOP_TIMER(t);
                dyn_timing = GET_TIMING(t);


                START_TIMER(t);
#ifdef HAVE_AVX2
                for(j = 0; j < timing_iter;j++){
                        bpm_score = bpm_256(a,b,len,len);
                }
#else
                bpm_score = dyn_score;
#endif
                STOP_TIMER(t);

                bpm_timing = GET_TIMING(t);

                ASSERT(dyn_score == bpm_score, "Scores differ: %d %d.",dyn_score, bpm_score);
                fprintf(stdout,"%f\t%f\t%f\n",dyn_timing,bpm_timing,  dyn_timing / bpm_timing);

                /* restore seq */
                for(j = 0;j < len;j++){
                        b[j] = a[j];
                }
        }

        MFREE(a);
        MFREE(b);
        MFREE(rng);
        return OK;
ERROR:
        return FAIL;
}

int mutate_seq(uint8_t* s, int len,int k,int L, struct rng_state* rng)
{
        int i,j;
        int r;

        for(i = 0; i < k;i++){
                r = tl_random_int(rng,len);

                j = r;
                r = tl_random_int(rng,L);
                s[j] = r;
        }
        return OK;

}


uint8_t dyn_256(const uint8_t* t,const uint8_t* p,int n,int m)
{
        uint8_t* prev = NULL;
        uint8_t* cur = NULL;

        uint8_t* tmp = NULL;
        int i,j,c;

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
                cur[0] = prev[0];
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

uint8_t dyn_256_print(const uint8_t* t,const uint8_t* p,int n,int m)
{
        uint8_t* prev = NULL;
        uint8_t* cur = NULL;

        uint8_t* tmp = NULL;
        int i,j,c;

        MMALLOC(prev, sizeof(uint8_t)* 257);
        MMALLOC(cur, sizeof(uint8_t)* 257);
        cur[0] = 0;
        fprintf(stdout,"%d ", cur[0]);
        for(j = 1; j <= m;j++){
                cur[j] = cur[j-1] +1;
                fprintf(stdout,"%d ", cur[j]);
        }
        fprintf(stdout,"\n");
        tmp  = cur;
        cur = prev;
        prev = tmp;

        for(i = 1; i <= n;i++){
                cur[0] = prev[0];
                fprintf(stdout,"%d ", cur[0]);
                for(j = 1; j < m;j++){
                        c = 1;
                        if(t[i-1] == p[j-1]){
                                c = 0;

                        }
                        cur[j] = prev[j-1] +c ;

                        cur[j] = MACRO_MIN(cur[j], prev[j]+1);
                        cur[j] = MACRO_MIN(cur[j], cur[j-1]+1);
                        fprintf(stdout,"%d ", cur[j]);
                }
                j =m;
                c = 1;
                if(t[i-1] == p[j-1]){
                        c = 0;

                }
                cur[j] = prev[j-1] +c ;

                cur[j] = MACRO_MIN(cur[j], prev[j]);
                cur[j] = MACRO_MIN(cur[j], cur[j-1]+1);



                fprintf(stdout,"%d ", cur[j]);

                fprintf(stdout,"\n");
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

#ifdef HAVE_AVX2
void print_256(__m256i X)
{
        alignas(32) uint64_t debug[4];
        _mm256_store_si256( (__m256i*)& debug,X);
        fprintf(stdout,"%lu ", debug[0]);
}


void print_256_all(__m256i X)
{
        alignas(32) uint64_t debug[4];
        _mm256_store_si256( (__m256i*)& debug,X);
        int i;
        for(i = 0; i < 4;i++){
                fprintf(stdout,"%lu ", debug[i]);
        }
        fprintf(stdout,"\n");
}
#endif

#endif


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

        for(i = 0; i < m;i++){
                B[p[i]] |= (1ul << i);
        }

        VP = (1ul << (m))-1 ;

        VN = 0ul;

        m--;
        MASK = 1ul << (m);

        //fprintf(stdout,"BEGINNING\t%lu %lu\n",VN,VP);
        for(i = 0; i < n;i++){
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

        alignas(32)  uint32_t f[21][8];
        //int ALIGNED_(64) f[8];
        if(m > 255){
                m = 255;
        }

        for(i = 0; i < 21;i++){
                for(j = 0;j < 8;j++){
                        f[i][j] =0u;
                }
        }

        for(i = 0; i < m;i++){
                f[p[i]][i/32] |= (1 << (i % 32));
        }

        for(i = 0; i < 13;i++){
                B[i] = _mm256_load_si256((__m256i const*) &f[i]);
        }

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

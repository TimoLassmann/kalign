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

#include "tldevel.h"
#include "tlrng.h"

#ifdef HAVE_AVX2
#include <immintrin.h>
#endif

#include "misc.h"
#include  <stdalign.h>
#include <string.h>
#include <immintrin.h>
#ifdef ITEST_MISC


#include "alphabet.h"



int hash_test(void);

//int  mutate_seq(uint8_t* s, int len,int k,int L, struct rng_state* rng);
int main(int argc, char *argv[])
{
        //int a = 2345;

        RUN(hash_test());

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}




int hash_test(void)
{
        struct alphabet* a = NULL;
        unsigned short hash = 0;
        unsigned short rolling = 0;
        int kmer_len = 10;
        int len = 0;
        int i;
        uint8_t* internal = NULL;

        RUNP(a = create_alphabet(ALPHA_defPROTEIN));

        char seq[] = "GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE";

        len = strlen(seq);
        MMALLOC(internal , sizeof(uint8_t) * len) ;

        for(i = 0;i < len;i++){
                internal[i] = a->to_internal[(int)seq[i]];
                //fprintf(stdout,"%c %d\n", seq[i], internal[i]);
        }
        fprintf(stdout,"LEN: %d\n", len);

        hash = circ_hash(internal , kmer_len);
        rolling = hash;
        for(i = 1; i < len - kmer_len;i++){
                rolling =circ_hash_next(internal+i, kmer_len, internal[i-1], rolling);
                hash = circ_hash(internal+i , kmer_len);
                //fprintf(stdout,"%d %d\n", hash, rolling);
        }
        MFREE(a);

        MFREE(internal);
        return OK;
ERROR:
        return FAIL;
}

#endif



int byg_count(char* pattern,char*text)
{
        int Tc;
        int count = 0;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }

        int m = strlen(pattern);
        int n = strlen (text);
        int mb = (1 << (m-1));

        for (i= 0;i < m;i++){
                T[(int)pattern[i]] |= (1 << i);
        }

        for (i = 0;i < n;i++){
                s <<= 1;
                s |= 1;
                Tc = T[(int)text[i]];
                s &= Tc;
                if(s & mb){
                        count++;
                }
        }
        return count;
}

int byg_end(char* pattern,char*text)
{
        int Tc;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }

        int m = strlen(pattern);
        int n = strlen (text);
        int mb = (1 << (m-1));

        for (i= 0;i < m;i++){
                T[(int)pattern[i]] |= (1 << i);
        }

        for (i = 0;i < n;i++){
                s <<= 1;
                s |= 1;
                if(!text[i]){
                        return -1;
                }
                Tc = T[(int)text[i]];
                s &= Tc;
                if(s & mb){
                        return i+1;
                }
        }
        return -1;
}


int byg_start(char* pattern,char*text)
{
        //LOG_MSG("searching for: %s",pattern);
        int Tc;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }

        int m = strlen(pattern);
        int n = strlen(text);
        int mb = (1 << (m-1));

        for (i= 0;i < m;i++){
                T[(int)pattern[i]] |= (1 << i);
        }

        for (i = 0;i < n;i++){
                s <<= 1;
                s |= 1;
                Tc = T[(int)text[i]];
                s &= Tc;
                if(s & mb){
                        return i-m+1;
                }
                //fprintf(stdout,"%c %d %d\n", text[i], i, s);
                if(i == 100){
                        return -1;
                }
        }
        return -1;
}









// (c) 2017 Johannes Soeding & Martin Steinegger, Gnu Public License version 3
// Rotate left macro: left circular shift by numbits within 16 bits
#define RoL(val, numbits) (val << numbits) ^ (val >> (16 - numbits))


// Transform each letter x[i] to a fixed random number RAND[x[i]]
// to ensure instantaneous mixing into the 16 bits
// Do XOR with RAND[x[i]] and 5-bit rotate left for each i from 1 to k
uint16_t circ_hash(const uint8_t* x, const uint8_t length)
{
        const uint16_t RAND[21] = {0x4567, 0x23c6, 0x9869, 0x4873, 0xdc51, 0x5cff, 0x944a, 0x58ec,
                                   0x1f29, 0x7ccd, 0x58ba, 0xd7ab, 0x41f2, 0x1efb, 0xa9e3, 0xe146,
                                   0x007c, 0x62c2, 0x0854, 0x27f8, 0x231b};// 16 bit random numbers
        register uint16_t h = 0x0;
        h = h^ RAND[x[0]];// XOR h and ki
        for (register uint8_t i = 1; i < length; ++i){
                h = RoL(h, 5);
                h ^= RAND[x[i]];// XOR h and ki
        }
        return h;
}

// Rolling hash variant for previous hash function:
// Computes hash value for next key x[0:length-1] from previous hash value
// hash( x[-1:length-2] ) and x_first = x[-1]

uint16_t circ_hash_next(const uint8_t * x,const uint8_t length,const uint8_t x_first, uint16_t h)
{
        const uint16_t RAND[21] = {0x4567, 0x23c6, 0x9869, 0x4873, 0xdc51, 0x5cff, 0x944a, 0x58ec,
                                   0x1f29, 0x7ccd, 0x58ba, 0xd7ab, 0x41f2, 0x1efb, 0xa9e3, 0xe146,
                                   0x007c, 0x62c2, 0x0854, 0x27f8, 0x231b};// 16 bit random numbers
// undo INITIAL_VALUE and first letter x[0] of old key
        h ^= RoL(RAND[x_first], (5*(length-1)) % 16);
// circularly permute all letters x[1:length-1] to 5 positions to left
        h =  RoL(h, 5);// add new, last letter of new key x[1:length]
        h ^= RAND[x[length-1]];
        return h;
}

int shuffle_arr_r(int* arr,int n, struct rng_state* rng)
{
        int r;
        int i,j;
        int tmp;
        for (i = 0; i < n - 1; i++) {
                r = tl_random_int(rng,n);
                //lrand48_r(randBuffer,&r);
                j = i +  r % (n-i);
                tmp = arr[j];
                arr[j] = arr[i];
                arr[i] = tmp;
        }
        return OK;
}

char* basename(const char* name)
{
        int i= 0;
        int c = 0;

        while(1){
                if(name[i] == '/'){
                        c = i+1;
                }
                if(!name[i]){
                        break;
                }
                i++;
        }
        return (char*)(name +c);
}

#include <immintrin.h>

#include "misc.h"
#include  <stdalign.h>

static __m256i bitShiftLeft256ymm (__m256i *data, int count);

static __m256i bitShiftRight256ymm (__m256i *data, int count);

#ifdef ITEST_MISC

#include "alphabet.h"



/*
__m256i shift_right256(__m256i A,const int N)
{
        if     ( N ==  0 ) return A;
        else if( N <  16 ) return _mm256_alignr_epi8(_mm256_permute2x128_si256(A, A, _MM_SHUFFLE(2, 0, 0, 1)), A, N);
        else if( N == 16 ) return                    _mm256_permute2x128_si256(A, A, _MM_SHUFFLE(2, 0, 0, 1));
        else               return _mm256_srli_si256 (_mm256_permute2x128_si256(A, A, _MM_SHUFFLE(2, 0, 0, 1)), N - 16);
}

 __m256i shift_left256(__m256i A,const  int N)
{
        if     ( N ==  0 ) return A;
        else if( N <  16 ) return _mm256_alignr_epi8(A, _mm256_permute2x128_si256(A, A, _MM_SHUFFLE(0, 0, 2, 0)), (uint8_t)(16 - N));
        else if( N == 16 ) return                       _mm256_permute2x128_si256(A, A, _MM_SHUFFLE(0, 0, 2, 0));
        else               return _mm256_slli_si256 (   _mm256_permute2x128_si256(A, A, _MM_SHUFFLE(0, 0, 2, 0)), (uint8_t)(N - 16));
}
*/

static __m256i bitShiftLeft256ymm (__m256i *data, int count);
//----------------------------------------------------------------------------
// bit shift left a 256-bit value using ymm registers
//          __m256i *data - data to shift
//          int count     - number of bits to shift
// return:  __m256i       - carry out bit(s)

static __m256i bitShiftLeft256ymm (__m256i *data, int count)
{
        __m256i innerCarry, carryOut, rotate;

        innerCarry = _mm256_srli_epi64 (*data, 64 - count);                        // carry outs in bit 0 of each qword
        rotate     = _mm256_permute4x64_epi64 (innerCarry, 0x93);                  // rotate ymm left 64 bits
        innerCarry = _mm256_blend_epi32 (_mm256_setzero_si256 (), rotate, 0xFC);   // clear lower qword
        *data    = _mm256_slli_epi64 (*data, count);                               // shift all qwords left
        *data    = _mm256_or_si256 (*data, innerCarry);                            // propagate carrys from low qwords
        carryOut   = _mm256_xor_si256 (innerCarry, rotate);                        // clear all except lower qword
        return carryOut;
}

static __m256i bitShiftRight256ymm (__m256i *data, int count)
{
        __m256i innerCarry, carryOut, rotate;


        innerCarry = _mm256_slli_epi64(*data, 64 - count);

        rotate =  _mm256_permute4x64_epi64 (innerCarry, 0x39);

        innerCarry = _mm256_blend_epi32 (_mm256_setzero_si256 (), rotate, 0x3F);

        *data = _mm256_srli_epi64(*data, count);

        *data = _mm256_or_si256(*data,  innerCarry);

        carryOut   = _mm256_xor_si256 (innerCarry, rotate);                        // clear all except lower qword
        return carryOut;
}



int hash_test(void);
int bpm_test(void);
int  mutate_seq(uint8_t* s, int len,int k,int L, struct drand48_data* b);
uint8_t bpm_256(const uint8_t* t,const uint8_t* p,int n,int m);
int main(int argc, char *argv[])
{
        int a = 2345;

        RUN(hash_test());
        RUN(bpm_test());
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int bpm_test(void)
{

        /* idea: start with identical sequences add errors run bpm */
        struct alphabet* alphabet = NULL;
        struct drand48_data randBuffer;
        long int r;
        int len = 0;
        int i,j,c;
        uint_fast8_t* a = NULL;
        uint_fast8_t* b = NULL;

        srand48_r(time(NULL), &randBuffer);

        RUNP(alphabet = create_alphabet(defPROTEIN));

        char seq[] = "GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE";

        len = strlen(seq);
        MMALLOC(a , sizeof(uint_fast8_t) * len) ;
        MMALLOC(b , sizeof(uint_fast8_t) * len) ;

        for(i = 0;i < len;i++){
                a[i] = alphabet->to_internal[(int)seq[i]];
                b[i] =a[i];
                fprintf(stdout,"%d %d\n", a[i],b[i]);
        }
        fprintf(stdout,"LEN: %d\n", len);

        bpm_256(a,b,len,len);
        exit(0);
        len = 8;
        uint64_t res= 0;

        for(i = 0; i < 10;i++){
                RUN(mutate_seq(b,len,i,alphabet->L,&randBuffer));
                res= 0;
                for(j = 0; j < len;j++){
                        if(a[j] != b[j]){
                                res++;
                        }
                }
                fprintf(stdout,"%d\tdiff:%ld",i,res);
                res= bpm(a,b,len,len);
                fprintf(stdout,"\t%ld\n",res);
                for(j = 0; j < len;j++){
                        fprintf(stdout,"%*d",3,a[j]);
                }


                fprintf(stdout,"\n");
                for(j = 0; j < len;j++){
                        fprintf(stdout,"%*d",3,b[j]);
                }
                fprintf(stdout,"\n");

                /* restore seq */
                for(j = 0;j < len;j++){
                        b[j] = a[j];
                }
        }
                //lrand48_r(randBuffer, &r);
        //r = r % num_samples;

        MFREE(a);
        MFREE(b);
        return OK;
ERROR:
        return FAIL;
}


int  mutate_seq(uint8_t* s, int len,int k,int L, struct drand48_data* b)
{
        int i,j;
        long int r;

        for(i = 0; i < k;i++){
                lrand48_r(b, &r);
                r = r % len;
                j = r;
                lrand48_r(b, &r);
                r = r % L;
                s[j] = r;
        }
        return OK;

}

int hash_test(void)
{
        struct alphabet* a = NULL;
        unsigned short hash = 0;
        unsigned short rolling = 0;
        int kmer_len = 10;
        int len = 0;
        int i;
        uint_fast8_t* internal = NULL;

        RUNP(a = create_alphabet(defPROTEIN));

        char seq[] = "GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE";

        len = strlen(seq);
        MMALLOC(internal , sizeof(uint_fast8_t) * len) ;

        for(i = 0;i < len;i++){
                internal[i] = a->to_internal[(int)seq[i]];
                fprintf(stdout,"%c %d\n", seq[i], internal[i]);
        }
        fprintf(stdout,"LEN: %d\n", len);

        hash = circ_hash(internal , kmer_len);
        rolling = hash;
        for(i = 1; i < len - kmer_len;i++){
                rolling =circ_hash_next(internal+i, kmer_len, internal[i-1], rolling);
                hash = circ_hash(internal+i , kmer_len);
                fprintf(stdout,"%d %d\n", hash, rolling);
        }


        MFREE(internal);
        return OK;
ERROR:
        return FAIL;
}

#endif

int byg_detect(uint8_t* text,int n)
{
        int Tc;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }
        int mb = 1;
        //char *unique_aa = "EFILPQXZ";//permissiv
        //ABCDEFGHIJKLMNOPQRSTUVWXYZ
        char *unique_aa = "BDEFHIJKLMNOPQRSVWYZ";//restrictive
        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        for (i= 0;i < 20;i++){
                T[(int)aacode[unique_aa[i]-65]] |= 1;
        }
        for (i = 0;i < n;i++){
                //	fprintf(stderr,"%d\n",text[i]);
                if(text[i] != -1){
                        s <<= 1;
                        s |= 1;
                        Tc = T[text[i]];
                        s &= Tc;
                        if(s & mb){
                                return 0;
                        }
                }
        }
        return 1;
}


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
        }
        return -1;
}



uint8_t bpm(const uint8_t* t,const uint8_t* p,int n,int m)
{
        register uint64_t VP,VN,D0,HN,HP,X;
        register uint64_t i;//,c;
        uint64_t MASK = 0;
        uint64_t diff;
        uint64_t B[21];
        uint8_t k = m;

        if(m > 64){
                m = 64;
        }
        diff = m;

        for(i = 0; i < 21;i++){
                B[i] = 0;
        }

        for(i = 0; i < m;i++){
                B[p[i]] |= (1ul << i);
        }

        VP = 0xFFFFFFFFFFFFFFFFul;
        VN = 0ul;
        m--;
        MASK = 1ul << m;
        for(i = 0; i < n;i++){
                X = (B[t[i]] | VN);
                D0 = ((VP+(X&VP)) ^ VP) | X ;
                HN = VP & D0;
                HP = VN | ~(VP | D0);
                X = HP << 1ul;
                VN = X & D0;
                VP = (HN << 1ul) | ~(X | D0);
                diff += (HP & MASK) >> m;
                diff -= (HN & MASK) >> m;
                if(diff < k){
                        k = diff;

                }
        }
        return k;
}



uint8_t bpm_256(const uint8_t* t,const uint8_t* p,int n,int m)
{
        //sss__m256i VP,VN,D0,HN,HP,X;
        int i,j,diff,k,c;

        __m256i MASK;

        __m256i B[21];




        alignas(32)  uint32_t f[21][8];
        //int ALIGNED_(64) f[8];
        if(m > 256){
                m = 256;
        }
        m = 10;
        for(i = 0; i < 21;i++){
                for(j = 0;j < 8;j++){
                        f[i][j] =0u;
                }

        }

        for(i = 0; i < m;i++){


                fprintf(stdout,"pos:%d %*d,",i,3,p[i]);
                fprintf(stdout,"\tb: %d\t",f[p[i]][i/32]);
                f[p[i]][i/32] |= (1 << (i % 32));
                fprintf(stdout,"set bit %d in arr %d[%d]", i%32,p[i],i/32);
                fprintf(stdout,"\ta: %d\n",f[p[i]][i/32]);
        }
        fprintf(stdout,"\n");
        for(i = 0; i < 21;i++){
                fprintf(stdout,"Letter: %d ",i);
                for(j = 0;j < 8;j++){
                        for(c = 0; c < 32;c++){
                                if(f[i][j] & (1 << c)){
                                        fprintf(stdout,"%d(%d),",c,j);
                                }

                        }
                }
                fprintf(stdout,"\n");

        }
        diff = m;

        for(i = 0; i < 21;i++){

                //B[i] = _mm256_setzero_si256();
                B[i] = _mm256_load_si256((__m256i const*) &f[i]);
                /*B[i] = _mm256_set_epi32(f[i][0],
                                 f[i][1],
                                 f[i][2],
                                 f[i][3],
                                 f[i][4],
                                 f[i][5],
                                 f[i][6],
                                 f[i][7]     );*/
        }


        for(i = 0; i < 21;i++){
                for(j = 0; j < 8;j++){
                        f[i][j] = 0;

                }
        }
        for(i = 0; i < 21;i++){
                //left_shift_by_one_256(B[i]);
                //bitShiftLeft256ymm
                //B[i] = shift_left256(B[i],1);
                for(j = 0; j < 137;j++){
                        bitShiftLeft256ymm(&B[i],1);
                }

                for(j = 0; j < 137;j++){
                        bitShiftRight256ymm(&B[i],1);
                }


        }

        for(i = 0; i< 21;i++){
                _mm256_store_si256((__m256i*) &f[0],B[i]);
                fprintf(stdout,"Letter %d: ",i);
                for(j = 7;j >= 0;j--){
                        //fprintf(stdout,"%x,",f[0][j]);
                        for(c =0;c < 32;c++){
                                if(f[0][j] & (1 << c)){
                                        fprintf(stdout,"%d(%d),",c,j);
                                }else{
                                        //fprintf(stdout,"----");
                                }
                                }
                }
                fprintf(stdout,"\n");
        }


/*
        rintf(stderr,"%c",*t + 65);
		X = _mm_or_si128 (*(nuc_p +( (int)(*t)  & 0x3u) ) , VN);
		//X = (B[(int) *t] | VN);
		xmm1 = _mm_and_si128(X, VP);
		xmm2 = _mm_add_epi32(VP ,xmm1);
		xmm1 = _mm_xor_si128 (xmm2, VP);
		D0 = _mm_or_si128(xmm1, X);
		//D0 = ((VP+(X&VP)) ^ VP) | X ;
		HN = _mm_and_si128(VP, D0);
		//HN = VP & D0;
		xmm1 = _mm_or_si128(VP, D0);
		xmm2 = _mm_andnot_si128 (xmm1,NOTONE);
		HP = _mm_or_si128(VN, xmm2);
		//HP = VN | ~(VP | D0);
		X = _mm_slli_epi32(HP,1);
		//X = HP << 1ul;
		VN = _mm_and_si128(X, D0);
		//VN = X & D0;
		xmm1 = _mm_slli_epi32(HN,1);
		xmm2 = _mm_or_si128(X, D0);
		xmm2 = _mm_andnot_si128 (xmm2,NOTONE);
		VP = _mm_or_si128(xmm1, xmm2);
		//VP = (HN << 1ul) | ~(X | D0);
		xmm1 = _mm_and_si128(HP, MASK);
		xmm2 = _mm_cmpgt_epi32(xmm1, zero);
		diff = _mm_add_epi32(diff , _mm_and_si128( xmm2, one));
		//diff += (HP & MASK) >> m;
		xmm1 = _mm_and_si128(HN, MASK);
		xmm2 = _mm_cmpgt_epi32(xmm1, zero);
		diff = _mm_sub_epi32(diff,  _mm_and_si128( xmm2, one));
		//diff -= (HN & MASK) >> m;
		xmm1 = _mm_cmplt_epi32(diff, K);
		xmm2 = _mm_and_si128(xmm1, diff);
		K = _mm_or_si128(xmm2, _mm_andnot_si128  (xmm1,K));
		xmm2 = _mm_and_si128(xmm1, _mm_set1_epi32(i));
		POS = _mm_or_si128(xmm2, _mm_andnot_si128  (xmm1,POS));
		t--;*/


        return OK;
}

/*
int validate_bpm_sse(struct qs_struct* qs,int* assignment,unsigned char* t,int n,int num,long int offset)
{
	int i,j;
	int len = 0;

	struct seq_info* si = 0;
	unsigned long int new;
	unsigned int _MM_ALIGN16 nuc[16];
	int _MM_ALIGN16 positions[4];
	int _MM_ALIGN16 errors[4];

	unsigned int _MM_ALIGN16 lengths[4];

	__m128i VP,VN,D0,HN,HP,X,MASK,K,NOTONE,POS,diff,zero,one;
	__m128i* nuc_p;
	__m128i xmm1,xmm2;
	//long int MASK = 0;

	for(i = 0; i < 16;i++){
		nuc[i] = 0ul;
	}

	for(i = 0; i < num;i++){
		si = qs->seq_info[abs(assignment[i])];
		len = si->len;
		if(len > 31){
			len = 31;
		}

		lengths[i] = len;


		if(assignment[i] >= 0){
			for(j = 0; j < len;j++){
				nuc[((int)(si->seq[j] & 0x3u) << 2) + i] |=  (1ul << (unsigned long int)(len-1-j));
	//			fprintf(stderr,"%c",si->seq[j]+65);
			}
		}else{
			for(j = 0; j < len;j++){
				nuc[((int)(si->reverse_seq[j] & 0x3u) << 2) + i] |=  (1ul << (unsigned long int)(len-1-j));
	//			fprintf(stderr,"%c",si->reverse_seq[j] + 65);
			}
		}
	//	fprintf(stderr,"\n" );
	}
	nuc_p = (__m128i*) nuc;
	//S = _mm_set1_epi32(0);
	zero = _mm_set1_epi32(0);
	one = _mm_set1_epi32(1);
	diff =  _mm_load_si128 ( (__m128i*) lengths );  // _mm_set1_epi32(m);
	VP =  _mm_set1_epi32(0xFFFFFFFFu);
	VN =  _mm_set1_epi32(0);
	NOTONE =  _mm_set1_epi32(0xFFFFFFFF);
	K =  _mm_set1_epi32(0x7FFFFFFF);
	POS = _mm_set1_epi32(0);

	//VP = 0xFFFFFFFFFFFFFFFFul;
	//VN = 0ul;

	for(i = 0; i< 4;i++){
		lengths[i]--;
		lengths[i] = 1 << lengths[i];
	//	fprintf(stderr,"%d	%d	LEN:%d\n",i,num,lengths[i]);
	}

	// m--;
	MASK =  _mm_load_si128 ( (__m128i*) lengths ); //  _mm_set1_epi32(1ul << m);
	for(i = 0; i < n ;i++){
		//fprintf(stderr,"%c",*t + 65);
		X = _mm_or_si128 (*(nuc_p +( (int)(*t)  & 0x3u) ) , VN);
		//X = (B[(int) *t] | VN);
		xmm1 = _mm_and_si128(X, VP);
		xmm2 = _mm_add_epi32(VP ,xmm1);
		xmm1 = _mm_xor_si128 (xmm2, VP);
		D0 = _mm_or_si128(xmm1, X);
		//D0 = ((VP+(X&VP)) ^ VP) | X ;
		HN = _mm_and_si128(VP, D0);
		//HN = VP & D0;
		xmm1 = _mm_or_si128(VP, D0);
		xmm2 = _mm_andnot_si128 (xmm1,NOTONE);
		HP = _mm_or_si128(VN, xmm2);
		//HP = VN | ~(VP | D0);
		X = _mm_slli_epi32(HP,1);
		//X = HP << 1ul;
		VN = _mm_and_si128(X, D0);
		//VN = X & D0;
		xmm1 = _mm_slli_epi32(HN,1);
		xmm2 = _mm_or_si128(X, D0);
		xmm2 = _mm_andnot_si128 (xmm2,NOTONE);
		VP = _mm_or_si128(xmm1, xmm2);
		//VP = (HN << 1ul) | ~(X | D0);
		xmm1 = _mm_and_si128(HP, MASK);
		xmm2 = _mm_cmpgt_epi32(xmm1, zero);
		diff = _mm_add_epi32(diff , _mm_and_si128( xmm2, one));
		//diff += (HP & MASK) >> m;
		xmm1 = _mm_and_si128(HN, MASK);
		xmm2 = _mm_cmpgt_epi32(xmm1, zero);
		diff = _mm_sub_epi32(diff,  _mm_and_si128( xmm2, one));
		//diff -= (HN & MASK) >> m;
		xmm1 = _mm_cmplt_epi32(diff, K);
		xmm2 = _mm_and_si128(xmm1, diff);
		K = _mm_or_si128(xmm2, _mm_andnot_si128  (xmm1,K));
		xmm2 = _mm_and_si128(xmm1, _mm_set1_epi32(i));
		POS = _mm_or_si128(xmm2, _mm_andnot_si128  (xmm1,POS));
		t--;
	}
	//fprintf(stderr,"\n");
	_mm_store_si128 ((__m128i*) positions, POS);
	//fprintf(stderr,"%d	%d
	%d	%d	",out1[0],out1[1],out1[2],out1[3]);
	_mm_store_si128 ((__m128i*) errors, K);
	//fprintf(stderr,"%d	%d	%d	%d\n",out2[0],out2[1],out2[2],out2[3]);

	for(i = 0; i < num;i++){
		si = qs->seq_info[abs(assignment[i])];
		//len = si->len;
		//lengths[i] = len;
	//	fprintf(stderr,"%d	%d	%d	%d	%d\n",  assignment[i] , si->len,num,out2[i],out1[i] + offset);
		if(assignment[i] >= 0){
			new = ((unsigned long int)(si->len - errors[i])) << 56ul;
			new |= ((unsigned long int) ((unsigned long int)offset - (unsigned long int) positions[i]) << 1ul);
			BinaryInsertionSortOne(si->match_units,LIST_STORE_SIZE ,new);
		}else{
			new = ((unsigned long int)(si->len - errors[i])) << 56ul;
			new |= ((unsigned long int) ((unsigned long int) offset - (unsigned long int) positions[i]) << 1ul);
			new |= 1ul ;
			BinaryInsertionSortOne(si->match_units,LIST_STORE_SIZE ,new);
		}

	}
	return 1;
}
*/


// (c) 2017 Johannes Soeding & Martin Steinegger, Gnu Public License version 3
// Rotate left macro: left circular shift by numbits within 16 bits
#define RoL(val, numbits) (val << numbits) ^ (val >> (16 - numbits))


// Transform each letter x[i] to a fixed random number RAND[x[i]]
// to ensure instantaneous mixing into the 16 bits
// Do XOR with RAND[x[i]] and 5-bit rotate left for each i from 1 to k
uint16_t circ_hash(const uint_fast8_t* x, const uint_fast8_t length)
{
        const uint16_t RAND[21] = {0x4567, 0x23c6, 0x9869, 0x4873, 0xdc51, 0x5cff, 0x944a, 0x58ec,
                                   0x1f29, 0x7ccd, 0x58ba, 0xd7ab, 0x41f2, 0x1efb, 0xa9e3, 0xe146,
                                   0x007c, 0x62c2, 0x0854, 0x27f8, 0x231b};// 16 bit random numbers
        register uint16_t h = 0x0;
        h = h^ RAND[x[0]];// XOR h and ki
        for (register uint_fast8_t i = 1; i < length; ++i){
                h = RoL(h, 5);
                h ^= RAND[x[i]];// XOR h and ki
        }
        return h;
}

// Rolling hash variant for previous hash function:
// Computes hash value for next key x[0:length-1] from previous hash value
// hash( x[-1:length-2] ) and x_first = x[-1]

uint16_t circ_hash_next(const uint_fast8_t * x,const uint_fast8_t length,const uint_fast8_t x_first, uint16_t h)
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

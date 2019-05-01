

#include "misc.h"

int byg_detect(int* text,int n)
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



unsigned long int bpm(const unsigned char* t,const unsigned char* p,int n,int m)
{
	register unsigned long int i;//,c;
	unsigned long int diff;
	unsigned long int B[5];
	if(m > 64){
		m = 64;
	}

	unsigned long int k = m;
	//static int counter = 0;
	register unsigned long int VP,VN,D0,HN,HP,X;

	long int MASK = 0;
	//c = 0;

	diff = m;

	for(i = 0; i < 5;i++){
		B[i] = 0;
	}

	for(i = 0; i < m;i++){
		B[(int)(p[i] & 0x3)] |= (1ul << i);
	}

	//c = 0;
	VP = 0xFFFFFFFFFFFFFFFFul;
	VN = 0ul;
	m--;
	MASK = 1ul << m;
	for(i = 0; i < n;i++){
		X = (B[(int)(t[i] &0x3)  ] | VN);
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
	//fprintf(stderr,"%d	%d	%d	%d	",out1[0],out1[1],out1[2],out1[3]);
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

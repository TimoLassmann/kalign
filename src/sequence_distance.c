#include <xmmintrin.h>
#include "sequence_distance.h"

#include "alphabet.h"
#include "alignment.h"

#include "misc.h"
#define NODESIZE 16

struct kmer{
        uint64_t kmer;
        int32_t seq;
        int32_t pos;
        int32_t hash;
};

/* small hash implementation */
struct bignode{
        struct bignode *next;
        unsigned int pos[NODESIZE];
        unsigned int num;
};

struct bignode* big_insert_hash(struct bignode *n,const unsigned int pos);
void big_remove_nodes(struct bignode *n);
void big_print_nodes(struct bignode *n);

float protein_wu_distance_calculation(struct bignode* hash[],const uint8_t* seq,const int seqlen,const int diagonals,const float mode);

float dna_distance_calculation(struct bignode* hash[], const uint8_t * p,const int seqlen,int diagonals,float mode);;

int sort_by_kmer_then_seq(const void *a, const void *b);
int sort_by_hash(const void *a, const void *b);

float** kmer_distance(struct alignment* aln, int* seeds, int num_seeds, int kmer_len)
{
        struct alphabet* alphabet = NULL;
        struct kmer** kmer_list = NULL;
        int kmer_idx;
        int kmer_old_idx;
        float** dm = NULL;
        uint8_t* s = NULL;
        int_fast8_t* t = NULL;
        void *p;
        int64_t code;
        int64_t mask;
        uint16_t hash;
        int len;
        int numseq;
        //int numprofiles;
        int i,j,a;
        ASSERT(aln != NULL,"No alignment");

        numseq = aln->numseq;
        //numprofiles = aln->num_profiles;

        alphabet =  create_alphabet( redPROTEIN);
        t = alphabet->to_internal;
        MMALLOC(dm, sizeof(float*)* numseq);

        mask = (1LL << (kmer_len *4)) - 1LL;
        //fprintf(stdout,"MASK: %lx\n", mask);
        a = num_seeds / 8;
        if( num_seeds%8){
                a++;
        }
        a = a << 3;

        for(i = 0; i < numseq;i++){
                dm[i] = NULL;
                dm[i] = _mm_malloc(sizeof(float) * a,32);
                for(j = 0; j < a;j++){
                        dm[i][j] = 0.0f;
                }

        }
        /* allocate kmerlist */
        len = 0;
        for(i = 0; i<  num_seeds;i++){
                len += aln->sl[seeds[i]] - (kmer_len-1);
        }
        for(i = 0; i < aln->numseq;i++){
                len += aln->sl[i] - (kmer_len-1);
        }

        MMALLOC(kmer_list, sizeof(struct kmer*) * len);
        p = NULL;
        MMALLOC(p, sizeof(struct kmer) * len);

        for(i = 0;i <len;i++){
                kmer_list[i] = p;
                p+= sizeof(struct kmer);
        }

        p-= sizeof(struct kmer)* len;

        /* read in seeds  */
        LOG_MSG("len:%d",len);
        int g = 0;
        kmer_idx = 0;
        for(i = 0; i < aln->numseq;i++){
                kmer_old_idx = kmer_idx;
                code = 0;
                s = (uint8_t*)aln->seq[i];

                len = aln->sl[i];
                for(j = 0; j < kmer_len;j++){
                        code = code << 4L;
                        code |= t[s[j]];

                }
                hash = circ_hash(s, kmer_len);
                fprintf(stdout,"%d HASH!!\n",hash);
                kmer_list[kmer_idx]->kmer = code & mask;
                kmer_list[kmer_idx]->pos = 0;
                kmer_list[kmer_idx]->seq = i;
                kmer_list[kmer_idx]->hash = hash;
                kmer_idx++;

                for(j = 1; j < len- kmer_len;j++){
                        hash = circ_hash_next( s+j, kmer_len, s[j-1], hash);
                        fprintf(stdout,"%d %d",hash, s[j-1]);
                        code = code << 4L;
                        code |= t[s[j+kmer_len]];

                        //                fprintf(stdout,"pattern: %d %d %*lx (%c %d)\n", i,j, kmer_len, code & mask,(char) s[j], t[s[j]]);
                        kmer_list[kmer_idx]->kmer = code & mask;
                        kmer_list[kmer_idx]->pos = j;
                        kmer_list[kmer_idx]->seq = i;
                        kmer_list[kmer_idx]->hash = hash;
                        kmer_idx++;
                }
                LOG_MSG("Sorting for seq %d between %d and %d",i, kmer_old_idx, kmer_idx);
                qsort(kmer_list + kmer_old_idx, kmer_idx - kmer_old_idx, sizeof(struct kmer*), sort_by_hash);
                for(j = 0; j < 20;j++){
                        fprintf(stdout,"%d %d: %d %*lx %d %d SEED\n", i,j,kmer_list[j]->hash, kmer_len, kmer_list[j]->kmer, kmer_list[j]->seq,kmer_list[j]->pos);
                }

        }
        LOG_MSG("len:%d",g);
        exit(0);

        qsort(kmer_list, g, sizeof(struct kmer*), sort_by_kmer_then_seq);

        code = 0UL;
        int c,f;


        for(i = 0; i < g;i++){
                if(kmer_list[i]->kmer != code){
                        code = kmer_list[i]->kmer;
                        j = i;
                        while(kmer_list[j]->kmer == code){
                                j++;
                                if(j == g){
                                        break;
                                }
                        }
                        if(j-i > 1){
                                if(kmer_list[i]->seq & 0x80000000){


                                        /*fprintf(stdout,"BLOCK %d-%d\n",i,j);
                                        for(c = i; c < j;c++){
                                                if(kmer_list[c]->seq &0x80000000){
                                                        fprintf(stdout,"%d: %*lx %d %d SEED\n", c, kmer_len, kmer_list[c]->kmer, kmer_list[c]->seq & 0x7FFFFFFF,kmer_list[c]->pos);
                                                }else{
                                                        fprintf(stdout,"%d: %*lx %d %d\n", c, kmer_len, kmer_list[c]->kmer, kmer_list[c]->seq,kmer_list[c]->pos);
                                                }
                                        }


                                        fprintf(stdout,"PAIRS: %d-%d\n",i,j);*/
                                        for(c = i; c < j-1;c++){
                                                if(kmer_list[c]->seq & 0x80000000){
                                                        for(f = c+1;f < j;f++){
                                                                if(!(kmer_list[f]->seq & 0x80000000)){
                                                                        //fprintf(stdout,"%d %d\n", kmer_list[c]->seq &0x7FFFFFFF,kmer_list[f]->seq );

                                                                        dm[kmer_list[f]->seq] [kmer_list[c]->seq &0x7FFFFFFF] += 1.0f / (float) (abs( kmer_list[f]->pos - kmer_list[c]->pos)+1.0f);
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                        i =j-1;


                }

        }


        MFREE(p);
        MFREE(kmer_list);
        MFREE(alphabet);

        return dm;
ERROR:
        return NULL;
}



int sort_by_kmer_then_seq(const void *a, const void *b)
{
        struct kmer* const *one = a;
        struct kmer* const *two = b;

        if((*one)->kmer > (*two)->kmer){
                return -1;
        }else if((*one)->kmer == (*two)->kmer){
                if((*one)->seq < (*two)->seq){
                        return -1;
                }else{
                        return 1;
                }
        }else{
                return 1;
        }
}

int sort_by_hash(const void *a, const void *b)
{
        struct kmer* const *one = a;
        struct kmer* const *two = b;
        if((*one)->hash < (*two)->hash){
                return -1;
        }else{
                return 1;
        }
}



float** protein_wu_distance(struct alignment* aln, float zlevel, int nj, int* seeds, int num_anchors)
{
        struct bignode* hash[1024];
        float** dm = NULL;
        uint8_t*p =0;
        int i,j,a;
        unsigned int hv;
        int numseq;
        int numprofiles;
        //float min;
        float cutoff;

        ASSERT(aln != NULL,"No alignment");

        numseq = aln->numseq;
        numprofiles = aln->num_profiles;

        for (i = 0;i < 1024;i++){
                hash[i] = 0;
        }

        if(num_anchors){
                MMALLOC(dm, sizeof(float*)* numseq);

                a = num_anchors / 8;
                if( num_anchors%8){
                        a++;
                }
                a = a << 3;

                for(i = 0; i < numseq;i++){
                        dm[i] = NULL;
                        dm[i] = _mm_malloc(sizeof(float) * a,32);
                       for(j = 0; j < a;j++){
                                dm[i][j] = 0.0f;
                        }

                }

                for(i = 0; i < aln->numseq;i++){
                        p = aln->s[i];

                        for (j = aln->sl[i]-2;j--;){
                                //for(j = 0; j < si->sl[i]-2;j++){
                                //hv = (p[j+1] << 5) + p[j+2];
                                //hash[hv] = big_insert_hash(hash[hv],j);
                                hv = (p[j] << 5) + p[j+1];
                                hash[hv] = big_insert_hash(hash[hv],j);
                                hv = (p[j] << 5) + p[j+2];
                                hash[hv] = big_insert_hash(hash[hv],j);
                        }

                        for(j = 0;j < num_anchors;j++){

                                cutoff = zlevel;
                                //cutoff = param->zlevel;
                                p = aln->s[seeds[j]];
                                dm[i][j] = protein_wu_distance_calculation(hash,p,aln->sl[seeds[j]],aln->sl[seeds[j]]+aln->sl[i],cutoff);
                        }

                        for (j = 1024;j--;){
                                if (hash[j]){
                                        big_remove_nodes(hash[j]);
                                        hash[j] = 0;
                                }
                        }

                }

        }else{
                i = numseq;
                if (nj){
                        i = numprofiles;
                }

                RUNP(dm = galloc(dm,i,i,0.0f));
                //fprintf(stderr,"Distance Calculation:\n");
                //b = (numseq*(numseq-1))/2;
                a = 1;

                for (i = 0; i < numseq-1;i++){
                        p = aln->s[i];

                        for (j = aln->sl[i]-2;j--;){
                                //for(j = 0; j < si->sl[i]-2;j++){
                                //hv = (p[j+1] << 5) + p[j+2];
                                //hash[hv] = big_insert_hash(hash[hv],j);
                                hv = (p[j] << 5) + p[j+1];
                                hash[hv] = big_insert_hash(hash[hv],j);
                                hv = (p[j] << 5) + p[j+2];
                                hash[hv] = big_insert_hash(hash[hv],j);
                        }
                        for (j = i+1; j < numseq;j++){
                                //min =  (aln->sl[i] > aln->sl[j]) ? aln->sl[j] :aln->sl[i];
                                cutoff = zlevel;
                                //cutoff = param->zlevel;
                                p = aln->s[j];
                                dm[i][j] = protein_wu_distance_calculation(hash,p,aln->sl[j],aln->sl[j]+aln->sl[i],cutoff);
                                //fprintf(stderr,"%d-%d:%f\n",i,j,dm[i][j]);
                                //exit(0);
                                //dm[i][j] /= min;
                                //dm[i][j] /= (si->sl[i] > si->sl[j]) ? si->sl[j] :si->sl[i];
                                dm[j][i] = dm[i][j];
                                //fprintf(stderr,"\r%8.0f percent done",(float)a /(float)b * 100);
                                a++;
                        }

                        for (j = 1024;j--;){
                                if (hash[j]){
                                        big_remove_nodes(hash[j]);
                                        hash[j] = 0;
                                }
                        }
                }
        }
        return dm;
ERROR:
        return NULL;
}


float protein_wu_distance_calculation(struct bignode* hash[],const uint8_t* seq,const int seqlen,const int diagonals,const float mode)
{

        struct bignode* node_p;
        unsigned int* d = NULL;
        unsigned int* tmp = NULL;
        float out = 0.0;
        register int i,j;
        register int c;
        register int num;
        register unsigned int hv;

        d = malloc(sizeof(unsigned int)*diagonals);
        //for (i = diagonals;i--;){
        for (i = 0;i < diagonals;i++){
                d[i] = 0;
        }
        for (i = seqlen-2;i--;){
                //for(i = 0; i < seqlen-2;i++){
                /*hv = (seq[i+1] << 5) + seq[i+2];

                  node_p = hash[hv];
                  while(node_p){
                  tmp = node_p->pos;
                  for(j = 0;j < node_p->num;j++){
                  d[tmp[j]]++;
                  }
                  node_p = node_p->next;
                  }*/
                hv = (seq[i] << 5) + seq[i+1];
                //printf("3:%d\n",hv);
                node_p = hash[hv];
                while(node_p){
                        tmp = node_p->pos;
                        num = node_p->num;
                        for(j = 0;j < num;j++){
                                c = tmp[j];
                                d[c]++;
                                c++;
                                d[c]++;
                        }
                        node_p = node_p->next;
                }
                /*hv = (seq[i] << 5) + seq[i+2];

                node_p = hash[hv];

                while(node_p){
                        tmp = node_p->pos;
                        num = node_p->num;
                        for(j = 0;j < num;j++){
                                c = tmp[j];
                                d[c]++;
                        }
                        node_p = node_p->next;
                        }*/
                d++;



        }
        //exit(0);
        d -= (seqlen-2);
        //unsigned int max = 0.0;
        for (i = diagonals;i--;){
                //      if(d[i] > max){
                //      max = d[i];
                //}
                //d[i] /= minlen;

                //fprintf(stderr,"%d ",d[i]);
                if(d[i] > mode){
                        out += d[i];
                        //	printf("%f	%d\n",d[i]/ minlen,d[i]);
                }
        }
        free(d);
        //return out;
        return  out;
}



float** dna_distance(struct alignment* aln, float zlevel, int nj)
{
        struct bignode* hash[1024];
        float** dm = NULL;
        uint8_t* p = NULL;
        int i,j,a;
        unsigned int hv;
        int numseq;
        int numprofiles;

        ASSERT(aln != NULL,"No alignment");

        numseq = aln->numseq;
        numprofiles = aln->num_profiles;

        //fprintf(stderr,"Distance Calculation:\n");


        for (i = 0;i < 1024;i++){
                hash[i] = 0;
        }

        i = numseq;
        if (nj){
                i = numprofiles;
        }

        RUNP(dm = galloc(dm,i,i,0.0f));


        //b = (numseq*(numseq-1))/2;
        a = 1;

        for (i = 0; i < numseq-1;i++){
                p = aln->s[i];
                for (j = aln->sl[i]-5;j--;){
                        hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+2]&3)<<4)  + ((p[j+3]&3)<<2) + (p[j+4]&3);//ABCDE
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+2]&3)<<4)  + ((p[j+3]&3)<<2) + (p[j+5]&3);//ABCDF
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+2]&3)<<4)  + ((p[j+4]&3)<<2) + (p[j+5]&3);//ABCEF
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+3]&3)<<4)  + ((p[j+4]&3)<<2) + (p[j+5]&3);//ABDEF
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = ((p[j]&3)<<8) + ((p[j+2]&3)<<6) + ((p[j+3]&3)<<4) + ((p[j+4]&3)<<2) + (p[j+5]&3);//ACDEF
                        hash[hv] = big_insert_hash(hash[hv],j);
                }
                for (j = i+1; j < numseq;j++){


                        //min =  (aln->sl[i] > aln->sl[j]) ?aln->sl[j] :aln->sl[i];
                        dm[i][j] = dna_distance_calculation(hash,aln->s[j],aln->sl[j],aln->sl[j]+aln->sl[i],zlevel);
                        dm[i][j] /= (aln->sl[i] > aln->sl[j]) ?aln->sl[j] :aln->sl[i];
                        dm[j][i] = dm[i][j];
                        //fprintf(stderr,"\r%8.0f percent done",(float)a /(float)b * 100);
                        a++;
                }

                for (j = 1024;j--;){
                        if (hash[j]){
                                big_remove_nodes(hash[j]);
                                hash[j] = 0;
                        }
                }
        }
        return dm;
ERROR:
        return NULL;
}

float dna_distance_calculation(struct bignode* hash[],const uint8_t * p,const int seqlen,int diagonals,float mode)
{

        struct bignode* node_p;
        float out = 0.0;
        unsigned int* tmp = NULL;
        unsigned int* d = NULL;
        int i,j;
        unsigned int hv;

        d = malloc(sizeof(int)*diagonals);
        for (i = 0;i < diagonals;i++){
                d[i] = 0;
        }
        for (i = seqlen-5;i--;){

                hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+2]&3)<<4)  + ((p[i+3]&3)<<2) + (p[i+4]&3);//ABCDE
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }


                hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+2]&3)<<4)  + ((p[i+3]&3)<<2) + (p[i+5]&3);//ABCDF
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }
                hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+2]&3)<<4)  + ((p[i+4]&3)<<2) + (p[i+5]&3);//ABCEF
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }
                hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+3]&3)<<4)  + ((p[i+4]&3)<<2) + (p[i+5]&3);//ABDEF
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }
                hv = ((p[i]&3)<<8) + ((p[i+2]&3)<<6) + ((p[i+3]&3)<<4) + ((p[i+4]&3)<<2) + (p[i+5]&3);//ACDEF
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }

                d++;
        }
        //exit(0);
        d -= (seqlen-5);

        for (i = diagonals;i--;){
                //d[i] /= minlen;

                //printf("%d ",d[i]);

                if(d[i] > mode){
                        //fprintf(stderr,"%f	%d\n",d[i]/ minlen,d[i]);
                        out += d[i];
                }
        }
        free(d);
        return out;
}


struct bignode* big_insert_hash(struct bignode *n,const unsigned int pos)
{
        struct bignode* p = NULL;
        if(n){
                if(n->num < NODESIZE){
                        n->pos[n->num] = pos;
                        n->num++;
                        return n;
                }else{
                        MMALLOC(p, sizeof(struct bignode));
                        p->pos[0] = pos;
                        p->num = 1;
                        p->next = n;
                }
        }else{
                MMALLOC(p, sizeof(struct bignode));
                p->pos[0] = pos;
                p->num = 1;
                p->next = n;
        }
        return p;
ERROR:
        return NULL;
}

void big_remove_nodes(struct bignode *n)
{
        struct bignode* p = NULL;
        while(n){
                p = n;
                n = n->next;
                MFREE(p);
        }
}

void big_print_nodes(struct bignode *n)
{
        int i;
        while(n){
                for (i = 0; i < n->num;i++){
                        fprintf(stderr,"%d ",n->pos[i]);
                }
                n = n->next;
        }
}

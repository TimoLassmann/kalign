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

#ifdef HAVE_AVX2
#include <xmmintrin.h>
#endif

#include <mm_malloc.h>
#include "sequence_distance.h"

#include "alphabet.h"
/* #include "alignment.h" */
#include "align_io.h"

#include "misc.h"
#include "bpm.h"

#define NODESIZE 16

/* small hash implementation */
struct bignode{
        struct bignode *next;
        unsigned int pos[NODESIZE];
        unsigned int num;
};

struct bignode* big_insert_hash(struct bignode *n,const unsigned int pos);
void big_remove_nodes(struct bignode *n);
void big_print_nodes(struct bignode *n);


float dna_distance_calculation(struct bignode* hash[],const uint8_t * p,const int seqlen,int diagonals,float mode);
float protein_wu_distance_calculation(struct bignode* hash[],const uint8_t* seq,const int seqlen,const int diagonals,const float mode);

float** d_estimation(struct msa* msa, int* samples, int num_samples,int pair)
{
        float** dm = NULL;
        uint8_t* seq_a;
        uint8_t* seq_b;

        float dist;

        int len_a;
        int len_b;

        int i,j;
#if HAVE_AVX2
        set_broadcast_mask();
#endif

        if(pair){

                RUN(galloc(&dm,num_samples,num_samples));
                for(i = 0; i < num_samples;i++){

                        seq_a = msa->sequences[samples[i]]->s;// aln->s[samples[i]];
                        len_a = msa->sequences[samples[i]]->len;//aln->sl[samples[i]];
                        for(j = 0;j < num_samples;j++){
                                //fprintf(stdout, "Working on %d %d\n", i,j);

                                seq_b = msa->sequences[samples[j]]->s; //aln->s[ samples[j]];
                                len_b = msa->sequences[samples[j]]->len;//aln->sl[selection[j]];
                                /*dm[i][j] = MACRO_MIN(len_a, len_b) - MACRO_MIN(
                                  bpm_256(seq_a, seq_b, len_a, len_b),
                                  bpm_256(seq_b, seq_a, len_b, len_a)
                                  );*/
                                dist = calc_distance(seq_a, seq_b, len_a, len_b,msa->L);
                                //dist = dist / (float) MACRO_MIN(len_a, len_b);
                                dm[i][j] = dist;// + (float)( i * num_samples + j) / (float) ( num_samples * num_samples);
                                dm[j][i] = dm[i][j];
                                //fprintf(stdout,"%f ", dm[i][j]);
                        }

                        //fprintf(stdout,"\n");
                }
        }else{
                int a;
                int numseq = msa->numseq;
                MMALLOC(dm, sizeof(float*)* numseq);
                //fprintf(stdout,"MASK: %lx\n", mask);
                a = num_samples / 8;
                if( num_samples%8){
                        a++;
                }
                a = a << 3;

                for(i = 0; i < numseq;i++){
                        dm[i] = NULL;

                        dm[i] = _mm_malloc(sizeof(float) * a,32);
                        for(j = 0; j < a;j++){
                                dm[i][j] = 0.0F;
                        }
                }


                struct msa_seq** s = msa->sequences;
#ifdef HAVE_OPENMP
#pragma omp parallel for shared(dm, s) private(i, j) collapse(2) schedule(static)
#endif
                for(i = 0; i < numseq;i++){
                        for(j = 0;j < num_samples;j++){
                                uint8_t* s1;
                                uint8_t* s2;
                                int l1;
                                int l2;
                                s1 = s[i]->s;
                                l1 = s[i]->len;
                                s2 = s[samples[j]]->s;
                                l2 = s[samples[j]]->len;
                                dm[i][j] = calc_distance(s1,
                                                         s2,
                                                         l1,
                                                         l2,
                                                         msa->L);

                                //dm[i][j] += (float)MACRO_MIN(l1, l2) / (float)MACRO_MAX(l1, l2);
                                //dm[i][j] = dm[i][j] / (float) MACRO_MIN(l1, l2);
                                //dm[i][j] = dist;
                        }
                }
                /* for(i = 0; i < numseq;i++){ */
                /*         seq_a = msa->sequences[i]->s;// aln->s[i]; */
                /*         len_a = msa->sequences[i]->len;//  aln->sl[i]; */
                /*         for(j = 0;j < num_samples;j++){ */
                /*                 seq_b = msa->sequences[samples[j]]->s;// aln->s[ seeds[j]]; */
                /*                 len_b = msa->sequences[samples[j]]->len;// aln->sl[seeds[j]]; */

                /*                 dist = calc_distance(seq_a, seq_b, len_a, len_b,msa->L); */
                /*                 dm[i][j] = dist; */
                /*         } */
                /* } */
        }


        return dm;
ERROR:
        return NULL;
}

float calc_distance(uint8_t* seq_a, uint8_t* seq_b, int len_a,int len_b, int L)
{
#ifdef HAVE_AVX2
        uint8_t dist;
        if(len_a > len_b){
                dist = bpm_256(seq_a, seq_b, len_a, len_b);
        }else{
                dist = bpm_256(seq_b, seq_a, len_b, len_a);
        }
        return (float)dist;
#else
        struct bignode* hash[1024];
        int i;
        float dist;
        unsigned int hv;
        for (i = 0;i < 1024;i++){
                hash[i] = 0;
        }
        /* Protein sequence  */
        if( L > ALPHA_defDNA){

                for (i = len_a-2;i--;){
                        hv = (seq_a[i] << 5) + seq_a[i+1];
                        hash[hv] = big_insert_hash(hash[hv],i);
                        hv = (seq_a[i] << 5) + seq_a[i+2];
                        hash[hv] = big_insert_hash(hash[hv],i);
                }

                dist = protein_wu_distance_calculation(hash,seq_b,len_b,len_a+len_b,58.9);
        }else{

                for (i = len_a-5;i--;){
                        hv = ((seq_a[i]&3)<<8) + ((seq_a[i+1]&3)<<6) + ((seq_a[i+2]&3)<<4)  + ((seq_a[i+3]&3)<<2) + (seq_a[i+4]&3);//ABCDE
                        hash[hv] = big_insert_hash(hash[hv],i);
                        hv = ((seq_a[i]&3)<<8) + ((seq_a[i+1]&3)<<6) + ((seq_a[i+2]&3)<<4)  + ((seq_a[i+3]&3)<<2) + (seq_a[i+5]&3);//ABCDF
                        hash[hv] = big_insert_hash(hash[hv],i);
                        hv = ((seq_a[i]&3)<<8) + ((seq_a[i+1]&3)<<6) + ((seq_a[i+2]&3)<<4)  + ((seq_a[i+4]&3)<<2) + (seq_a[i+5]&3);//ABCEF
                        hash[hv] = big_insert_hash(hash[hv],i);
                        hv = ((seq_a[i]&3)<<8) + ((seq_a[i+1]&3)<<6) + ((seq_a[i+3]&3)<<4)  + ((seq_a[i+4]&3)<<2) + (seq_a[i+5]&3);//ABDEF
                        hash[hv] = big_insert_hash(hash[hv],i);
                        hv = ((seq_a[i]&3)<<8) + ((seq_a[i+2]&3)<<6) + ((seq_a[i+3]&3)<<4) + ((seq_a[i+4]&3)<<2) + (seq_a[i+5]&3);//ACDEF
                        hash[hv] = big_insert_hash(hash[hv],i);
                }
                dist = dna_distance_calculation(hash,seq_b,len_b,len_a+len_b, 61.08);
        }


        for (i = 1024;i--;){
                if (hash[i]){
                        big_remove_nodes(hash[i]);
                        hash[i] = 0;
                }
        }
        return dist;
#endif

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
                hv = (seq[i] << 5) + seq[i+2];

                node_p = hash[hv];

                while(node_p){
                        tmp = node_p->pos;
                        num = node_p->num;
                        for(j = 0;j < num;j++){
                                c = tmp[j];
                                d[c]++;
                        }
                        node_p = node_p->next;
                        }
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
                                for(j = 0;j < (int) node_p->num;j++){
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
                                for(j = 0;j < (int)node_p->num;j++){
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
                                for(j = 0;j < (int)node_p->num;j++){
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
                                for(j = 0;j < (int)node_p->num;j++){
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
                                for(j = 0;j < (int)node_p->num;j++){
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
                for (i = 0; i < (int)n->num;i++){
                        fprintf(stderr,"%d ",n->pos[i]);
                }
                n = n->next;
        }
}

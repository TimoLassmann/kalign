#include <stdint.h>

#include "tldevel.h"

#include "alignment_parameters.h"

#include "aln_struct.h"
#define ALN_SEQSEQ
#include "aln_seqseq.h"
#define MAX(a, b) (a > b ? a : b)
#define MAX3(a,b,c) MAX(MAX(a,b),c)

int aln_seqseq_foward(struct aln_mem* m)
{
        struct states* s = m->f;
        const uint8_t* seq1 = m->seq1;
        const uint8_t* seq2 = m->seq2;
        float *subp = 0;
        const int starta = m->starta;
        const int enda = m->enda;
        const int startb =m->startb;
        const int endb = m->endb;
        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;
        register float xa = 0;
        register float xga = 0;
        register int i = 0;
        register int j = 0;

        const float gpo = m->ap->gpo;
        const float gpe = m->ap->gpe;
        const float tgpe = m->ap->tgpe;
        float** subm = m->ap->subm;

        s[startb].a = s[0].a;
        s[startb].ga = s[0].ga;
        s[startb].gb = s[0].gb;
        if(startb){
                for (j = startb+1; j < endb;j++){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j-1].ga - gpe,s[j-1].a-gpo);
                        s[j].gb = -FLT_MAX;
                }
        }else{
                for (j = startb+1; j < endb;j++){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j-1].ga,s[j-1].a)-tgpe;
                        s[j].gb = -FLT_MAX;
                }
        }
        s[endb].a = -FLT_MAX;
        s[endb].ga = -FLT_MAX;
        s[endb].gb = -FLT_MAX;

        seq2--;
        for (i = starta;i < enda;i++){
                subp = subm[seq1[i]];

                pa = s[startb].a;
                pga = s[startb].ga;
                pgb = s[startb].gb;
                s[startb].a = -FLT_MAX;
                s[startb].ga = -FLT_MAX;

                xa = s[startb].a;
                xga = s[startb].ga;

                if(startb){
                        s[startb].gb = MAX(pgb - gpe,pa - gpo);
                }else{
                        s[startb].gb = MAX(pgb,pa) - tgpe;
                }
                for (j = startb+1; j < endb;j++){
                        ca = s[j].a;
                        pa = MAX3(pa,pga-gpo,pgb-gpo);
                        pa += subp[seq2[j]];

                        s[j].a = pa;

                        pga = s[j].ga;
                        //s[j].ga = MAX(s[j-1].ga-gpe,s[j-1].a-gpo);
                        s[j].ga = MAX(xga-gpe,xa-gpo);

                        pgb = s[j].gb;
                        s[j].gb = MAX(pgb-gpe ,ca-gpo);

                        pa = ca;

                        xa = s[j].a;
                        xga = s[j].ga;

                }
                ca = s[j].a;
                pa = MAX3(pa,pga-gpo,pgb-gpo);
                pa += subp[seq2[j]];

                s[j].a = pa;

                s[j].ga = -FLT_MAX;//MAX(s[j-1].ga-gpe,s[j-1].a-gpo);
                if (endb != m->len_b){
                        s[j].gb = MAX(s[j].gb-gpe ,ca-gpo);
                }else{
                        s[j].gb = MAX(s[j].gb,ca)-tgpe;
                }

        }
        return OK;
}

int aln_seqseq_backward(struct aln_mem* m)
{

        struct states* s = m->b;
        const uint8_t* seq1 = m->seq1;
        const uint8_t* seq2 = m->seq2;
        float *subp = NULL;
        const int starta = m->starta_2;
        const int enda = m->enda_2;
        const int startb =m->startb;
        const int endb = m->endb;

        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;

        register float xa = 0;
        register float xga = 0;

        register int i = 0;
        register int j = 0;

        const float gpo = m->ap->gpo;
        const float gpe = m->ap->gpe;
        const float tgpe = m->ap->tgpe;
        float** subm = m->ap->subm;


        s[endb].a = s[0].a ;
        s[endb].ga = s[0].ga;
        s[endb].gb = s[0].gb;


        //init of first row;

        //j = endb-startb;
        if(endb != m->len_b){
                for(j = endb-1;j > startb;j--){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);
                        s[j].gb = -FLT_MAX;
                }
        }else{
                for(j = endb-1;j > startb;j--){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j+1].ga,s[j+1].a)-tgpe;
                        s[j].gb = -FLT_MAX;
                }
        }


        s[startb].a = -FLT_MAX;
        s[startb].ga = -FLT_MAX;
        s[startb].gb = -FLT_MAX;

        i = enda-starta;
        seq1+= starta;
        while(i--){
                subp = subm[seq1[i]];
                pa = s[endb].a;
                pga = s[endb].ga;
                pgb = s[endb].gb;
                s[endb].a = -FLT_MAX;
                s[endb].ga = -FLT_MAX;

                xa = s[endb].a;
                xga = s[endb].ga;

                if(endb != m->len_b){
                        s[endb].gb = MAX(pgb-gpe,pa-gpo);
                }else{
                        s[endb].gb = MAX(pgb,pa)-tgpe;
                }

                for(j = endb-1;j > startb;j--){

                        ca = s[j].a;

                        pa = MAX3(pa,pga - gpo,pgb-gpo);

                        pa += subp[seq2[j]];

                        s[j].a = pa;

                        pga = s[j].ga;

                        //s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);

                        s[j].ga = MAX(xga-gpe,xa-gpo);

                        pgb = s[j].gb;
                        s[j].gb = MAX(pgb-gpe,ca-gpo);

                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;
                }
                ca = s[j].a;

                pa = MAX3(pa,pga - gpo,pgb-gpo);

                pa += subp[seq2[j]];

                s[j].a = pa;

                s[j].ga = -FLT_MAX;//MAX(s[j+1].ga-gpe,s[j+1].a-gpo);

                if(startb){
                        s[j].gb = MAX(s[j].gb-gpe,ca-gpo);
                }else{
                        s[j].gb = MAX(s[j].gb,ca)-tgpe;
                }
        }
        return OK;
}


int aln_seqseq_meetup(struct aln_mem* m,int old_cor[],int* meet,int* t,float* score)
{
        struct states* f = m->f;
        struct states* b = m->b;


        const float gpo = m->ap->gpo;
        const float gpe = m->ap->gpe;
        const float tgpe = m->ap->tgpe;
        int i;
        int c;
        int transition = -1;

        //code:
        // a -> a = 1
        // a -> ga = 2
        // a -> gb = 3
        // ga ->ga = 4
        // ga -> a = 5
        //gb->gb = 6;
        //gb->a = 7;

        //int max = -FLT_MAX;
        float max = -FLT_MAX;
        //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
        float middle =  (float)(old_cor[3] - old_cor[2])/2.0F + (float)old_cor[2];
        float sub = 0.0F;

        //i = hm->startb;

        c = -1;
        //for(i = hm->startb; i < hm->endb;i++){
        for(i = old_cor[2]; i < old_cor[3];i++){
                sub = fabsf(middle - (float)i);
                sub /= 1000.0F;
                //	fprintf(stderr,"%d-%d	%f\n",hm->startb,hm->endb,sub);
                if(f[i].a+b[i].a-sub > max){
                        max = f[i].a+b[i].a-sub;
                        //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
                        transition = 1;
                        c = i;
                }
                if(f[i].a+b[i].ga-gpo-sub > max){
                        max = f[i].a+b[i].ga-gpo-sub;
                        //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
                        transition = 2;
                        c = i;
                }
                if(f[i].a+b[i].gb -gpo-sub > max){
                        max = f[i].a+b[i].gb - gpo-sub;
                        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                        transition = 3;
                        c = i;
                }
                if(f[i].ga+b[i].a - gpo-sub > max){
                        max = f[i].ga+b[i].a - gpo-sub;
                        //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
                        transition = 5;
                        c = i;
                }


                if(m->startb == 0){
                        if(f[i].gb+b[i].gb - tgpe-sub > max){
                                max = f[i].gb+b[i].gb -tgpe-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }else{
                        if(f[i].gb+b[i].gb - gpe -sub> max){
                                max = f[i].gb+b[i].gb - gpe-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }
                if(f[i].gb+b[i].a - gpo-sub > max){
                        max = f[i].gb+b[i].a - gpo-sub;
                        //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
                        transition = 7;
                        c = i;
                }
        }
        //i = hm->endb;
        i = old_cor[3];
        sub = fabsf(middle - (float)i);
        sub /= 1000.0F;

        if(f[i].a+b[i].gb-gpo-sub > max){
                max = f[i].a+b[i].gb - gpo-sub;
                //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                transition = 3;
                c = i;
        }
        if(m->endb == m->len_b){
                if(f[i].gb+b[i].gb -tgpe-sub > max){
                        max = f[i].gb+b[i].gb - tgpe-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }else{
                if(f[i].gb+b[i].gb - gpe-sub > max){
                        max = f[i].gb+b[i].gb - gpe-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }
        *meet = c;
        *t = transition;
        *score = max;
        return OK;
}

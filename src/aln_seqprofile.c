#include <stdint.h>

#include "tldevel.h"

#include "alignment_parameters.h"

#include "aln_struct.h"
#define ALN_SEQPROF_IMPORT
#include "aln_seqprofile.h"


#define MAX(a, b) (a > b ? a : b)
#define MAX3(a,b,c) MAX(MAX(a,b),c)


int aln_seqprofile_foward(struct aln_mem* m)
{
        struct states* s = m->f;
        const float* prof1 = m->prof1;
        const uint8_t* seq2 = m->seq2;
        const int sip = m->sip;

        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;

        register float xa = 0;
        register float xga = 0;

        register int i = 0;
        register int j = 0;

        const float open = m->ap->gpo  * (float)sip;
        const float ext  = m->ap->gpe  * (float)sip;
        const float text = m->ap->tgpe * (float)sip;


        prof1 += (m->starta)<< 6;
        s[m->startb].a = s[0].a;
        s[m->startb].ga = s[0].ga;
        s[m->startb].gb = s[0].gb;
        if(m->startb){
                for (j = m->startb +1; j < m->endb;j++){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
                        s[j].gb = -FLT_MAX;
                }
        }else{
                for (j = m->startb+1; j < m->endb;j++){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;
                        s[j].gb = -FLT_MAX;
                }

        }


        s[m->endb].a = -FLT_MAX;
        s[m->endb].ga = -FLT_MAX;
        s[m->endb].gb = -FLT_MAX;
        seq2--;

        for (i = m->starta;i < m->enda;i++){
                prof1 += 64;

                pa = s[m->startb].a;
                pga = s[m->startb].ga;
                pgb = s[m->startb].gb;
                s[m->startb].a = -FLT_MAX;
                s[m->startb].ga = -FLT_MAX;

                xa = s[m->startb].a;
                xga = s[m->startb].ga;


                if(m->startb){
                        s[m->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
                }else{
                        s[m->startb].gb = MAX(pgb,pa)+prof1[29];
                }
                for (j = m->startb+1; j < m->endb;j++){
                        ca = s[j].a;

                        pa = MAX3(pa,pga -open,pgb + prof1[-37]);

                        pa += prof1[32 + seq2[j]];


                        s[j].a = pa;

                        pga = s[j].ga;

                        //s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
                        s[j].ga = MAX(xga-ext,xa-open);


                        pgb = s[j].gb;

                        s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);

                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;

                }
                ca = s[j].a;

                pa = MAX3(pa,pga -open,pgb + prof1[-37]);

                pa += prof1[32 + seq2[j]];


                s[j].a = pa;

                s[j].ga = -FLT_MAX;//MAX(s[j-1].ga-ext,s[j-1].a-open);

                if (m->endb != m->len_b){
                        s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);
                }else{
                        s[j].gb = MAX(s[j].gb,ca)+ prof1[29];
                }

        }
        //prof1 -= m->enda << 6;
        return OK;
}

int aln_seqprofile_backward(struct aln_mem* m)
{
        struct states* s = m->b;
        const float* prof1 = m->prof1;
        const uint8_t* seq2 = m->seq2;
        const int sip = m->sip;
        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;

        register float xa = 0;
        register float xga = 0;

        register int i = 0;
        register int j = 0;

        const float open = m->ap->gpo  * (float)sip;
        const float ext  = m->ap->gpe  * (float)sip;
        const float text = m->ap->tgpe * (float)sip;


        prof1 += (m->enda_2 +1) << 6;

        s[m->endb].a = s[0].a;
        s[m->endb].ga = s[0].ga;
        s[m->endb].gb = s[0].gb;

        if(m->endb != m->len_b){
                for(j = m->endb-1;j > m->startb ;j--){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
                        s[j].gb = -FLT_MAX;
                }
        }else{
                for(j = m->endb-1;j > m->startb;j--){
                        s[j].a = -FLT_MAX;
                        s[j].ga = MAX(s[j+1].ga,s[j+1].a)-text;
                        s[j].gb = -FLT_MAX;
                }
        }

        s[m->startb].a = -FLT_MAX;
        s[m->startb].ga = -FLT_MAX;
        s[m->startb].gb = -FLT_MAX;

        i = m->enda_2 -m->starta_2;
        while(i--){
                prof1 -= 64;
                pa = s[m->endb].a;
                pga = s[m->endb].ga;
                pgb = s[m->endb].gb;
                s[m->endb].a = -FLT_MAX;
                s[m->endb].ga = -FLT_MAX;

                xa = s[m->endb].a;
                xga = s[m->endb].ga;


                if(m->endb != m->len_b){
                        s[m->endb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
                }else{
                        s[m->endb].gb = MAX(pgb,pa) +prof1[29];
                }

                for(j = m->endb-1;j > m->startb;j--){
                        ca = s[j].a;

                        pa = MAX3(pa,pga - open,pgb +prof1[91]);
                        pa += prof1[32 + seq2[j]];

                        s[j].a = pa;

                        pga = s[j].ga;

                        //s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
                        s[j].ga = MAX(xga-ext,xa-open);

                        pgb = s[j].gb;

                        s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);

                        pa = ca;
                        xa = s[j].a;
                        xga = s[j].ga;


                }
                ca = s[j].a;

                pa = MAX3(pa,pga - open,pgb +prof1[91]);
                pa += prof1[32 + seq2[j]];

                s[j].a = pa;


                s[j].ga = -FLT_MAX;//MAX(s[j+1].ga-ext,s[j+1].a-open);
                if(m->startb){
                        s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);
                }else{
                        s[j].gb = MAX(s[j].gb,ca)+prof1[29];
                }

        }
        return OK;
}

int aln_seqprofile_meetup(struct aln_mem* m,int old_cor[],int* meet,int* t,float* score)
{
        struct states* f = m->f;
        struct states* b = m->b;
        const float* prof1 = m->prof1;
        int sip = m->sip;
        int i;
        int c;
        int transition = -1;

        const float open = m->ap->gpo * (float)sip;


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
        //float middle =  (m->endb - m->startb)/2 + m->startb;
        float middle = (float)(old_cor[3] - old_cor[2])/2.0F + (float)old_cor[2];
        float sub = 0.0F;


        prof1+= ((old_cor[4]+1)<<6);

        //i = m->startb;

        c = -1;
        //for(i = m->startb; i < m->endb;i++){
        for(i = old_cor[2]; i < old_cor[3];i++){
                sub = fabsf(middle -i);
                sub /= 1000.0F;
                if(f[i].a+b[i].a-sub> max){
                        max = f[i].a+b[i].a-sub;
                        //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
                        transition = 1;
                        c = i;
                }
                if(f[i].a+b[i].ga-open-sub > max){
                        max = f[i].a+b[i].ga-open-sub;
                        //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
                        transition = 2;
                        c = i;
                }
                if(f[i].a+b[i].gb+prof1[27]-sub > max){
                        max = f[i].a+b[i].gb+prof1[27]-sub;
                        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                        transition = 3;
                        c = i;
                }
                if(f[i].ga+b[i].a-open-sub > max){
                        max = f[i].ga+b[i].a-open-sub;
                        //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
                        transition = 5;
                        c = i;
                }


                if(m->startb == 0){
                        if(f[i].gb+b[i].gb+prof1[29]-sub > max){
                                max = f[i].gb+b[i].gb+prof1[29]-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }else{
                        if(f[i].gb+b[i].gb+prof1[28]-sub > max){
                                max = f[i].gb+b[i].gb+prof1[28]-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }
                if(f[i].gb+b[i].a+prof1[-37]-sub > max){
                        max = f[i].gb+b[i].a+prof1[-37]-sub;
                        //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
                        transition = 7;
                        c = i;
                }
        }
        //i = m->endb;
        i = old_cor[3];

        sub = fabsf(middle -i);
        sub /= 1000.0F;
        if(f[i].a+b[i].gb+prof1[27]-sub > max){
                max = f[i].a+b[i].gb+prof1[27]-sub;
                //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                transition = 3;
                c = i;
        }
        if(m->endb == m->len_b){
                if(f[i].gb+b[i].gb+prof1[29]-sub > max){
                        max = f[i].gb+b[i].gb+prof1[29]-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }else{
                if(f[i].gb+b[i].gb+prof1[28]-sub > max){
                        max = f[i].gb+b[i].gb+prof1[28]-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }


        //prof1-= ((ap->old_cor[4]+1)<<6);

        *meet = c;
        *t = transition;
        *score = max;
        return OK;
}

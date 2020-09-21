#include <stdint.h>

#include "tldevel.h"

#include "alignment_parameters.h"

#include "aln_struct.h"
#define ALN_SEQSEQ
#include "aln_seqseq.h"

static int aln_seqseq_foward(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct aln_mem* m);
static int aln_seqseq_backward(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct aln_mem* m);

static int aln_seqseq_join(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct aln_mem* m,int* path,float input_states[],int old_cor[]);

int aln_seqseq_runner(const struct aln_param* ap, const uint8_t* seq1,const uint8_t* seq2,struct aln_mem* m, int* path)
{
        int mid = ((m->enda - m->starta) / 2)+ m->starta;
        float input_states[6] = {
                m->f[0].a,
                m->f[0].ga,
                m->f[0].gb,
                m->b[0].a,
                m->b[0].ga,
                m->b[0].gb};
        int old_cor[5] = {
                m->starta,
                m->enda,
                m->startb,
                m->endb,
                mid};

        if(m->starta >= m->enda){
                return OK;//hirsch_path;
        }
        if(m->startb >= m->endb){
                return OK;///hirsch_path;
        }


        m->enda = mid;

        //fprintf(stderr,"Forward:%d-%d	%d-%d\n",m->starta,m->enda,m->startb,m->endb);
        aln_seqseq_foward(ap,seq1,seq2,m);

        m->starta = mid;
        m->enda = old_cor[1];
        //fprintf(stderr,"Backward:%d-%d	%d-%d\n",m->starta,m->enda,m->startb,m->endb);
        aln_seqseq_backward(ap,seq1,seq2,m);


        aln_seqseq_join(ap,seq1,seq2,m,path,input_states,old_cor);
        return  OK;
}

int aln_seqseq_foward(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct aln_mem* m)
{
        struct states* s = m->f;
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

        float gpo,gpe,tgpe;
        float** subm = NULL;
        gpo = ap->gpo;
        gpe = ap->gpe;
        tgpe = ap->tgpe;
        subm = ap->subm;


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

int aln_seqseq_backward(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct aln_mem* m)
{

        struct states* s = m->b;
        float** subm = NULL;
        float *subp = NULL;
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
        float gpo,gpe,tgpe;

        gpo = ap->gpo;
        gpe = ap->gpe;
        tgpe = ap->tgpe;
        subm = ap->subm;

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




int aln_seqseq_join(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct aln_mem* m,int* path,float input_states[],int old_cor[])
{
        struct states* f = m->f;
        struct states* b = m->b;
        int i,c;
        int transition = -1;

        float gpo,gpe,tgpe;

        gpo = ap->gpo;
        gpe = ap->gpe;
        tgpe = ap->tgpe;


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
        float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
        float sub = 0.0;

        //i = hm->startb;
        i = old_cor[2];
        c = -1;
        //for(i = hm->startb; i < hm->endb;i++){
        for(i = old_cor[2]; i < old_cor[3];i++){

                sub = fabsf(middle - i);
                sub /= 1000;
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
        sub = fabsf(middle -i);
        sub /= 1000;

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


        //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
        //LOG_MSG("MAX: %f",max);
        //j = hirsch_path[0];
        switch(transition){
        case 1: //a -> a = 1

                path[old_cor[4]] = c;
                path[old_cor[4]+1] = c+1;

                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                m->f[0].a = input_states[0];
                m->f[0].ga = input_states[1];
                m->f[0].gb = input_states[2];
                m->b[0].a = 0.0;
                m->b[0].ga = -FLT_MAX;
                m->b[0].gb = -FLT_MAX;
                //		fprintf(stderr,"Using this for start:%d	%d	%d\n",m->f[0].a,m->f[0].ga,m->f[0].gb);

                m->starta = old_cor[0];
                m->enda = old_cor[4]-1;

                m->startb = old_cor[2];
                m->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,m->starta,m->enda,m->startb,m->endb);
                aln_seqseq_runner(ap,seq1,seq2,m,path);
                //backward:
                m->starta = old_cor[4]+1;
                m->enda = old_cor[1];
                m->startb = c+1;
                m->endb = old_cor[3];
                m->f[0].a = 0.0;
                m->f[0].ga = -FLT_MAX;
                m->f[0].gb = -FLT_MAX;
                m->b[0].a = input_states[3];
                m->b[0].ga = input_states[4];
                m->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,m->starta,m->enda,m->startb,m->endb);
                aln_seqseq_runner(ap,seq1,seq2,m,path);
                break;
        case 2:// a -> ga = 2
                path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                m->f[0].a = input_states[0];
                m->f[0].ga = input_states[1];
                m->f[0].gb = input_states[2];
                m->b[0].a = 0.0;
                m->b[0].ga = -FLT_MAX;
                m->b[0].gb = -FLT_MAX;


                m->starta = old_cor[0];
                m->enda = old_cor[4]-1;

                m->startb = old_cor[2];
                m->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,m->starta,m->enda,m->startb,m->endb);
                aln_seqseq_runner(ap,seq1,seq2,m,path);

                //backward:
                m->starta = old_cor[4];
                m->enda = old_cor[1];
                m->startb = c+1;
                m->endb = old_cor[3];
                m->f[0].a = -FLT_MAX;
                m->f[0].ga = 0.0;
                m->f[0].gb = -FLT_MAX;
                m->b[0].a = input_states[3];
                m->b[0].ga = input_states[4];
                m->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,m->starta,m->enda,m->startb,m->endb);
                aln_seqseq_runner(ap,seq1,seq2,m,path);
                break;
        case 3:// a -> gb = 3
                path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                m->f[0].a = input_states[0];
                m->f[0].ga = input_states[1];
                m->f[0].gb = input_states[2];
                m->b[0].a = 0.0;
                m->b[0].ga = -FLT_MAX;
                m->b[0].gb = -FLT_MAX;

                m->starta = old_cor[0];
                m->enda = old_cor[4]-1;

                m->startb = old_cor[2];
                m->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,m->starta,m->enda,m->startb,m->endb);
                aln_seqseq_runner(ap,seq1,seq2,m,path);
                //backward:
                m->starta = old_cor[4]+1;
                m->enda = old_cor[1];
                m->startb = c;
                m->endb = old_cor[3];
                m->f[0].a = -FLT_MAX;
                m->f[0].ga = -FLT_MAX;
                m->f[0].gb = 0.0;
                m->b[0].a = input_states[3];
                m->b[0].ga = input_states[4];
                m->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                aln_seqseq_runner(ap,seq1,seq2,m,path);


                break;
        case 5://ga -> a = 5
                path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

                //foward:
                m->f[0].a = input_states[0];
                m->f[0].ga = input_states[1];
                m->f[0].gb = input_states[2];
                m->b[0].a = -FLT_MAX;
                m->b[0].ga = 0.0;
                m->b[0].gb = -FLT_MAX;

                m->starta = old_cor[0];
                m->enda = old_cor[4];

                m->startb = old_cor[2];
                m->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,m->starta,m->enda,m->startb,m->endb);
                aln_seqseq_runner(ap,seq1,seq2,m,path);
                //backward:
                m->starta = old_cor[4]+1;
                m->enda = old_cor[1];
                m->startb = c+1;
                m->endb = old_cor[3];
                m->f[0].a = 0.0;
                m->f[0].ga = -FLT_MAX;
                m->f[0].gb = -FLT_MAX;
                m->b[0].a = input_states[3];
                m->b[0].ga = input_states[4];
                m->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                aln_seqseq_runner(ap,seq1,seq2,m,path);
                break;
        case 6://gb->gb = 6;

                //foward:
                m->f[0].a = input_states[0];
                m->f[0].ga = input_states[1];
                m->f[0].gb = input_states[2];
                m->b[0].a = -FLT_MAX;
                m->b[0].ga = -FLT_MAX;
                m->b[0].gb = 0.0;

                m->starta = old_cor[0];
                m->enda = old_cor[4]-1;
                m->startb = old_cor[2];
                m->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,m->starta,m->enda,m->startb,m->endb);

                aln_seqseq_runner(ap,seq1,seq2,m,path);
                //backward:
                m->starta = old_cor[4]+1;
                m->enda = old_cor[1];
                m->startb = c;
                m->endb = old_cor[3];
                m->f[0].a = -FLT_MAX;
                m->f[0].ga = -FLT_MAX;
                m->f[0].gb = 0.0;
                m->b[0].a = input_states[3];
                m->b[0].ga = input_states[4];
                m->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+
                aln_seqseq_runner(ap,seq1,seq2,m,path);
                break;
        case 7://gb->a = 7;

                path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                m->f[0].a = input_states[0];
                m->f[0].ga = input_states[1];
                m->f[0].gb = input_states[2];
                m->b[0].a = -FLT_MAX;
                m->b[0].ga = -FLT_MAX;
                m->b[0].gb = 0.0;

                m->starta = old_cor[0];
                m->enda = old_cor[4]-1;
                m->startb = old_cor[2];
                m->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,m->starta,m->enda,m->startb,m->endb);

                aln_seqseq_runner(ap,seq1,seq2,m,path);
                //backward:
                m->starta = old_cor[4]+1;
                m->enda = old_cor[1];
                m->startb = c+1;
                m->endb = old_cor[3];
                m->f[0].a = 0.0;
                m->f[0].ga = -FLT_MAX;
                m->f[0].gb = -FLT_MAX;
                m->b[0].a = input_states[3];
                m->b[0].ga = input_states[4];
                m->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                aln_seqseq_runner(ap,seq1,seq2,m,path);
                break;
        default:
                break;

        }
        return OK;
}

#include <stdint.h>
#include <math.h>

#include "tldevel.h"

#include "aln_param.h"

#include "aln_struct.h"

#define ALN_SEQSEQ_IMPORT
#include "aln_seqseq.h"
#define MAX(a, b) (a > b ? a : b)
#define MAX3(a,b,c) MAX(MAX(a,b),c)

int aln_seqseq_foward(struct aln_mem* m)
{
        struct states* s = m->f;
        const uint8_t* seq1 = m->seq1;
        const uint8_t* seq2 = m->seq2;
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

        const float gpo = m->ap->gpo;
        const float gpe = m->ap->gpe;
        const float tgpe = m->ap->tgpe;
        float** subm = m->ap->subm;
        const float soff = m->ap->subm_offset;

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
                        pa += subp[seq2[j]] - soff;
                        if(m->consistency){
                                pa += m->consistency[i * m->consistency_stride + j];
                        }

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
                pa += subp[seq2[j]] - soff;
                if(m->consistency){
                        pa += m->consistency[i * m->consistency_stride + j];
                }

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
        const float soff = m->ap->subm_offset;

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

                        pa += subp[seq2[j]] - soff;
                        if(m->consistency){
                                pa += m->consistency[(starta + i) * m->consistency_stride + j];
                        }

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

                pa += subp[seq2[j]] - soff;
                if(m->consistency){
                        pa += m->consistency[(starta + i) * m->consistency_stride + j];
                }

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
        int c2 = -1;
        int transition = -1;
        int transition2 = -1;
        float s_tmp;

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
        float max2 = -FLT_MAX;
        //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
        float middle =  (float)(old_cor[3] - old_cor[2])/2.0F + (float)old_cor[2];
        float sub = 0.0F;

        //i = hm->startb;

        c = -1;
        //for(i = hm->startb; i < hm->endb;i++){
        for(i = old_cor[2]; i < old_cor[3];i++){
                sub = fabsf(middle - (float)i);
                sub /= 1000.0F;

                s_tmp = f[i].a+b[i].a-sub;
                if(s_tmp > max){
                        max2 = max; c2 = c; transition2 = transition;
                        max = s_tmp;
                        transition = 1;
                        c = i;
                }else if(s_tmp > max2){
                        max2 = s_tmp; c2 = i; transition2 = 1;
                }

                s_tmp = f[i].a+b[i].ga-gpo-sub;
                if(s_tmp > max){
                        max2 = max; c2 = c; transition2 = transition;
                        max = s_tmp;
                        transition = 2;
                        c = i;
                }else if(s_tmp > max2){
                        max2 = s_tmp; c2 = i; transition2 = 2;
                }

                s_tmp = f[i].a+b[i].gb-gpo-sub;
                if(s_tmp > max){
                        max2 = max; c2 = c; transition2 = transition;
                        max = s_tmp;
                        transition = 3;
                        c = i;
                }else if(s_tmp > max2){
                        max2 = s_tmp; c2 = i; transition2 = 3;
                }

                s_tmp = f[i].ga+b[i].a-gpo-sub;
                if(s_tmp > max){
                        max2 = max; c2 = c; transition2 = transition;
                        max = s_tmp;
                        transition = 5;
                        c = i;
                }else if(s_tmp > max2){
                        max2 = s_tmp; c2 = i; transition2 = 5;
                }


                if(m->startb == 0){
                        s_tmp = f[i].gb+b[i].gb-tgpe-sub;
                }else{
                        s_tmp = f[i].gb+b[i].gb-gpe-sub;
                }
                if(s_tmp > max){
                        max2 = max; c2 = c; transition2 = transition;
                        max = s_tmp;
                        transition = 6;
                        c = i;
                }else if(s_tmp > max2){
                        max2 = s_tmp; c2 = i; transition2 = 6;
                }

                s_tmp = f[i].gb+b[i].a-gpo-sub;
                if(s_tmp > max){
                        max2 = max; c2 = c; transition2 = transition;
                        max = s_tmp;
                        transition = 7;
                        c = i;
                }else if(s_tmp > max2){
                        max2 = s_tmp; c2 = i; transition2 = 7;
                }
        }
        //i = hm->endb;
        i = old_cor[3];
        sub = fabsf(middle - (float)i);
        sub /= 1000.0F;

        s_tmp = f[i].a+b[i].gb-gpo-sub;
        if(s_tmp > max){
                max2 = max; c2 = c; transition2 = transition;
                max = s_tmp;
                transition = 3;
                c = i;
        }else if(s_tmp > max2){
                max2 = s_tmp; c2 = i; transition2 = 3;
        }

        if(m->endb == m->len_b){
                s_tmp = f[i].gb+b[i].gb-tgpe-sub;
        }else{
                s_tmp = f[i].gb+b[i].gb-gpe-sub;
        }
        if(s_tmp > max){
                max2 = max; c2 = c; transition2 = transition;
                max = s_tmp;
                transition = 6;
                c = i;
        }else if(s_tmp > max2){
                max2 = s_tmp; c2 = i; transition2 = 6;
        }

        /* Accumulate confidence margin and record per-meetup margins */
        if(max2 > -FLT_MAX){
                float _margin = max - max2;
                if(m->flip_margins != NULL && m->margin_count < m->flip_margin_alloc){
                        m->flip_margins[m->margin_count] = _margin;
                }
                m->margin_sum += _margin;
                m->margin_count++;
        }

        /* Perturbation: flip uncertain midpoints to second-best choice.
           Three modes: individual targeting (MCTS), stride bitmask, or round-robin. */
        if(m->flip_threshold > 0.0F && c2 >= 0 && max2 > -FLT_MAX){
                float margin = max - max2;
                if(margin < m->flip_threshold){
                        if(m->flip_bit_map != NULL){
                                /* Individual midpoint targeting (MCTS) */
                                if(m->flip_counter < m->flip_n_uncertain){
                                        int bit = m->flip_bit_map[m->flip_counter];
                                        if(bit >= 0 && ((1U << bit) & m->flip_mask)){
                                                c = c2;
                                                transition = transition2;
                                        }
                                }
                        }else if(m->flip_mask != 0){
                                /* Stride-based bitmask mode */
                                if((1U << (m->flip_counter % m->flip_stride)) & m->flip_mask){
                                        c = c2;
                                        transition = transition2;
                                }
                        }else if(m->flip_trial > 0){
                                /* Round-robin mode */
                                if(m->flip_counter % m->flip_stride == m->flip_trial - 1){
                                        c = c2;
                                        transition = transition2;
                                }
                        }
                        m->flip_counter++;
                }
        }

        *meet = c;
        *t = transition;
        *score = max;
        return OK;
}

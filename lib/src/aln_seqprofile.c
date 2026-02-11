#include <stdint.h>
#include <math.h>

#include "tldevel.h"

#include "aln_param.h"
#include "aln_struct.h"
#define ALN_SEQPROFILE_IMPORT
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
        int c2 = -1;
        int transition = -1;
        int transition2 = -1;
        float s_tmp;

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
        float max2 = -FLT_MAX;
        //float middle =  (m->endb - m->startb)/2 + m->startb;
        float middle = (float)(old_cor[3] - old_cor[2])/2.0F + (float)old_cor[2];
        float sub = 0.0F;


        prof1+= ((old_cor[4]+1)<<6);

        //i = m->startb;

        c = -1;
        //for(i = m->startb; i < m->endb;i++){
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

                s_tmp = f[i].a+b[i].ga-open-sub;
                if(s_tmp > max){
                        max2 = max; c2 = c; transition2 = transition;
                        max = s_tmp;
                        transition = 2;
                        c = i;
                }else if(s_tmp > max2){
                        max2 = s_tmp; c2 = i; transition2 = 2;
                }

                s_tmp = f[i].a+b[i].gb+prof1[27]-sub;
                if(s_tmp > max){
                        max2 = max; c2 = c; transition2 = transition;
                        max = s_tmp;
                        transition = 3;
                        c = i;
                }else if(s_tmp > max2){
                        max2 = s_tmp; c2 = i; transition2 = 3;
                }

                s_tmp = f[i].ga+b[i].a-open-sub;
                if(s_tmp > max){
                        max2 = max; c2 = c; transition2 = transition;
                        max = s_tmp;
                        transition = 5;
                        c = i;
                }else if(s_tmp > max2){
                        max2 = s_tmp; c2 = i; transition2 = 5;
                }


                if(m->startb == 0){
                        s_tmp = f[i].gb+b[i].gb+prof1[29]-sub;
                }else{
                        s_tmp = f[i].gb+b[i].gb+prof1[28]-sub;
                }
                if(s_tmp > max){
                        max2 = max; c2 = c; transition2 = transition;
                        max = s_tmp;
                        transition = 6;
                        c = i;
                }else if(s_tmp > max2){
                        max2 = s_tmp; c2 = i; transition2 = 6;
                }

                s_tmp = f[i].gb+b[i].a+prof1[-37]-sub;
                if(s_tmp > max){
                        max2 = max; c2 = c; transition2 = transition;
                        max = s_tmp;
                        transition = 7;
                        c = i;
                }else if(s_tmp > max2){
                        max2 = s_tmp; c2 = i; transition2 = 7;
                }
        }
        //i = m->endb;
        i = old_cor[3];

        sub = fabsf(middle - (float)i);
        sub /= 1000.0F;

        s_tmp = f[i].a+b[i].gb+prof1[27]-sub;
        if(s_tmp > max){
                max2 = max; c2 = c; transition2 = transition;
                max = s_tmp;
                transition = 3;
                c = i;
        }else if(s_tmp > max2){
                max2 = s_tmp; c2 = i; transition2 = 3;
        }

        if(m->endb == m->len_b){
                s_tmp = f[i].gb+b[i].gb+prof1[29]-sub;
        }else{
                s_tmp = f[i].gb+b[i].gb+prof1[28]-sub;
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

#include "tldevel.h"

#include <stdint.h>

#include "aln_struct.h"
#include "alignment_parameters.h"

#include "aln_mem.h"

#define ALN_SETUP_IMPORT
#include "aln_setup.h"

int init_alnmem(struct aln_mem* m)
{
        int i;
        int g;
        m->starta = 0;
        m->startb = 0;
        m->enda = m->len_a;
        m->endb = m->len_b;

        m->f[0].a  = 0.0F;
        m->f[0].ga = -FLT_MAX;
        m->f[0].gb = -FLT_MAX;
        m->b[0].a  = 0.0F;
        m->b[0].ga = -FLT_MAX;
        m->b[0].gb = -FLT_MAX;

        RUN(resize_aln_mem(m));

        g = MACRO_MAX(m->len_a, m->len_b) + 2;
        for(i = 0 ;i  < g ;i++){
                m->path[i] = -1;
        }
        return OK;
ERROR:
        return FAIL;
}

int make_profile_n(struct aln_param* ap,const uint8_t* seq,const int len, float** p)
{
        float** subm = NULL;
        float* prof = NULL;
        float gpo;
        float gpe;
        float tgpe;
        int i;
        int j;
        int c;
        gpo = ap->gpo;
        gpe = ap->gpe;
        tgpe = ap->tgpe;
        subm = ap->subm;

        MMALLOC(prof,sizeof(float)*(len+2)*64);
        prof +=  (64 *(len+1));

        for (i = 0;i < 64;i++){
                prof[i] = 0;
        }
        prof[23+32] = -gpo;
        prof[24+32] = -gpe;
        prof[25+32] = -tgpe;

        i = len;
        while(i--){
                prof -= 64;

                for (j = 0;j < 64;j++){
                        prof[j] = 0;
                }
                c = seq[i];

                prof[c] += 1;

                prof += 32;

                for(j = 23;j--;){
                        prof[j] = subm[c][j];
                }
                prof[23] = -gpo;
                prof[24] = -gpe;
                prof[25] = -tgpe;

                prof -= 32;
        }
        prof -= 64;
        for (i = 0;i < 64;i++){
                prof[i] = 0;
        }
        prof[23+32] = -gpo;
        prof[24+32] = -gpe;
        prof[25+32] = -tgpe;
        *p = prof;
        return OK;
ERROR:
        return FAIL;
}

int set_gap_penalties_n(float* prof,int len,int nsip)
{
        int i;

        prof +=  (64 *(len+1));

        prof[27] = prof[55] * (float)nsip;//gap open or close  23
        prof[28] = prof[56] * (float)nsip;//gap extention 24
        prof[29] = prof[57] * (float)nsip;//gap open or close 25

        i = len+1;
        while(i--){
                prof -= 64;
                prof[27] = prof[55] * (float)nsip;//gap open or close
                prof[28] = prof[56] * (float)nsip;//gap extention
                prof[29] = prof[57] * (float)nsip;//gap open or close
        }
        return OK;
}

int add_gap_info_to_path_n(struct aln_mem* m)
{
        int* path = NULL;
        int* o_path = NULL;
        int* tmp_path = NULL;
        int i,j;
        int a = 0;
        int b = 0;
        int len_a;
        int len_b;

        len_a = m->len_a;
        len_b = m->len_b;
        path = m->path;
        o_path = m->tmp_path;

        for(i = 0; i < len_a+len_b+2;i++){
                o_path[i] = 0;
        }

        j = 1;
        b = -1;
        if(path[1] == -1){
                o_path[j] = 2;
                j++;
        }else{
                if(path[1] != 1){
                        for ( a = 0;a < path[1] -1;a++){
                                o_path[j] = 1;
                                j++;
                        }
                        o_path[j] = 0;
                        j++;
                }else{
                        o_path[j] = 0;
                        j++;
                }
        }
        b = path[1];

        for(i = 2; i <= len_a;i++){

                if(path[i] == -1){
                        o_path[j] = 2;
                        j++;
                }else{
                        if(path[i]-1 != b && b != -1){
                                for ( a = 0;a < path[i] - b-1;a++){
                                        o_path[j] = 1;
                                        j++;
                                }
                                o_path[j] = 0;
                                j++;
                        }else{
                                o_path[j] = 0;
                                j++;
                        }
                }
                b = path[i];
        }

        if(path[len_a] < len_b && path[len_a] != -1){
                //	fprintf(stderr,"WARNING:%d	%d\n",path[len_a],len_b);
                for ( a = 0;a < len_b - path[len_a];a++){
                        o_path[j] = 1;
                        j++;
                }
        }

        o_path[0] = j-1;
        o_path[j] = 3;

        //add gap info..
        i = 2;
        while(o_path[j] != 3){
                if ((o_path[i-1] &3) && !(o_path[i] & 3)){
                        if(o_path[i-1] & 8){
                                o_path[i-1] += 8;
                        }else{
                                o_path[i-1] |= 16;
                        }
                }else if (!(o_path[i-1] & 3) &&(o_path[i] &3)){
                        o_path[i] |= 4;
                }else if ((o_path[i-1] & 1) && (o_path[i] & 1)){
                        o_path[i] |= 8;
                }else if ((o_path[i-1] & 2) && (o_path[i] & 2)){
                        o_path[i] |= 8;
                }
                i++;
        }
        //add terminal gap...
        i = 1;
        while(o_path[i] != 0){
                o_path[i] |= 32;
                i++;
        }
        /* j = i; */
        i = o_path[0];
        while(o_path[i] != 0){
                o_path[i] |= 32;
                i--;
        }

        tmp_path = m->path;
        m->path = m->tmp_path;
        m->tmp_path = tmp_path;
        return OK;
}

int update_n(const float* profa, const float* profb,float* newp, struct aln_param*ap, int* path,int sipa,int sipb)
{
        float gp;
        int i;
        int j;
        int c;
        for (i = 64; i--;){
                newp[i] = profa[i] + profb[i];
        }

        profa += 64;
        profb += 64;
        newp += 64;

        c = 1;

        while(path[c] != 3){
                //Idea: limit the 'virtual' number of residues of one type to x.
                // i.e. only allow a maximum of 10 alanines to be registered in each column
                // the penalty for aligning a 'G' to this column will stay stable even when many (>10) alanines are present.
                // the difference in score between the 'correct' (all alanine) and incorrect (alanines + glycine) will not increase
                // with the number of sequences. -> see Durbin pp 140

                if (!path[c]){
                        //fprintf(stderr,"Align	%d\n",c);
                        for (i = 64; i--;){
                                newp[i] = profa[i] + profb[i];
                        }
                        profa += 64;
                        profb += 64;
                }

                if (path[c] & 1){
                        //fprintf(stderr,"Gap_A:%d\n",c);
                        //printf("open:%d	ext:%d	%d	%d\n",si->nsip[a] * gpo,si->nsip[a] * gpe,si->nsip[a] * profb[41],si->nsip[a] * profb[46]);
                        for (i = 64; i--;){
                                newp[i] = profb[i];
                        }
                        profb += 64;
                        if(!(path[c] & 20)){
                                if(path[c] & 32){
                                        newp[25] += (float)sipa;//1;
                                        gp = ap->tgpe*(float)sipa;
                                }else{
                                        newp[24] += (float)sipa;//1;
                                        gp = ap->gpe*(float)sipa;
                                }

                                for (j = 32; j < 55;j++){
                                        newp[j] -= gp;
                                }
                        }else{
                                if (path[c] & 16){
                                        //			fprintf(stderr,"close_open");
                                        if(path[c] & 32){
                                                newp[25] += (float)sipa;//1;
                                                gp = ap->tgpe*(float)sipa;
                                                newp[23] += (float)sipa;//1;
                                                gp += ap->gpo*(float)sipa;
                                        }else{
                                                newp[23] += (float)sipa;//1;
                                                gp = ap->gpo*(float)sipa;
                                        }

                                        for (j = 32; j < 55;j++){
                                                newp[j] -= gp;
                                        }
                                }
                                if (path[c] & 4){
                                        //			fprintf(stderr,"Gap_open");
                                        if(path[c] & 32){
                                                newp[25] += (float)sipa;//1;
                                                gp = ap->tgpe*(float)sipa;
                                                newp[23] += (float)sipa;//1;
                                                gp += ap->gpo*(float)sipa;
                                        }else{
                                                newp[23] += (float)sipa;//1;
                                                gp = ap->gpo*(float)sipa;
                                        }
                                        for (j = 32; j < 55;j++){
                                                newp[j] -= gp;
                                        }
                                }
                        }
                }
                if (path[c] & 2){
                        //fprintf(stderr,"Gap_B:%d\n",c);
                        //printf("open:%d	ext:%d	%d	%d\n",si->nsip[b] * gpo,si->nsip[b] * gpe,profa[26],profa[27]);
                        for (i = 64; i--;){
                                newp[i] = profa[i];
                        }
                        profa+=64;
                        if(!(path[c] & 20)){
                                if(path[c] & 32){
                                        newp[25] += sipb;//1;
                                        gp = ap->tgpe*sipb;
                                }else{
                                        newp[24] += sipb;//1;
                                        gp = ap->gpe*sipb;
                                }
                                for (j = 32; j < 55;j++){
                                        newp[j] -= gp;
                                }
                        }else{
                                if (path[c] & 16){
                                        //			fprintf(stderr,"close_open");
                                        if(path[c] & 32){
                                                newp[25] += (float)sipb;//1;
                                                gp =  ap->tgpe*(float)sipb;
                                                newp[23] += (float)sipb;//1;
                                                gp +=  ap->gpo*(float)sipb;
                                        }else{
                                                newp[23] += (float)sipb;//1;
                                                gp =  ap->gpo*(float)sipb;
                                        }
                                        for (j = 32; j < 55;j++){
                                                newp[j] -= gp;
                                        }
                                }
                                if (path[c] & 4){
                                        //			fprintf(stderr,"Gap_open");
                                        if(path[c] & 32){
                                                newp[25] += (float)sipb;//1;
                                                gp = ap->tgpe*(float)sipb;
                                                newp[23] += (float)sipb;//1;
                                                gp += ap->gpo*(float)sipb;
                                        }else{
                                                newp[23] += (float)sipb;//1;
                                                gp = ap->gpo*(float)sipb;
                                        }

                                        for (j = 32; j < 55;j++){
                                                newp[j] -= gp;
                                        }
                                }
                        }
                }
                newp += 64;
                c++;
        }
        for (i = 64; i--;){
                newp[i] =  profa[i] + profb[i];
        }
        //newp -= (path[0]+1) *64;
        return OK;
}

int mirror_path_n(struct aln_mem* m,int len_a,int len_b)
{
        int* apath = NULL;
        int* opath = NULL;
        int* tmppath = NULL;
        int i;

        apath = m->path;
        opath = m->tmp_path;

        for(i =0; i < len_a+2;i++){
                opath[i] = -1;
        }

        for(i = 1; i <= len_b;i++){
                if(apath[i] != -1){
                        opath[apath[i]] = i;
                }
        }

        tmppath = m->path;
        m->path = m->tmp_path;
        m->tmp_path = tmppath;
        return OK;
}

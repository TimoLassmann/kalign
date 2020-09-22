#include "tldevel.h"

#include <stdint.h>

#include "aln_struct.h"
#include "alignment_parameters.h"


#define ALN_SETUP_IMPORT
#include "aln_setup.h"


int init_alnmem(struct aln_mem* m, int len_a, int len_b)
{
        m->starta = 0;
        m->startb = 0;
        m->enda = len_a;
        m->endb = len_b;
        m->len_a = len_a;
        m->len_b = len_b;

        m->f[0].a  = 0.0F;
        m->f[0].ga = -FLT_MAX;
        m->f[0].gb = -FLT_MAX;
        m->b[0].a  = 0.0F;
        m->b[0].ga = -FLT_MAX;
        m->b[0].gb = -FLT_MAX;
        return OK;
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

                for(j = 21;j--;){
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

        prof[27] = prof[55]*(float)nsip;//gap open or close  23
        prof[28] = prof[56]*(float)nsip;//gap extention 24
        prof[29] = prof[57]*(float)nsip;//gap open or close 25

        i = len+1;
        while(i--){
                prof -= 64;
                prof[27] = prof[55]*(float)nsip;//gap open or close
                prof[28] = prof[56]*(float)nsip;//gap extention

                prof[29] = prof[57]*(float)nsip;//gap open or close
        }
        return OK;
}

int add_gap_info_to_path_n(int** p_path,int len_a,int len_b)
{
        int* path = NULL;
        int i,j;
        int a = 0;
        int b = 0;

        int* np = NULL;

        path = *p_path;

        MMALLOC(np,sizeof(int)*(len_a+len_b+2));
        for(i =0; i < len_a+len_b+2;i++){
                np[i] = 0;
        }

        j = 1;
        b = -1;
        if(path[1] == -1){
                np[j] = 2;
                j++;
        }else{
                if(path[1] != 1){
                        for ( a = 0;a < path[1] -1;a++){
                                np[j] = 1;
                                j++;
                        }
                        np[j] = 0;
                        j++;
                }else{
                        np[j] = 0;
                        j++;
                }
        }
        b = path[1];

        /*for ( i= 0;i <= len_a;i++){
          fprintf(stderr,"%d,",path[i]);
          }
          fprintf(stderr,"\n");*/

        for(i = 2; i <= len_a;i++){

                if(path[i] == -1){
                        np[j] = 2;
                        j++;
                }else{
                        if(path[i]-1 != b && b != -1){
                                for ( a = 0;a < path[i] - b-1;a++){
                                        np[j] = 1;
                                        j++;
                                }
                                np[j] = 0;
                                j++;
                        }else{
                                np[j] = 0;
                                j++;
                        }
                }
                b = path[i];
        }





        if(path[len_a] < len_b && path[len_a] != -1){
                //	fprintf(stderr,"WARNING:%d	%d\n",path[len_a],len_b);
                for ( a = 0;a < len_b - path[len_a];a++){
                        np[j] = 1;
                        j++;
                }
        }
        np[0] = j-1;
        np[j] = 3;

        MREALLOC(np,sizeof(int)* (np[0]+2));
        //for ( i= 0;i <= np[0];i++){
        //	fprintf(stderr,"%d,",np[i]);
        //}
        //fprintf(stderr,"\n");

        MFREE(path);

        //add gap info..
        i = 2;
        while(np[i] != 3){
                if ((np[i-1] &3) && !(np[i] & 3)){
                        if(np[i-1] & 8){
                                np[i-1] += 8;
                        }else{
                                np[i-1] |= 16;
                        }
                }else if (!(np[i-1] & 3) &&(np[i] &3)){
                        np[i] |= 4;
                }else if ((np[i-1] & 1) && (np[i] & 1)){
                        np[i] |= 8;
                }else if ((np[i-1] & 2) && (np[i] & 2)){
                        np[i] |= 8;
                }
                i++;
        }
        //add terminal gap...
        i = 1;
        while(np[i] != 0){
                np[i] |= 32;
                i++;
        }
        j = i;
        i = np[0];
        while(np[i] != 0){
                np[i] |= 32;
                i--;
        }
        //for ( i= 0;i <= np[0];i++){
        //	fprintf(stderr,"%d,",np[i]);
        //}
        //fprintf(stderr,"\n");
        *p_path = np;
        return OK;
ERROR:
        if(np){
                MFREE(np);
        }
        return FAIL;
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

int mirror_path_n(int** p_path,int len_a,int len_b)
{
        int* path = NULL;
        int* np =NULL;

        int i;

        path = *p_path;

        MMALLOC(np,sizeof(int)*(len_a+2));
        for(i =0; i < len_a+2;i++){
                np[i] = -1;
        }

        for(i = 1; i <= len_b;i++){
                if(path[i] != -1){
                        np[path[i]] = i;
                }
        }

        MFREE(path);
        *p_path = np;
        return OK;
ERROR:
        return FAIL;
}

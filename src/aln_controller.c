#include <stdint.h>

#include "tldevel.h"

#include "alignment_parameters.h"

#include "aln_struct.h"

/* The dynamic programming modules */
#include "aln_seqseq.h"
#include "aln_seqprofile.h"
#include "aln_profileprofile.h"

#define ALN_CONTROLLER_IMPORT
#include "aln_controller.h"

//static int aln_continue(struct aln_mem* m,struct aln_param* ap,int* path,int meet,int transition);

static int aln_continue(struct aln_mem* m,struct aln_param* ap,float input_states[],int old_cor[],int* path,int meet,int transition);

int aln_runner(struct aln_mem* m,struct aln_param* ap, int* path)
{
        float input_states[6];
        int old_cor[5];
        float score;
        int mid;
        int meet;
        int transition;

        input_states[0] = m->f[0].a;
        input_states[1] = m->f[0].ga;
        input_states[2] = m->f[0].gb;
        input_states[3] = m->b[0].a;
        input_states[4] = m->b[0].ga;
        input_states[5] = m->b[0].gb;

        mid = ((m->enda - m->starta) / 2)+ m->starta;

        old_cor[0] = m->starta;
        old_cor[1] = m->enda;
        old_cor[2] = m->startb;
        old_cor[3] = m->endb;
        old_cor[4] = mid;

        if(m->starta >= m->enda){
                return OK;//hirsch_path;
        }
        if(m->startb >= m->endb){
                return OK;///hirsch_path;
        }

        //fprintf(stderr,"Forward:%d-%d	%d-%d\n",m->starta,m->enda,m->startb,m->endb);

        m->enda = mid;

        //fprintf(stderr,"Forward:%d-%d	%d-%d\n",m->starta,m->enda,m->startb,m->endb);
        if(ap->seq1){
                aln_seqseq_foward(m,ap);
        }else if(ap->prof2){
                aln_profileprofile_foward(m,ap);
        }else{
                aln_seqprofile_foward(m,ap);
        }

        m->starta = mid;
        m->enda = old_cor[1];
        //fprintf(stderr,"Backward:%d-%d	%d-%d\n",m->starta,m->enda,m->startb,m->endb);


        if(ap->seq1){
                aln_seqseq_backward(m,ap);
                aln_seqseq_meetup(m,ap,old_cor,&meet,&transition,&score);
        }else if(ap->prof2){
                aln_profileprofile_backward(m,ap);
                aln_profileprofile_meetup(m,ap,old_cor,&meet,&transition,&score);
        }else{
                aln_seqprofile_backward(m,ap);
                aln_seqprofile_meetup(m,ap,old_cor,&meet,&transition,&score);
        }

        if(ap->mode == ALN_MODE_SCORE_ONLY){
                ap->score = score;
        }else{
                aln_continue(m, ap,input_states,old_cor,path, meet, transition);
        }
        return OK;
}

int aln_continue(struct aln_mem* m,struct aln_param* ap,float input_states[],int old_cor[],int* path,int meet,int transition)
{
        //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
        //LOG_MSG("MAX: %f",max);
        //j = hirsch_path[0];
        switch(transition){
        case 1: //a -> a = 1

                path[old_cor[4]] = meet;
                path[old_cor[4]+1] = meet+1;

                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                m->f[0].a = input_states[0];
                m->f[0].ga = input_states[1];
                m->f[0].gb = input_states[2];
                m->b[0].a = 0.0F;
                m->b[0].ga = -FLT_MAX;
                m->b[0].gb = -FLT_MAX;
                //		fprintf(stderr,"Using this for start:%d	%d	%d\n",m->f[0].a,m->f[0].ga,m->f[0].gb);

                m->starta = old_cor[0];
                m->enda = old_cor[4]-1;

                m->startb = old_cor[2];
                m->endb = meet-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,m->starta,m->enda,m->startb,m->endb);
                aln_runner(m,ap,path);
                //backward:
                m->starta = old_cor[4]+1;
                m->enda = old_cor[1];
                m->startb = meet+1;
                m->endb = old_cor[3];
                m->f[0].a = 0.0F;
                m->f[0].ga = -FLT_MAX;
                m->f[0].gb = -FLT_MAX;
                m->b[0].a = input_states[3];
                m->b[0].ga = input_states[4];
                m->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,m->starta,m->enda,m->startb,m->endb);
                aln_runner(m,ap,path);
                break;
        case 2:// a -> ga = 2
                path[old_cor[4]] = meet;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                m->f[0].a = input_states[0];
                m->f[0].ga = input_states[1];
                m->f[0].gb = input_states[2];
                m->b[0].a = 0.0F;
                m->b[0].ga = -FLT_MAX;
                m->b[0].gb = -FLT_MAX;


                m->starta = old_cor[0];
                m->enda = old_cor[4]-1;

                m->startb = old_cor[2];
                m->endb = meet-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,m->starta,m->enda,m->startb,m->endb);
                aln_runner(m,ap,path);

                //backward:
                m->starta = old_cor[4];
                m->enda = old_cor[1];
                m->startb = meet+1;
                m->endb = old_cor[3];
                m->f[0].a = -FLT_MAX;
                m->f[0].ga = 0.0F;
                m->f[0].gb = -FLT_MAX;
                m->b[0].a = input_states[3];
                m->b[0].ga = input_states[4];
                m->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,m->starta,m->enda,m->startb,m->endb);
                aln_runner(m,ap,path);
                break;
        case 3:// a -> gb = 3
                path[old_cor[4]] = meet;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                m->f[0].a = input_states[0];
                m->f[0].ga = input_states[1];
                m->f[0].gb = input_states[2];
                m->b[0].a = 0.0F;
                m->b[0].ga = -FLT_MAX;
                m->b[0].gb = -FLT_MAX;

                m->starta = old_cor[0];
                m->enda = old_cor[4]-1;

                m->startb = old_cor[2];
                m->endb = meet-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,m->starta,m->enda,m->startb,m->endb);
                aln_runner(m,ap,path);
                //backward:
                m->starta = old_cor[4]+1;
                m->enda = old_cor[1];
                m->startb = meet;
                m->endb = old_cor[3];
                m->f[0].a = -FLT_MAX;
                m->f[0].ga = -FLT_MAX;
                m->f[0].gb = 0.0;
                m->b[0].a = input_states[3];
                m->b[0].ga = input_states[4];
                m->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                aln_runner(m,ap,path);


                break;
        case 5://ga -> a = 5
                path[old_cor[4]+1] = meet+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

                //foward:
                m->f[0].a = input_states[0];
                m->f[0].ga = input_states[1];
                m->f[0].gb = input_states[2];
                m->b[0].a = -FLT_MAX;
                m->b[0].ga = 0.0F;
                m->b[0].gb = -FLT_MAX;

                m->starta = old_cor[0];
                m->enda = old_cor[4];

                m->startb = old_cor[2];
                m->endb = meet-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,m->starta,m->enda,m->startb,m->endb);
                aln_runner(m,ap,path);
                //backward:
                m->starta = old_cor[4]+1;
                m->enda = old_cor[1];
                m->startb = meet+1;
                m->endb = old_cor[3];
                m->f[0].a = 0.0F;
                m->f[0].ga = -FLT_MAX;
                m->f[0].gb = -FLT_MAX;
                m->b[0].a = input_states[3];
                m->b[0].ga = input_states[4];
                m->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                aln_runner(m,ap,path);
                break;
        case 6://gb->gb = 6;

                //foward:
                m->f[0].a = input_states[0];
                m->f[0].ga = input_states[1];
                m->f[0].gb = input_states[2];
                m->b[0].a = -FLT_MAX;
                m->b[0].ga = -FLT_MAX;
                m->b[0].gb = 0.0F;

                m->starta = old_cor[0];
                m->enda = old_cor[4]-1;
                m->startb = old_cor[2];
                m->endb = meet;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,m->starta,m->enda,m->startb,m->endb);

                aln_runner(m,ap,path);
                //backward:
                m->starta = old_cor[4]+1;
                m->enda = old_cor[1];
                m->startb = meet;
                m->endb = old_cor[3];
                m->f[0].a = -FLT_MAX;
                m->f[0].ga = -FLT_MAX;
                m->f[0].gb = 0.0F;
                m->b[0].a = input_states[3];
                m->b[0].ga = input_states[4];
                m->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+
                aln_runner(m,ap,path);
                break;
        case 7://gb->a = 7;

                path[old_cor[4]+1] = meet+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                m->f[0].a = input_states[0];
                m->f[0].ga = input_states[1];
                m->f[0].gb = input_states[2];
                m->b[0].a = -FLT_MAX;
                m->b[0].ga = -FLT_MAX;
                m->b[0].gb = 0.0F;

                m->starta = old_cor[0];
                m->enda = old_cor[4]-1;
                m->startb = old_cor[2];
                m->endb = meet;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,m->starta,m->enda,m->startb,m->endb);

                aln_runner(m,ap,path);
                //backward:
                m->starta = old_cor[4]+1;
                m->enda = old_cor[1];
                m->startb = meet+1;
                m->endb = old_cor[3];
                m->f[0].a = 0.0F;
                m->f[0].ga = -FLT_MAX;
                m->f[0].gb = -FLT_MAX;
                m->b[0].a = input_states[3];
                m->b[0].ga = input_states[4];
                m->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                aln_runner(m,ap,path);
                break;
        default:
                break;

        }
        return OK;
}

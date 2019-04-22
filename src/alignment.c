#include "alignment.h"




int align(struct alignment* aln, struct aln_param* ap)
{

        ASSERT(aln!= NULL, "No alignment.");
        ASSERT(ap != NULL, "No parameters.");

        if(aln->dna){

        }

        return OK;
ERROR:
        return FAIL;
}

int** dna_alignment(struct alignment* aln, struct aln_param* ap)
{
        struct hirsch_mem* hm = NULL;
        int i,j,g,a,b,c;
        int len_a;
        int len_b;
        float** profile = NULL;

        profile = malloc(sizeof(float*)*numprofiles);
        for ( i = 0;i< numprofiles;i++){
                profile[i] = 0;
        }

        map = malloc(sizeof(int*)*numprofiles);
        for ( i = 0;i < numprofiles;i++){
                map[i] = 0;
        }

        hm = hirsch_mem_alloc(hm,1024);
        fprintf(stderr,"\nAlignment:\n");
        for (i = 0; i < (numseq-1);i++){
                a = tree[i*3];
                b = tree[i*3+1];
                c = tree[i*3+2];
                fprintf(stderr,"\r%8.0f percent done",(float)(i) /(float)numseq * 100);
                //fprintf(stderr,"Aligning:%d %d->%d	done:%0.2f\n",a,b,c,((float)(i+1)/(float)numseq)*100);
                len_a = aln->sl[a];
                len_b = aln->sl[b];


                g = (len_a > len_b)? len_a:len_b;
                map[c] = malloc(sizeof(int) * (g+2));
                if(g > hm->size){
                        hm = hirsch_mem_realloc(hm,g);
                }

                for (j = 0; j < (g+2);j++){
                        map[c][j] = -1;
                }

                if (a < numseq){
                        profile[a] = dna_make_profile(profile[a],aln->s[a],len_a,submatrix);
                }
                if (b < numseq){
                        profile[b] = dna_make_profile(profile[b],aln->s[b],len_b,submatrix);
                }
                fprintf(stderr,"Saving mem...\n");

                dna_set_gap_penalties(profile[a],len_a,aln->nsip[b],strength,aln->nsip[a]);
                dna_set_gap_penalties(profile[b],len_b,aln->nsip[a],strength,aln->nsip[b]);

                hm->starta = 0;
                hm->startb = 0;
                hm->enda = len_a;
                hm->endb = len_b;
                hm->len_a = len_a;
                hm->len_b = len_b;

                hm->f[0].a = 0.0;
                hm->f[0].ga =  -FLT_MAX;
                hm->f[0].gb = -FLT_MAX;
                hm->b[0].a = 0.0;
                hm->b[0].ga =  -FLT_MAX;
                hm->b[0].gb =  -FLT_MAX;
                //	fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
                if(a < numseq){
                        if(b < numseq){
                                map[c] = hirsch_dna_ss_dyn(submatrix,aln->s[a],aln->s[b],hm,map[c]);
                        }else{
                                hm->enda = len_b;
                                hm->endb = len_a;
                                hm->len_a = len_b;
                                hm->len_b = len_a;
                                map[c] = hirsch_dna_ps_dyn(profile[b],aln->s[a],hm,map[c],aln->nsip[b]);
                                map[c] = mirror_hirsch_path(map[c],len_a,len_b);
                        }
                }else{
                        if(b < numseq){
                                map[c] = hirsch_dna_ps_dyn(profile[a],aln->s[b],hm,map[c],aln->nsip[a]);
                        }else{
                                if(len_a < len_b){
                                        map[c] = hirsch_dna_pp_dyn(profile[a],profile[b],hm,map[c]);
                                }else{
                                        hm->enda = len_b;
                                        hm->endb = len_a;
                                        hm->len_a = len_b;
                                        hm->len_b = len_a;
                                        map[c] = hirsch_dna_pp_dyn(profile[b],profile[a],hm,map[c]);
                                        map[c] = mirror_hirsch_path(map[c],len_a,len_b);
                                }
                        }
                }
                map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

                if(i != numseq-2){
                        profile[c] = malloc(sizeof(float)*22*(map[c][0]+2));
                        profile[c] = dna_update(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);
                }

                aln->sl[c] = map[c][0];

                aln->nsip[c] = aln->nsip[a] + aln->nsip[b];
                aln->sip[c] = malloc(sizeof(int)*(aln->nsip[a] + aln->nsip[b]));
                g =0;
                for (j = aln->nsip[a];j--;){
                        aln->sip[c][g] = aln->sip[a][j];
                        g++;
                }
                for (j = aln->nsip[b];j--;){
                        aln->sip[c][g] = aln->sip[b][j];
                        g++;
                }

                free(profile[a]);
                free(profile[b]);
        }
        //free(profile[numprofiles-1]);
        free(profile);
        hirsch_mem_free(hm);
        return map;
}


struct hirsch_mem* hirsch_mem_alloc(struct hirsch_mem* hm,int x)
{

        struct hirsch_mem* hm = NULL;
        // a=((typeof(a))(((int)(((void *)malloc(c+15))+15))&-16)).
        MMALLOC(hm,sizeof(struct hirsch_mem));
        hm->starta = 0;
        hm->startb = 0;
        hm->enda = 0;
        hm->endb = 0;
        hm->size = x+1;
        hm->len_a = 0;
        hm->len_b = 0;
        hm->f = NULL;
        hm->b = NULL;
        MMALLOC(hm->f,sizeof(struct states)* (x+1));
        MMALLOC(hm->b,sizeof(struct states)* (x+1));
        return hm;
ERROR:
        return NULL;
}

int hirsch_mem_realloc(struct hirsch_mem* hm,int x)
{

        if((x+1) > hm->size){
                hm->size = x+1;

                MREALLOC(hm->f,sizeof(struct states)* (x+1));
                MREALLOC(hm->b,sizeof(struct states)* (x+1));
        }
        return OK;
ERROR:
        return FAIL;
}

void hirsch_mem_free(struct hirsch_mem* hm)
{
        if(hm){
                MFREE(hm->f);
                MFREE(hm->b);
                MFREE(hm);
        }
}

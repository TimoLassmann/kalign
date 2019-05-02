#include "weave_alignment.h"

int clean_aln(struct alignment* aln);
struct alignment* make_seq(struct alignment* aln,int a,int b,int* path);
int update_gaps(int old_len,int*gis,int new_len,int *newgaps);

int weave(struct alignment* aln, int** map, int* tree)
{
        int i;
        int a,b;

        RUN(clean_aln(aln));

        for (i = 0; i < (aln->numseq-1)*3;i +=3){
                a = tree[i];
                b = tree[i+1];
                aln = make_seq(aln,a,b,map[tree[i+2]]);
        }

        for (i = 0; i < aln->numseq;i++){
                aln->nsip[i] = i;
        }
        return OK;
ERROR:
        return FAIL;
}

int clean_aln(struct alignment* aln)
{
        int i,j;
        int* p = NULL;
        for (i = 0; i < aln->numseq;i++){
                p = aln->gaps[i];
                for (j = 0; j <= aln->sl[i];j++){
                        p[j] = 0;
                }
        }
        return OK;

}

struct alignment* make_seq(struct alignment* aln,int a,int b,int* path)
{
        int* gap_a = NULL;
        int* gap_b = NULL;

        int c;
        int i;
        int posa = 0;
        int posb = 0;


        MMALLOC(gap_a,(path[0]+1)*sizeof(int));
        MMALLOC(gap_b,(path[0]+1)*sizeof(int));

        for (i = path[0]+1;i--;){
                gap_a[i] = 0;
                gap_b[i] = 0;
        }
        c = 1;
        while(path[c] != 3){
                if (!path[c]){
                        posa++;
                        posb++;
                }else
                if (path[c] & 1){
                        gap_a[posa] += 1;
                        posb++;
                }else
                if (path[c] & 2){
                        gap_b[posb] += 1;
                        posa++;
                }
                c++;
        }
        for (i = aln->nsip[a];i--;){
                RUN(update_gaps(aln->sl[aln->sip[a][i]],aln->gaps[aln->sip[a][i]],path[0],gap_a));
        }
        for (i = aln->nsip[b];i--;){
                RUN(update_gaps(aln->sl[aln->sip[b][i]],aln->gaps[aln->sip[b][i]],path[0],gap_b));
        }
        MFREE(gap_a);
        MFREE(gap_b);
        return aln;
ERROR:
        return NULL;
}

int update_gaps(int old_len,int*gis,int new_len,int *newgaps)
{
        unsigned int i,j;
        int add = 0;
        int rel_pos = 0;
        for (i = 0; i <= old_len;i++){
                add = 0;
                for (j = rel_pos;j <= rel_pos + gis[i];j++){
                        if (newgaps[j] != 0){
                                add += newgaps[j];
                        }
                }
                rel_pos += gis[i]+1;
                gis[i] += add;
        }
        return OK;
}



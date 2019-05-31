
#include "counts_from_random_trees.h"

struct counts{
        double emit[26][26];
        double back[26];
        double MM;
        double GPO;
        double GPE;
        double TM;
        double TGPE;
        double tau;
        double eta;
        double num_alignments;
        double num_seq;
        double id_threshold;
        int L;
};
struct counts* init_counts(void);
int clean_counts(struct counts* aln_counts);
int normalize_counts(struct counts*aln_counts);
int fill_counts(struct counts* aln_counts, struct alignment* aln);

int pair_fill(struct counts* aln_counts, uint8_t*a,uint8_t*b,int len,double id_threshold);


int counts_from_random_trees(struct alignment* aln, struct aln_param* ap, int num_iter)
{

        int** map = NULL;       /* holds all alignment paths  */

        struct counts* aln_counts = NULL;
        int i,j,outer_iter;

        double sum;

        int mpos = 0;
        RUNP(aln_counts = init_counts());
        aln_counts->id_threshold = 0.1;
        int outer = 5;

        for(outer_iter = 0; outer_iter < outer; outer_iter++){
                clean_counts(aln_counts);
                for(i = 0; i < num_iter;i++){
                        LOG_MSG("Alignment %d",i);
                        random_tree(ap, aln->numseq);
                        RUNP(map = hirschberg_alignment(aln, ap));
                        RUN(weave(aln , map, ap->tree));

                        /* collect counts */

                        RUN(fill_counts(aln_counts, aln));
                        clean_aln(aln);

                        for(j = 0; j < aln->num_profiles ;j++){
                                if(map[j]){
                                        MFREE(map[j]);
                                }

                        }
                        MFREE(map);
                        map = NULL;


                }
                aln_counts->L = aln->L;
                /* normalize  */
                RUN(normalize_counts(aln_counts));
                /* add to ap */


                for(i = 0; i < aln_counts->L;i++){
                        //fprintf(stdout,"%d",i);
                        for(j = 0; j <= i;j++){
                                sum = log2(aln_counts->emit[i][j] / ( aln_counts->back[i] * aln_counts->back[j])) + log2(aln_counts->MM/((aln_counts->eta)*(aln_counts->eta)));
                                //fprintf(stdout,"old: %f new:%f\n", ap->subm[i][j], sum);
                                ap->subm[i][j] = sum;
                                ap->subm[j][i] = sum;
                                //ap->
                                //fprintf(stdout," %f,",sum);;
                        }
                        //fprintf(stdout,"\n");
                }


                sum = -1.0 * log2( (aln_counts->GPO * aln_counts->TM) / ((aln_counts->eta) * aln_counts->MM));
                ap->gpo = sum;
                //fprintf(stdout,"GPO:%f %f\n", sum,ap->gpo);
                sum = -1.0 *log2(aln_counts->GPE/(1.0 - aln_counts->tau));
                ap->gpe = sum;
                //fprintf(stdout,"GPE:%f %f\n", sum,ap->gpe);
                sum = -1.0 *log2(aln_counts->TGPE/(1.0 - aln_counts->tau));
                ap->tgpe = sum;
                //fprintf(stdout,"TGPE:%f %f\n", sum ,ap->tgpe);
        }
        MFREE(aln_counts);
        return OK;
ERROR:
        if(aln_counts){
                MFREE(aln_counts);
        }
        if(map){
                for(j = 0; j < aln->num_profiles ;j++){
                        if(map[j]){
                                MFREE(map[j]);
                        }

                }
                MFREE(map);
                map = NULL;

        }

        return FAIL;
}




int fill_counts(struct counts* aln_counts, struct alignment* aln)
{
        int i,j;
        int aln_len;

        uint8_t* aligned_a = NULL;
        uint8_t* aligned_b = NULL;

        aln_len = 0;
        for(i = 0; i <= aln->sl[0];i++){
                aln_len += aln->gaps[0][i];
        }
        aln_len += aln->sl[0];
        LOG_MSG("Aln len: %d.",aln_len);

        MMALLOC(aligned_a, sizeof(uint8_t) * aln_len);
        MMALLOC(aligned_b, sizeof(uint8_t) * aln_len);

        for(i = 0; i < aln->numseq;i++){
                aln_counts->num_seq++;
                aln_counts->eta += aln->sl[i];
                RUN(make_aliged_seq(aligned_a, aln->s[i], aln->gaps[i], aln->sl[i]));
                for(j = i+1; j < aln->numseq;j++){
                        //ASSERT(aln->sl[i] == aln->sl[j], "Sequences not aligned?");


                        RUN(make_aliged_seq(aligned_b, aln->s[j], aln->gaps[j], aln->sl[j]));

                        RUN(pair_fill(aln_counts, aligned_a, aligned_b, aln_len, aln_counts->id_threshold));
                }
        }
        MFREE(aligned_a);
        MFREE(aligned_b);
        return OK;
ERROR:
        return FAIL;
}


int pair_fill(struct counts* aln_counts, uint8_t*a,uint8_t*b,int len,double id_threshold)
{
        int i;
        int state = 1;          /* match - 0 is a galn_counts */
        int p_aln_len = 0;
        int begin = 0;
        int end = 0;
        int len_a;
        int len_b;
        double sim = 1.0;



        for(i = 0;i < len;i++){
                aln_counts->TGPE += sim;
                if(a[i] != 255 && b[i] != 255){

                        begin = i;
                        break;
                }
        }

        for(i = len-1; i>= 0;i--){
                aln_counts->TGPE += sim;
                if(a[i] != 255 && b[i] != 255){

                        end = i;
                        break;
                }

        }

        for(i = begin;i < end;i++){
                if(a[i] != 255 && b[i] != 255){
                        if(a[i] == b[i]){
                                sim += 1.0;
                        }
                }
        }
        len_a = 0;
        len_b = 0;
        for(i = 0;i < len;i++){

                if(a[i] != 255){
                        len_a++;
                }
                if(b[i]  != 255){
                        len_b++;
                }
        }

        sim =   sim / (double) (MACRO_MIN(len_a, len_b));
        //sim = sim * sim * sim * sim;
        //sim = 1.0;
        //fprintf(stdout,"Sim :%f\n", sim);

        if(sim > id_threshold){
                sim = 1.0;
        }else{

                return OK;
        }


        for(i = begin;i < end;i++){
                if(a[i] == 255 && b[i] == 255){ /* two galn_countss do nothing */

                }else if(a[i] != 255 && b[i] != 255){ /* aligned residues  */

                        if(state == 0){
                                aln_counts->TM += sim;
                        }else{
                                aln_counts->MM += sim;
                        }
                        if(a[i] == b[i]){
                                aln_counts->emit[a[i]][b[i]] += sim + sim;
                        }else if(a[i] > b[i]){
                                aln_counts->emit[a[i]][b[i]] += sim;
                        }else{
                                aln_counts->emit[b[i]][a[i]] += sim;
                        }

                        aln_counts->back[a[i]] += sim;
                        aln_counts->back[b[i]] += sim;
                        p_aln_len++;
                        state = 1;
                }else if(a[i] != 255 && b[i] == 255){
                        if(state == 0){
                                aln_counts->GPE += sim;
                        }else{
                                aln_counts->GPO += sim;
                        }
                        aln_counts->back[a[i]]++;
                        p_aln_len++;

                        state = 0;

                }else if(a[i] == 255 && b[i] != 255){
                        if(state == 0){
                                aln_counts->GPE += sim;
                        }else{
                                aln_counts->GPO += sim;
                        }
                        aln_counts->back[b[i]]+= sim;
                        p_aln_len++;
                        state = 0;
                }
        }
        aln_counts->num_alignments++;
        aln_counts->tau += p_aln_len;
        return OK;

}


struct counts* init_counts(void)
{

        struct counts* aln_counts = NULL;
        MMALLOC(aln_counts, sizeof(struct counts));


        RUN(clean_counts(aln_counts));
        return aln_counts;
ERROR:
        return NULL;
}


int clean_counts(struct counts* aln_counts)
{
        int i,j;

        double pseudocount = 1.0;
        for(i = 0; i < 26;i++){

                aln_counts->back[i] = pseudocount;
                for(j = 0; j < 26;j++){
                        aln_counts->emit[i][j] = pseudocount;
                }
        }
        aln_counts->MM = pseudocount;
        aln_counts->TM = pseudocount;
        aln_counts->TGPE = pseudocount;
        aln_counts->GPE = pseudocount;
        aln_counts->GPO = pseudocount;
        aln_counts->tau = 0.0;
        aln_counts->eta = 0.0;
        aln_counts->num_seq = 0.0;
        aln_counts->num_alignments = 0.0;
        return OK;
}


int normalize_counts(struct counts*aln_counts)
{
        int i,j;
        double sum;
        //double tmp;

        sum = 0.0;

        for(i = 0; i < aln_counts->L;i++){
                //fprintf(stdout,"BAK:%d %f\n",i,aln_counts->back[i]);
                sum += aln_counts->back[i];
        }
        for(i = 0; i < aln_counts->L;i++){
                aln_counts->back[i] /= sum;
                //fprintf(stdout,"BAK:%d %f\n",i,aln_counts->back[i]);
        }
        sum = 0.0;
        for(i = 0; i < aln_counts->L;i++){
                for(j = 0; j <= i;j++){
                        sum += aln_counts->emit[i][j];
                }
        }
        for(i = 0; i < aln_counts->L;i++){
                for(j = 0; j <= i;j++){
                        aln_counts->emit[i][j]/= sum;
                        aln_counts->emit[j][i] = aln_counts->emit[i][j];
                }
        }

        aln_counts->tau = 1.0 / (aln_counts->tau / aln_counts->num_alignments);

        sum = aln_counts->TGPE;
        sum += aln_counts->num_alignments*2.0; /* we need to get in an out of the alignment */

        aln_counts->TGPE = aln_counts->TGPE / sum;


        //tmp = (aln_counts->num_alignments*2.0) / sum;
        //fprintf(stdout,"%f %f ", aln_counts->TGPE, tmp);
        //aln_counts->TGPE = -log2 (1.0-aln_counts->TGPE);
        //fprintf(stdout,"-> %f\n", aln_counts->TGPE);
        sum = aln_counts->MM;
        sum += aln_counts->GPO;
        sum += aln_counts->GPO;

        aln_counts->MM = aln_counts->MM / sum;
        aln_counts->GPO = aln_counts->GPO / sum;


        /* add in tau */
        sum = aln_counts->MM;
        sum += aln_counts->GPO;
        sum += aln_counts->GPO;
        sum += aln_counts->tau;

        aln_counts->MM = aln_counts->MM / sum;
        aln_counts->GPO = aln_counts->GPO / sum;


        sum = aln_counts->GPE;
        sum += aln_counts->TM;

        aln_counts->GPE = aln_counts->GPE / sum;
        aln_counts->TM = aln_counts->TM / sum;

        sum = aln_counts->GPE;
        sum += aln_counts->TM;
        sum += aln_counts->tau;
        aln_counts->GPE = aln_counts->GPE / sum;
        aln_counts->TM = aln_counts->TM / sum;

        /* eta  */
        aln_counts->eta = 1.0 - 1.0 / (aln_counts->eta / aln_counts->num_seq);
        return OK;
}

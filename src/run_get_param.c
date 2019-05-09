

#include "global.h"
#include <getopt.h>
#include "align_io.h"

#include "alphabet.h"

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
        int L;
};
struct counts* init_counts(void);
int normalize_counts(struct counts*ap);
int print_counts(struct counts* ap);
int fill_counts(struct counts* ap, struct alignment* aln);
int pair_fill(struct counts* ap, uint8_t*a,uint8_t*b,int len);

int print_probabilies(struct counts*ap);


int main(int argc, char *argv[])
{
        struct alignment* aln = NULL;

        struct counts* ap = NULL;

        int i;
        int num_infiles = 0;
        char** infile = NULL;
        int c = 0;
        //int help = 0;
        while (1){
                static struct option long_options[] ={
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"h",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 'h':
                        //help = 1;
                        break;
                default:
                        ERROR_MSG("not recognized");
                        break;
                }
        }


        if (optind < argc){

                //fprintf(stderr,"EXTRA :%d\n",argc - optind);
                num_infiles = argc-optind;
                MMALLOC(infile, sizeof(char*) * num_infiles);
                c = 0;
                while (optind < argc){
                        infile[c] =  argv[optind++];
                        c++;
                }
        }

        RUNP(ap = init_counts());

        for(i = 0; i < num_infiles;i++){
                fprintf(stdout,"%s\n",infile[i]);
                RUNP(aln = read_alignment(infile[i]));

                RUN(convert_alignment_to_internal(aln, defPROTEIN));
                fprintf(stdout,"%d L \n", aln->L);
                RUN(fill_counts(ap, aln));
                ap->L = aln->L;
                free_aln(aln);
        }
        normalize_counts(ap);
        print_counts(ap);
        print_probabilies(ap);

        if(infile){
                MFREE(infile);
        }
        MFREE(ap);
        return EXIT_SUCCESS;
ERROR:

        if(infile){
                MFREE(infile);
        }

        return EXIT_FAILURE;
}

int fill_counts(struct counts* ap, struct alignment* aln)
{
        int i,j;

        for(i = 0; i < aln->numseq;i++){
                ap->num_seq++;
                ap->eta += aln->sl[i];
                for(j = i+1; j < aln->numseq;j++){
                        ASSERT(aln->sl[i] == aln->sl[j], "Sequences not aligned?");
                        RUN(pair_fill(ap, aln->s[i], aln->s[j], aln->sl[i]));
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int pair_fill(struct counts* ap, uint8_t*a,uint8_t*b,int len)
{
        int i;
        int state = 1;          /* match - 0 is a gap */
        int p_aln_len = 0;
        int begin = 0;
        int end = 0;

        double sim = 0;
        for(i = 0;i < len;i++){
                ap->TGPE += sim;
                if(a[i] != 255 && b[i] != 255){

                        begin = i;
                        break;
                }
        }

        for(i = len-1; i>= 0;i--){
                ap->TGPE += sim;
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

        sim =   sim / (double) (end - begin);
        sim = sim * sim * sim * sim;
        //sim = 1.0;
        fprintf(stdout,"Sim :%f\n", sim);


        //sim = 1.0;


        for(i = begin;i < end;i++){
                if(a[i] == 255 && b[i] == 255){ /* two gaps do nothing */

                }else if(a[i] != 255 && b[i] != 255){ /* aligned residues  */

                        if(state == 0){
                                ap->TM += sim;
                        }else{
                                ap->MM += sim;
                        }
                        if(a[i] == b[i]){
                                ap->emit[a[i]][b[i]] += sim + sim;
                        }else if(a[i] > b[i]){
                                ap->emit[a[i]][b[i]] += sim;
                        }else{
                                ap->emit[b[i]][a[i]] += sim;
                        }

                        ap->back[a[i]] += sim;
                        ap->back[b[i]] += sim;
                        p_aln_len++;
                        state = 1;
                }else if(a[i] != 255 && b[i] == 255){
                        if(state == 0){
                                ap->GPE += sim;
                        }else{
                                ap->GPO += sim;
                        }
                        ap->back[a[i]]++;
                        p_aln_len++;

                        state = 0;

                }else if(a[i] == 255 && b[i] != 255){
                        if(state == 0){
                                ap->GPE += sim;
                        }else{
                                ap->GPO += sim;
                        }
                        ap->back[b[i]]+= sim;
                        p_aln_len++;
                        state = 0;
                }
        }
        ap->num_alignments++;
        ap->tau += p_aln_len;
        return OK;

}


struct counts* init_counts(void)
{
        int i,j;
        struct counts* ap = NULL;
        MMALLOC(ap, sizeof(struct counts));

        for(i = 0; i < 26;i++){
                ap->back[i] = 1.0;
                for(j = 0; j < 26;j++){
                        ap->emit[i][j] = 1.0;
                }
        }
        ap->MM = 0.0;
        ap->TM = 0.0;
        ap->TGPE = 0.0;
        ap->GPE = 0.0;
        ap->GPO = 0.0;
        ap->tau = 0.0;
        ap->eta = 0.0;
        ap->num_seq = 0.0;
        ap->num_alignments = 0.0;
        return ap;
ERROR:
        return NULL;
}


int normalize_counts(struct counts*ap)
{
        int i,j;
        double sum;
        double tmp;

        sum = 0.0;

        for(i = 0; i < ap->L;i++){
                sum += ap->back[i];
        }
        for(i = 0; i < ap->L;i++){
                ap->back[i] /= sum;
        }
        sum = 0.0;
        for(i = 0; i < ap->L;i++){
                for(j = 0; j <= i;j++){
                        sum += ap->emit[i][j];
                }
        }
        for(i = 0; i < ap->L;i++){
                for(j = 0; j <= i;j++){
                        ap->emit[i][j]/= sum;
                        ap->emit[j][i] = ap->emit[i][j];
                }
        }

        ap->tau = 1.0 / (ap->tau / ap->num_alignments);

        sum = ap->TGPE;
        sum += ap->num_alignments*2.0; /* we need to get in an out of the alignment */

        ap->TGPE = ap->TGPE / sum;

        tmp = (ap->num_alignments*2.0) / sum;
        fprintf(stdout,"%f %f ", ap->TGPE, tmp);
        ap->TGPE = -log2 (1.0-ap->TGPE);
        fprintf(stdout,"-> %f\n", ap->TGPE);
        sum = ap->MM;
        sum += ap->GPO;
        sum += ap->GPO;

        ap->MM = ap->MM / sum;
        ap->GPO = ap->GPO / sum;


        /* add in tau */
        sum = ap->MM;
        sum += ap->GPO;
        sum += ap->GPO;
        sum += ap->tau;

        ap->MM = ap->MM / sum;
        ap->GPO = ap->GPO / sum;


        sum = ap->GPE;
        sum += ap->TM;

        ap->GPE = ap->GPE / sum;
        ap->TM = ap->TM / sum;

        sum = ap->GPE;
        sum += ap->TM;
        sum += ap->tau;
        ap->GPE = ap->GPE / sum;
        ap->TM = ap->TM / sum;

        /* eta  */

        return OK;
}


int print_probabilies(struct counts*ap)
{

        int i,j;
        fprintf(stdout,"float prior_back[%d] = {\n",ap->L);
        for(i = 0; i < ap->L;i++){
                fprintf(stdout,"%f,\n", prob2scaledprob(ap->back[i]));
        }
        fprintf(stdout,"%f};\n", prob2scaledprob(ap->back[ap->L]));

        fprintf(stdout,"float prior_m[%d][%d] = {\n", ap->L,ap->L);
        for(i = 0; i < ap->L;i++){
                j = 0;
                fprintf(stdout,"{%f", prob2scaledprob(ap->emit[i][j]));
                for(j = 1; j < ap->L;j++){
                        fprintf(stdout,",%f", prob2scaledprob(ap->emit[i][j]));
                }
                fprintf(stdout,"}");
                if(i != ap->L){
                        fprintf(stdout,",");
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"};\n");

        fprintf(stdout,"float prior_MM = %f;\n",prob2scaledprob(ap->MM));
        fprintf(stdout,"float prior_GPO = %f;\n",prob2scaledprob(ap->GPO));
        fprintf(stdout,"float prior_GPE = %f;\n",prob2scaledprob(ap->GPE));
        fprintf(stdout,"float prior_TM = %f;\n",prob2scaledprob(ap->TM ));


        return OK;
}

int print_counts(struct counts* ap)
{
        int i,j;

        double sum = 0.0;

        ap->eta = 1.0 - 1.0 / (ap->eta / ap->num_seq);
        /*for(i = 0; i < 26;i++){
                fprintf(stdout,"%*d ",3,i);
        }
        fprintf(stdout,"\n");
        */
        /*for(i = 0; i < 26;i++){
                fprintf(stdout,"%*.2f ",3, ap->back[i]);
        }
        fprintf(stdout,"\n");
        fprintf(stdout,"\n");
        */
        fprintf(stdout,"float balimt[]={\n");

        for(i = 0; i < ap->L;i++){
                //fprintf(stdout,"%d",i);
                for(j = 0; j <= i;j++){
                        sum = log2(ap->emit[i][j] / ( ap->back[i] * ap->back[j])) + log2(ap->MM/((ap->eta)*(ap->eta)));
                        fprintf(stdout," %f,",sum);;
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"};\n");
        /*fprintf(stdout,"%f\tMM\n", ap->MM);
        fprintf(stdout,"%f\tGPO\n", ap->GPO);
        fprintf(stdout,"%f\tGPE\n", ap->GPE);
        fprintf(stdout,"%f\tTM\n", ap->TM);*/
        /* taushould be 1/ average length */
        fprintf(stdout,"%f\ttau\n", ap->tau);
        fprintf(stdout,"%f\teta\n", ap->eta);
        sum = 0.0;

        sum = -1.0 * log2( (ap->GPO * ap->TM) / ((ap->eta) * ap->MM));
        fprintf(stdout,"ap->gpo = %f;\n", sum);
        sum = -1.0 *log2(ap->GPE/(1.0 - ap->tau));
        fprintf(stdout,"ap->gpe =  %f;\n", sum);
        fprintf(stdout,"ap->tgpe =  %f;\n", ap->TGPE);
        //;
        //sum = ap->MM / ((ap->eta)*(ap->eta));
        //fprintf(stdout,"%f\n",sum);
        return OK;
}

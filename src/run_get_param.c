
#include "global.h"
#include <getopt.h>
#include "align_io.h"

struct counts{
        double emit[26][26];
        double back[26];
        double MM;
        double GPO;
        double GPE;
        double TM;
        double tau;
        double num_alignments;
};
struct counts* init_counts(void);
int print_counts(struct counts* ap);
int fill_counts(struct counts* ap, struct alignment* aln);
int pair_fill(struct counts* ap, int*a,int*b,int len);

int main(int argc, char *argv[])
{
        struct alignment* aln = NULL;

        struct counts* ap = NULL;

        int i;
        int num_infiles;
        char** infile = NULL;
        int c = 0;
        int help = 0;
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
                        help = 1;
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

                RUN(fill_counts(ap, aln));
                free_aln(aln);
        }
        print_counts(ap);

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
                for(j = i+1; j < aln->numseq;j++){
                        ASSERT(aln->sl[i] == aln->sl[j], "Sequences not aligned?");
                        RUN(pair_fill(ap, aln->s[i], aln->s[j], aln->sl[i]));
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int pair_fill(struct counts* ap, int*a,int*b,int len)
{
        int i;
        int state = 1;          /* match - 0 is a gap */
        int p_aln_len = 0;
        for(i = 0;i < len;i++){
                if(a[i] == -1 && b[i] == -1){ /* two gaps do nothing */

                }else if(a[i] != -1 && b[i] != -1){ /* aligned residues  */
                        if(state == 0){
                                ap->TM++;
                        }else{
                                ap->MM++;
                        }
                        if(a[i] < b[i]){
                                ap->emit[a[i]][b[i]]++;
                        }else{
                                ap->emit[b[i]][a[i]]++;
                        }

                        ap->back[a[i]]++;
                        ap->back[b[i]]++;
                        p_aln_len++;
                        state = 1;
                }else if(a[i] != -1 && b[i] == -1){
                        if(state == 0){
                                ap->GPE++;
                        }else{
                                ap->GPO++;
                        }
                        ap->back[a[i]]++;
                        p_aln_len++;

                        state = 0;

                }else if(a[i] == -1 && b[i] != -1){
                        if(state == 0){
                                ap->GPE++;
                        }else{
                                ap->GPO++;
                        }
                        ap->back[b[i]]++;
                        p_aln_len++;
                        state = 0;
                }
        }
        return OK;

}


struct counts* init_counts(void)
{
        int i,j;
        struct counts* ap = NULL;
        MMALLOC(ap, sizeof(struct counts));

        for(i = 0; i < 26;i++){
                ap->back[i] = 0.0;
                for(j = 0; j < 26;j++){
                        ap->emit[i][j] = 0.0;
                }
        }
        ap->MM = 0.0;
        ap->TM = 0.0;
        ap->GPE = 0.0;
        ap->GPO = 0.0;
        ap->tau = 0.0;
        ap->num_alignments = 0.0;
        return ap;
ERROR:
        return NULL;
}

int print_counts(struct counts* ap)
{
        int i,j;


        for(i = 0; i < 26;i++){
                fprintf(stdout,"%*d ",3,i);
        }
        fprintf(stdout,"\n");

        for(i = 0; i < 26;i++){
                fprintf(stdout,"%*.0f ",3, ap->back[i]);
        }
        fprintf(stdout,"\n");
        fprintf(stdout,"\n");


        for(i = 0; i < 26;i++){
                fprintf(stdout,"%d",i);
                for(j = 0; j < 26;j++){
                        fprintf(stdout," %*.0f",3,ap->emit[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        return OK;
}

/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019 Timo Lassmann

    This file is part of kalign.

    Kalign is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/




#include "global.h"
#include "msa.h"
#include <getopt.h>
//#include "align_io.h"

#include "alphabet.h"
#include "kmeans.h"


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
int clean_counts(struct counts* ap);
int normalize_counts(struct counts*ap);
int print_counts(struct counts* ap);
int print_counts_to_file(struct counts* ap,char* out);
int print_counts_flat(struct counts* ap,double* raw);
int fill_counts(struct counts* ap, struct msa* msa);
int make_aliged_seq(uint8_t* aligned, uint8_t* unaligned, int* gaps,int len);
int pair_fill(struct counts* ap, uint8_t*a,uint8_t*b,int len,double id_threshold);
int print_probabilies(struct counts*ap);



int main(int argc, char *argv[])
{
        //struct alignment* aln = NULL;
        struct msa* msa = NULL;

        struct counts* ap = NULL;

        int* assignment = NULL;

        double** raw = NULL;
        double id_threshold;
        int i,j;
        int num_infiles = 0;
        char** infile = NULL;
        char* outfile = NULL;
        int c = 0;

        int print_all = 0;
        //int help = 0;
        id_threshold = 0.0;
        while (1){
                static struct option long_options[] ={
                        {"thres",  required_argument, 0, 't'},
                        {"out",  required_argument, 0, 'o'},
                        {"all",0,0,'a'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"hao:t:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 't':
                        id_threshold = atof(optarg);
                        break;
                case 'o':
                        outfile = optarg;
                        break;

                case 'a':
                        print_all = 1;
                        break;

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
        ap->id_threshold = id_threshold;
        if(print_all){
                /* get L */
                RUNP(msa = read_input(infile[0],msa));
                //RUNP(aln = read_alignment(infile[0]));

                RUN(convert_msa_to_internal(msa, defPROTEIN));
                //RUN(convert_alignment_to_internal(aln, defPROTEIN));
                ap->L = msa->L;
                int len = (ap->L *( ap->L-1)) / 2 + ap->L + 3;

                free_msa(msa);



                MMALLOC(assignment, sizeof(int) * num_infiles);
                for(i = 0; i < num_infiles;i++){
                        assignment[i] = 0;
                }

                //fprintf(stdout,"LEN:%d\n",ap->L);

                raw = galloc(raw,num_infiles,len,0.0);
                for(i = 0; i < num_infiles;i++){
                        clean_counts(ap);
                        //fprintf(stdout,"%s\n",infile[i]);
                        RUNP(msa = read_input(infile[0],msa));
                        //RUNP(aln = read_alignment(infile[0]));

                        RUN(convert_msa_to_internal(msa, defPROTEIN));

                        //RUNP(aln = read_alignment(infile[i]));

                        //RUN(convert_alignment_to_internal(aln, defPROTEIN));
                        //fprintf(stdout,"%d L \n", aln->L);
                        RUN(fill_counts(ap, msa));
                        //ap->L = aln->L;
                        normalize_counts(ap);
                        print_counts_flat(ap,raw[i]);
                        free_msa(msa);
                }
                /*
                for(i = 0; i < num_infiles;i++){
                        for(j = 0; j < len;j++){
                                fprintf(stdout,"%f,",raw[i][j]);
                        }
                        fprintf(stdout,"\n");
                }
                fprintf(stdout,"\n");
                */


                double** means = NULL;
                int k =4;
                means = kmeans(raw, assignment,num_infiles,len, k);
                //means = galloc(means,1,len,0.0);
                for(j = 0; j < k;j++){
                        clean_counts(ap);
                        for(i = 0; i < num_infiles;i++){
                                if(assignment[i] == j){
                                        //                fprintf(stdout,"%s %d\n",infile[i],assignment[i]);
                                        //RUNP(aln = read_alignment(infile[i]));
                                        RUNP(msa = read_input(infile[i],msa));
                                        RUN(convert_msa_to_internal(msa, defPROTEIN));

                                        //RUN(convert_alignment_to_internal(aln, defPROTEIN));
                                        //fprintf(stdout,"%d L \n", aln->L);
                                        RUN(fill_counts(ap, msa));
                                        //ap->L = aln->L;
                                        free_msa(msa);

                                        //free_aln(aln);
                                }


                        }

                        normalize_counts(ap);
                        print_counts_flat(ap,means[j]);
                }

                fprintf(stdout,"MEANs:\n");
                fprintf(stdout,"double defprot_set10[%d][%d] = {",k, len);
                for(i = 0; i < k;i++){
                        //means[i][len-2] /= 2.0;
                        fprintf(stdout, "{ %f ", means[i][0]);
                        for(j = 1; j < len;j++){
                                fprintf(stdout,",%f",means[i][j]);
                        }
                        if(i == k-1){
                                fprintf(stdout,"}\n");
                        }else{
                                fprintf(stdout,"},\n");
                        }
                }
                fprintf(stdout,"};\n");
        }else{
                for(i = 0; i < num_infiles;i++){
                        fprintf(stdout,"%s\n",infile[i]);
                        RUNP(msa = read_input(infile[i],msa));
                        //RUNP(aln = read_alignment(infile[i]));
                        //dealign(aln);
                        RUN(convert_msa_to_internal(msa, defPROTEIN));
                        //RUN(convert_alignment_to_internal(aln, defPROTEIN));
                        //fprintf(stdout,"%d L \n", aln->L);
                        RUN(fill_counts(ap, msa));
                        ap->L = msa->L;
                        free_msa(msa);
                        //free_aln(aln);
                }
                normalize_counts(ap);
                if(outfile){
                        print_counts_to_file(ap, outfile);
                }
                print_counts(ap);
                print_probabilies(ap);
        }
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

int fill_counts(struct counts* ap, struct msa* msa)
{
        int i,j;
        int aln_len;

        uint8_t* aligned_a = NULL;
        uint8_t* aligned_b = NULL;

        aln_len = 0;
        for(i = 0; i <= msa->sequences[0]->len;i++){
                aln_len += msa->sequences[0]->gaps[i];//  aln->gaps[0][i];
        }
        aln_len += msa->sequences[0]->len;//  aln->sl[0];
        LOG_MSG("Aln len: %d.",aln_len);

        MMALLOC(aligned_a, sizeof(uint8_t) * aln_len);
        MMALLOC(aligned_b, sizeof(uint8_t) * aln_len);

        for(i = 0; i < msa->numseq;i++){
                ap->num_seq++;
                ap->eta +=  msa->sequences[i]->len;//  aln->sl[i];
                RUN(make_aliged_seq(aligned_a, msa->sequences[i]->s, msa->sequences[i]->gaps, msa->sequences[i]->len));

                //RUN(make_aliged_seq(aligned_a, aln->s[i], aln->gaps[i], aln->sl[i]));
                for(j = i+1; j < msa->numseq;j++){
                        //ASSERT(aln->sl[i] == aln->sl[j], "Sequences not aligned?");


                        RUN(make_aliged_seq(aligned_b, msa->sequences[j]->s, msa->sequences[j]->gaps, msa->sequences[j]->len));
                        //RUN(make_aliged_seq(aligned_b, aln->s[j], aln->gaps[j], aln->sl[j]));

                        RUN(pair_fill(ap, aligned_a, aligned_b, aln_len, ap->id_threshold));
                }
        }
        MFREE(aligned_a);
        MFREE(aligned_b);
        return OK;
ERROR:
        return FAIL;
}


int pair_fill(struct counts* ap, uint8_t*a,uint8_t*b,int len,double id_threshold)
{
        int i;
        int state = 1;          /* match - 0 is a gap */
        int p_aln_len = 0;
        int begin = 0;
        int end = 0;
        int len_a;
        int len_b;
        double sim = 1.0;

        double tmp_tgpe = 0.0;


        for(i = 0;i < len;i++){
                tmp_tgpe += sim;
                //ap->TGPE += sim;
                if(a[i] != 255 && b[i] != 255){

                        begin = i;
                        break;
                }
        }

        for(i = len-1; i>= 0;i--){
                tmp_tgpe += sim;

                //ap->TGPE += sim;
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
        fprintf(stdout,"Sim :%f\n", sim);

        if(sim > id_threshold){
                ap->TGPE += tmp_tgpe;
                sim = 1.0;
        }else{

                return OK;
        }


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

        struct counts* ap = NULL;
        MMALLOC(ap, sizeof(struct counts));


        RUN(clean_counts(ap));
        return ap;
ERROR:
        return NULL;
}


int clean_counts(struct counts* ap)
{
        int i,j;

        double pseudocount = 1.0;
        for(i = 0; i < 26;i++){

                ap->back[i] = pseudocount;
                for(j = 0; j < 26;j++){
                        ap->emit[i][j] = pseudocount;
                }
        }
        ap->MM = pseudocount;
        ap->TM = pseudocount;
        ap->TGPE = pseudocount;
        ap->GPE = pseudocount;
        ap->GPO = pseudocount;
        ap->tau = 0.0;
        ap->eta = 0.0;
        ap->num_seq = 0.0;
        ap->num_alignments = 0.0;
        return OK;
}


int normalize_counts(struct counts*ap)
{
        int i,j;
        double sum;
        //double tmp;

        sum = 0.0;

        for(i = 0; i < ap->L;i++){
                //fprintf(stdout,"BAK:%d %f\n",i,ap->back[i]);
                sum += ap->back[i];
        }
        for(i = 0; i < ap->L;i++){
                ap->back[i] /= sum;
                //fprintf(stdout,"BAK:%d %f\n",i,ap->back[i]);
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


        //tmp = (ap->num_alignments*2.0) / sum;
        //fprintf(stdout,"%f %f ", ap->TGPE, tmp);
        //ap->TGPE = -log2 (1.0-ap->TGPE);
        //fprintf(stdout,"-> %f\n", ap->TGPE);
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
        ap->eta = 1.0 - 1.0 / (ap->eta / ap->num_seq);
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


int print_counts_flat(struct counts* ap,double* raw)
{
        int i,j;

        double sum = 0.0;
        int c = 0;
        for(i = 0; i < ap->L;i++){
                //fprintf(stdout,"%d",i);
                for(j = 0; j <= i;j++){

                        //fprintf(stdout,"\n%d %f\n",i,ap->back[i]);
                        sum = log2(ap->emit[i][j] / ( ap->back[i] * ap->back[j])) + log2(ap->MM/((ap->eta)*(ap->eta)));
                        raw[c] = sum;
                        //                    fprintf(stdout,"%f,",sum);
                        c++;
                }
        }

        sum = -1.0 * log2( (ap->GPO * ap->TM) / ((ap->eta) * ap->MM));
        //fprintf(stdout,"%f,", sum);
        raw[c] = sum /2.0;
        c++;
        sum = -1.0 *log2(ap->GPE/(1.0 - ap->tau));
        //fprintf(stdout,"%f   %d\n", sum,c);
        raw[c] = sum;
        c++;

        sum = -1.0 *log2(ap->TGPE/(1.0 - ap->tau));
        raw[c] = sum;
        c++;

        return OK;
}

int print_counts_to_file(struct counts* ap,char* out)
{
        FILE* f_ptr = NULL;
        int i,j;
        double sum;
        LOG_MSG("Open %s",out);
        RUNP( f_ptr = fopen(out, "w"));
        for(i = 0; i < ap->L;i++){
                //fprintf(stdout,"%d",i);
                for(j = 0; j <= i;j++){
                        sum = log2(ap->emit[i][j] / ( ap->back[i] * ap->back[j])) + log2(ap->MM/((ap->eta)*(ap->eta)));
                        fprintf(f_ptr,"%f\n",sum);;
                }

        }
        sum = 0.0;

        sum = -1.0 * log2( (ap->GPO * ap->TM) / ((ap->eta) * ap->MM));
        fprintf(f_ptr,"%f\n", sum/ 2.0);
        sum = -1.0 *log2(ap->GPE/(1.0 - ap->tau));
        fprintf(f_ptr,"%f\n", sum);
        sum = -1.0 *log2(ap->TGPE/(1.0 - ap->tau));
        fprintf(f_ptr,"%f\n", sum);

        fclose(f_ptr);
        return OK;
ERROR:
        return FAIL;
}

int print_counts(struct counts* ap)
{
        int i,j;

        double sum = 0.0;

        //ap->eta = 1.0 - 1.0 / (ap->eta / ap->num_seq);
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
        fprintf(stdout,"%f\tTGPE\n", ap->TGPE);
        fprintf(stdout,"%f\tTM\n", ap->TM);*/
        /* taushould be 1/ average length */
        //fprintf(stdout,"%f\ttau\n", ap->tau);
        //fprintf(stdout,"%f\teta\n", ap->eta);
        sum = 0.0;

        sum = -1.0 * log2( (ap->GPO * ap->TM) / ((ap->eta) * ap->MM));
        fprintf(stdout,"ap->gpo = %f;\n", sum/ 2.0);
        sum = -1.0 *log2(ap->GPE/(1.0 - ap->tau));
        fprintf(stdout,"ap->gpe =  %f;\n", sum);
        sum = -1.0 *log2(ap->TGPE/(1.0 - ap->tau));
        fprintf(stdout,"ap->tgpe =  %f;\n", sum);
        //;
        //sum = ap->MM / ((ap->eta)*(ap->eta));
        //fprintf(stdout,"%f\n",sum);
        return OK;
}


int make_aliged_seq(uint8_t* aligned, uint8_t* unaligned, int* gaps,int len)
{
        int c;
        int i;
        int tmp;
        c= 0;
        for(i = 0; i < len;i++){
                tmp = gaps[i];
                while(tmp){
                        aligned[c] = 255;
                        tmp--;
                        c++;
                }
                aligned[c] = unaligned[i];
                c++;

        }
        tmp = gaps[len];
        while(tmp){
                aligned[c] = 255;
                tmp--;
                c++;
        }
        return OK;
}

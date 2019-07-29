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
#include "alphabet.h"


/* only local; */
struct line_buffer{
        struct out_line** lines;
        int max_line_len;
        int alloc_num_lines;
        int num_line;
};

struct out_line{
        char* line;
        int block;
        int seq_id;
};




struct msa* read_fasta(char* infile, struct msa* msa);
struct msa* read_msf(char* infile, struct msa* msa);
struct msa* read_clu(char* infile, struct msa* msa);

int write_msa_fasta(struct msa* msa,char* outfile);
int write_msa_clustal(struct msa* msa,char* outfile);
int write_msa_msf(struct msa* msa,char* outfile);

/* memory functions  */
struct msa* alloc_msa(void);
int resize_msa(struct msa* msa);


struct msa_seq* alloc_msa_seq(void);
int resize_msa_seq(struct msa_seq* seq);
void free_msa_seq(struct msa_seq* seq);



struct line_buffer* alloc_line_buffer(int max_line_len);
int resize_line_buffer(struct line_buffer* lb);
void free_line_buffer(struct line_buffer* lb);

/* local helper functions  */

static int detect_alignment_format(char* infile,int* type);
static int detect_aligned(struct msa* msa);
static int detect_alphabet(struct msa* msa);
static int set_sip_nsip(struct msa* msa);
static int null_terminate_sequences(struct msa* msa);
static int sort_out_lines(const void *a, const void *b);
static int make_linear_sequence(struct msa_seq* seq, char* linear_seq);
int GCGMultchecksum(struct msa* msa);
/* Taken from squid library by Sean Eddy  */
int GCGchecksum(char *seq, int len);


#ifdef RWALIGN_TEST
int print_msa(struct msa* msa);
int main(int argc, char *argv[])
{
        struct msa* msa = NULL;
        char buffer[BUFFER_LEN];
        LOG_MSG("Start io tests.");

        char* datadir;

        datadir = getenv("testdatafiledir");

        snprintf(buffer,BUFFER_LEN, "%s/%s", datadir, argv[1]);
        LOG_MSG("reading: %s", buffer);
        RUNP(msa = read_input(buffer,NULL));
//print_msa(msa);
        write_msa_clustal(msa,"rwtest.clu");

        write_msa_fasta(msa, "rwtest.fasta");

        write_msa_msf(msa,"rwtest.msf");
        free_msa(msa);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int print_msa(struct msa* msa)
{
        int i;
        int j;
        int c;
        LOG_MSG("MSA:");
        for(i = 0; i < msa->numseq;i++){
                fprintf(stdout,"%s\n", msa->sequences[i]->name);
                for(j = 0;j < msa->sequences[i]->len;j++){
                        for(c = 0;c < msa->sequences[i]->gaps[j];c++){
                                fprintf(stdout,"-");
                        }
                        fprintf(stdout,"%c", msa->sequences[i]->seq[j]);
                }
                for(c = 0;c < msa->sequences[i]->gaps[ msa->sequences[i]->len];c++){
                        fprintf(stdout,"-");
                }
                fprintf(stdout,"\n");
        }
        return OK;
}
#endif

struct msa* read_input(char* infile,struct msa* msa)
{

        int type;

        ASSERT(infile != NULL,"No input file");
        /* sanity checks  */
        if(!my_file_exists(infile)){
                ERROR_MSG("File: %s does not exist.",infile);
        }

        DECLARE_TIMER(timer);

        START_TIMER(timer);
        RUN(detect_alignment_format(infile, &type));


        if(type == FORMAT_FA){
                RUNP(msa = read_fasta(infile,msa));
        }else if(type == FORMAT_MSF){
                RUNP(msa = read_msf(infile,msa));
        }else if(type == FORMAT_CLU){
                RUNP(msa = read_clu(infile,msa));
        }

        RUN(detect_alphabet(msa));
        RUN(detect_aligned(msa));

        RUN(set_sip_nsip(msa));

        STOP_TIMER(timer);
        LOG_MSG("Done reading input sequences in %f seconds.", GET_TIMING(timer));
        return msa;
ERROR:
        return NULL;
}

int write_msa(struct msa* msa, char* outfile, int type)
{

        ASSERT(msa!= NULL, "No alignment");

        if(type == FORMAT_FA){
                RUN(write_msa_fasta(msa, outfile));
        }else if(type == FORMAT_MSF){
                RUN(write_msa_msf(msa, outfile));
        }else if(type == FORMAT_CLU){
                RUN(write_msa_clustal(msa, outfile));
        }else{
                ERROR_MSG("Output format not recognized.");
        }

        return OK;
ERROR:
        return FAIL;
}

/* detect alignment format  */

int detect_alignment_format(char* infile,int* type)
{
        FILE* f_ptr = NULL;
        char line[BUFFER_LEN];
        int hints[3];
        int line_len;
        int line_number;
        int set;
        int i;
        ASSERT(infile != NULL,"No input file");
        /* sanity checks  */
        if(!my_file_exists(infile)){
                ERROR_MSG("File: %s does not exist.",infile);
        }

        line_number = 0;
        for(i = 0; i < 3; i++){
                hints[i] =0;
        }
        RUNP(f_ptr = fopen(infile, "r"));

        /* scan through first line header  */
        while(fgets(line, BUFFER_LEN, f_ptr)){
                line_len = strnlen(line, BUFFER_LEN);
                line[line_len-1] = 0;

                line_len--;
                if(line[0] == '>'){
                        hints[0]++; /* fasta */
                }

                if(strstr(line, "multiple sequence alignment")){
                        hints[2]++; /* clustal format  */
                }
                if(strstr(line, "CLUSTAL W")){
                        hints[2]++; /* clustal format  */
                }

                if(strstr(line, "CLUSTAL O")){
                        hints[2]++; /* clustal format  */
                }
                if(strstr(line, "!!AA_MULTIPLE_ALIGNMENT")){
                        hints[1]++;
                }

                if(strstr(line, "!!NA_MULTIPLE_ALIGNMENT")){
                        hints[1]++;
                }
                if(strstr(line, "MSF:")){
                        hints[1]++;
                }
                line_number++;
                if(line_number == 100){
                        break;
                }
        }

        fclose(f_ptr);
        set = 0;

        for(i = 0; i < 3;i++){
                if(hints[i]){
                        set++;
                }
        }
        if(set == 0){
                ERROR_MSG("Input alignment format could not be detected.");
        }
        if(set > 1){
                ERROR_MSG("Input format could not be unambiguously detected");
        }


        if(hints[0]){
                *type = FORMAT_FA;
        }
        if(hints[1]){
                *type = FORMAT_MSF;
        }
        if(hints[2]){
                *type = FORMAT_CLU;
        }
        //fprintf(stdout,"fa: %d msf:%d clu:%d", hints[0],hints[1],hints[2]);

        return OK;
ERROR:
        return FAIL;

}

/* detect alphabet */
int detect_alphabet(struct msa* msa)
{
        int i;
        int min,c;
        uint8_t DNA[128];
        uint8_t protein[128];
        int diff[3];
        char DNA_letters[]= "acgtuACGTUnN";
        char protein_letters[] = "acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY";

        ASSERT(msa != NULL, "No alignment");

        for(i = 0; i <128;i++){
                DNA[i] = 0;
                protein[i] = 0;
        }

        for(i = 0 ; i < strlen(DNA_letters);i++){
                DNA[(int) DNA_letters[i]] = 1;
        }

        for(i = 0 ; i < strlen(protein_letters);i++){
                protein[(int) protein_letters[i]] = 1;
        }

        diff[0] = 0;
        diff[1] = 0;
        for(i = 0; i < 128;i++){
                if((msa->letter_freq[i]) && (!DNA[i])){
                        diff[0]++;
                }
                if((msa->letter_freq[i]) && (!protein[i])){
                        diff[1]++;
                }
        }

        if( diff[0] + diff[1] == 0){
                ERROR_MSG("Could not detect any AA or nucleotides.");
        }
        c = -1;
        min = 2147483647;
        for(i = 0; i < 2;i++){
                if(diff[i] < min){
                        min = diff[i];
                        c = i;
                }
        }
        if(c == 0){
                LOG_MSG("Detected DNA sequences.");
                msa->L = defDNA;
                RUN(convert_msa_to_internal(msa, defDNA));
        }else if(c == 1){
                LOG_MSG("Detected protein sequences.");
                msa->L = redPROTEIN;
                RUN(convert_msa_to_internal(msa, redPROTEIN));
        }else{
                ERROR_MSG("Alphabet not recognized.");
        }
        return OK;
ERROR:
        return FAIL;
}

int detect_aligned(struct msa* msa)
{
        /* Here I look for gap character in the character frequency vector
           if they are there we know the input sequences were aligned. This is not strictly
           necessary for .clu and .msf infiles - but better safe than sorry. */
        int min_len;
        int max_len;
        int i;
        int gaps = 0;
        /* assume that sequences are not aligned */
        msa->aligned = 0;

        gaps += msa->letter_freq[ (int) '-'];
        gaps += msa->letter_freq[ (int) '.'];
        gaps += msa->letter_freq[ (int) '_'];
        gaps += msa->letter_freq[ (int) '~'];



        if(gaps){
                msa->aligned = 1;
        }
        /* if all sequences are the same length they could be considered
           aligned (at least for the purpose of writing out a .msf / .clu
           file) */
        min_len = INT32_MAX;
        max_len = INT32_MIN;
        for(i = 0; i < msa->numseq;i++){
                min_len = MACRO_MIN(min_len, msa->sequences[i]->len);
                max_len = MACRO_MAX(max_len, msa->sequences[i]->len);
        }
        //LOG_MSG("%d %d",min_len,max_len);
        if(min_len == max_len){
                msa->aligned = 1;
        }
        return OK;
}


int dealign_msa(struct msa* msa)
{
        struct msa_seq* seq = NULL;
        int i;
        int j;

        for(i = 0; i < msa->numseq;i++){
                seq = msa->sequences[i];
                for(j = 0; j <=  seq->len;j++){
                        seq->gaps[j] = 0;
                }
        }
        return OK;
}


int set_sip_nsip(struct msa* msa)
{
        int i;
        ASSERT(msa!= NULL, "No msa");
        if(msa->plen){
                for (i = msa->num_profiles;i--;){
                        if(msa->sip[i]){
                                MFREE(msa->sip[i]);
                        }
                }
                if(msa->plen){
                        MFREE(msa->plen);
                }
                if(msa->sip){
                        MFREE(msa->sip);
                }
                if(msa->nsip){
                        MFREE(msa->nsip);
                }
                msa->plen = NULL;
                msa->sip = NULL;
                msa->nsip = NULL;
        }

        msa->num_profiles = (msa->numseq << 1 )-1;

        MMALLOC(msa->sip,sizeof(int*)* msa->num_profiles);
        MMALLOC(msa->nsip,sizeof(int)* msa->num_profiles);
        MMALLOC(msa->plen,sizeof(int)* msa->num_profiles);


        for (i =0;i < msa->num_profiles;i++){
                msa->sip[i] = NULL;
                msa->nsip[i] = 0;

        }

        for(i = 0;i < msa->numseq;i++){

                MMALLOC(msa->sip[i],sizeof(int));
                msa->nsip[i] = 1;
                msa->sip[i][0] = i;
                msa->plen[i] = 0;
        }
        return OK;
ERROR:
        return FAIL;
}

int convert_msa_to_internal(struct msa* msa, int type)
{
        struct alphabet* a = NULL;
        struct msa_seq* seq = NULL;
        int8_t* t = NULL;
        int i,j;

        RUNP(a = create_alphabet(type));

        t = a->to_internal;
        msa->L = a->L;
        for(i = 0; i <  msa->numseq;i++){
                seq = msa->sequences[i];
                for(j =0 ; j < seq->len;j++){
                        if(t[(int) seq->seq[j]] == -1){
                                WARNING_MSG("there should be no character not matching the alphabet");
                                WARNING_MSG("offending character: >>>%c<<<", seq->seq[j]);
                        }else{
                                seq->s[j] = t[(int) seq->seq[j]];
                        }
                }

        }
        MFREE(a);
        return OK;
ERROR:
        if(a){
                MFREE(a);
        }
        return FAIL;
}

/* rw functions; I wand fasta, msf and clustal */

int write_msa_fasta(struct msa* msa,char* outfile)
{
        FILE* f_ptr = NULL;
        int i,j,c,f;

        if(!outfile){
                f_ptr = stdout;
        }else{
                RUNP(f_ptr = fopen(outfile, "w"));
        }

        for(i = 0; i < msa->numseq;i++){
                fprintf(f_ptr,">%s\n", msa->sequences[i]->name);
                f = 0;
                for(j = 0;j < msa->sequences[i]->len;j++){

                        for(c = 0;c < msa->sequences[i]->gaps[j];c++){
                                fprintf(f_ptr,"-");
                                f++;
                                if(f == 60){
                                        fprintf(f_ptr, "\n");
                                        f = 0;
                                }
                        }
                        fprintf(f_ptr,"%c", msa->sequences[i]->seq[j]);
                        f++;
                        if(f == 60){
                                        fprintf(f_ptr, "\n");
                                        f = 0;
                                }
                }
                for(c = 0;c < msa->sequences[i]->gaps[ msa->sequences[i]->len];c++){
                        fprintf(f_ptr,"-");
                        f++;
                        if(f == 60){
                                fprintf(f_ptr, "\n");
                                f = 0;
                        }
                }
                if(f){
                        fprintf(f_ptr,"\n");
                }
        }
        if(outfile){
                fclose(f_ptr);
        }

        return OK;
ERROR:
        return FAIL;
}

int write_msa_msf(struct msa* msa,char* outfile)
{
        struct line_buffer* lb = NULL;
        struct out_line* ol= NULL;
        struct msa_seq* seq = NULL;
        time_t now;			/* current time as a time_t */
        char   date[64];		/* today's date in GCG's format "October 3, 1996 15:57" */
        char* linear_seq = NULL;
        char* ptr;
        FILE* f_ptr = NULL;

        int aln_len;
        int i,j,c,f;
        int block;
        int max_name_len = 0;
        int line_length;
        int header_index;

        if(!msa->aligned){
                ERROR_MSG("Sequences appear to be unaligned!");
        }

        if(!outfile){
                f_ptr = stdout;
        }else{
                RUNP(f_ptr = fopen(outfile, "w"));
        }

        for(i = 0; i < msa->numseq;i++){
                max_name_len = MACRO_MAX(max_name_len, strnlen( msa->sequences[i]->name,MSA_NAME_LEN));
        }

        aln_len = 0;
        for (j = 0; j <= msa->sequences[0]->len;j++){
                aln_len+=  msa->sequences[0]->gaps[j];
        }
        aln_len += msa->sequences[0]->len;

        MMALLOC(linear_seq, sizeof(char)* (aln_len+1));


        /* length of write line buffer should be:
           max_name_len + 5 +
           60 (for sequences) +
           1 for newline character +
           6 for spaces every 10 letters
        */
        line_length = max_name_len +5 + 60 + 2+6 ;

        RUNP(lb = alloc_line_buffer(line_length));

        /* print header line
           here I will use seq index -1 to make sure this ends up on top.
         */
        /*
!!AA_MULTIPLE_ALIGNMENT 1.0

  stdout MSF:  131 Type: P 16/01/02 CompCheck: 3003 ..

  Name: IXI_234 Len: 131  Check: 6808 Weight: 1.00
  Name: IXI_235 Len: 131  Check: 4032 Weight: 1.00
  Name: IXI_236 Len: 131  Check: 2744 Weight: 1.00
  Name: IXI_237 Len: 131  Check: 9419 Weight: 1.00
//
        */
        header_index = -1 * (msa->numseq+10);
        ol = lb->lines[lb->num_line];
        //LOG_MSG("Alphabet : %d", msa->L);
        if(msa->L == defPROTEIN){
                snprintf(ol->line, line_length,"!!AA_MULTIPLE_ALIGNMENT 1.0");
        }else if(msa->L == redPROTEIN){
                snprintf(ol->line, line_length,"!!AA_MULTIPLE_ALIGNMENT 1.0");
        }else if(msa->L == defDNA){
                snprintf(ol->line, line_length,"!!NA_MULTIPLE_ALIGNMENT 1.0");
        }else{
                snprintf(ol->line, line_length,"!!NA_MULTIPLE_ALIGNMENT 1.0");
        }
        ol->block = -1;
        ol->seq_id = header_index;
        header_index++;
        lb->num_line++;
        /* a space */
        ol = lb->lines[lb->num_line];
        ol->line[0] = 0;
        ol->block = -1;
        ol->seq_id = header_index;
        header_index++;
        lb->num_line++;

        /* The msf line*/
        now = time(NULL);
        if (strftime(date, 64, "%B %d, %Y %H:%M", localtime(&now)) == 0){
                ERROR_MSG("time failed???");
        }
        ol = lb->lines[lb->num_line];
        snprintf(ol->line, line_length," %s  MSF: %d  Type: %c  %s  Check: %d  ..", outfile == NULL ? "kalign" : outfile,aln_len, msa->L == defPROTEIN ? 'P' : 'N', date, GCGMultchecksum(msa));

        ol->block = -1;
        ol->seq_id = header_index;
        header_index++;
        lb->num_line++;

        /* another space */
        ol = lb->lines[lb->num_line];
        ol->line[0] = 0;
        ol->block = -1;
        ol->seq_id = header_index;
        header_index++;
        lb->num_line++;

        //Name: IXI_234 Len: 131  Check: 6808 Weight: 1.00
        for(i = 0; i < msa->numseq;i++){
                if(lb->alloc_num_lines == lb->num_line){
                        resize_line_buffer(lb);
                }

                ol = lb->lines[lb->num_line];

                snprintf(ol->line, line_length," Name: %-*.*s  Len:  %5d  Check: %4d  Weight: %.2f",
                         max_name_len,max_name_len,
                         msa->sequences[i]->name ,
                         aln_len,
                         GCGchecksum(msa->sequences[i]->seq, msa->sequences[i]->len),
                         1.0);

                ol->block = -1;
                ol->seq_id = header_index;
                header_index++;
                lb->num_line++;

        }
        /* another space */
        if(lb->alloc_num_lines == lb->num_line){
                resize_line_buffer(lb);
        }
        ol = lb->lines[lb->num_line];
        ol->line[0] = 0;
        ol->block = -1;
        ol->seq_id = header_index;
        header_index++;
        lb->num_line++;

        /* header section finished */
        if(lb->alloc_num_lines == lb->num_line){
                resize_line_buffer(lb);
        }
        ol = lb->lines[lb->num_line];
        snprintf(ol->line, line_length,"//");
        ol->block = -1;
        ol->seq_id = header_index;
        header_index++;
        lb->num_line++;

        /* another space */
        if(lb->alloc_num_lines == lb->num_line){
                resize_line_buffer(lb);
        }
        ol = lb->lines[lb->num_line];
        ol->line[0] = 0;
        ol->block = -1;
        ol->seq_id = header_index;
        header_index++;
        lb->num_line++;


        /* now the actual sequence lines  */
        for(i = 0; i < msa->numseq;i++){
                block = 0;
                seq = msa->sequences[i];
                RUN(make_linear_sequence(seq,linear_seq));

                f = 0;
                while(1){
                        if(lb->alloc_num_lines == lb->num_line){
                                resize_line_buffer(lb);
                        }
                        ol = lb->lines[lb->num_line];
                        c = strnlen(seq->name, MSA_NAME_LEN);


                        for(j = 0;j < c;j++){
                                ol->line[j] = seq->name[j];
                        }
                        for(j = c;j <  max_name_len+5;j++){
                                ol->line[j] = ' ';
                        }
                        ptr = ol->line + max_name_len+5;
                        for(j = 0;j < 60;j++){
                                if(f == aln_len){

                                        ptr[j] = 0;
                                        break;
                                }
                                ptr[j] = linear_seq[f];
                                f++;
                        }
                        ptr[j] = 0;



                        ol->block = block;
                        ol->seq_id = i;

                        lb->num_line++;

                        if(i == 0){
                                if(lb->alloc_num_lines == lb->num_line){
                                        resize_line_buffer(lb);
                                }

                                ol = lb->lines[lb->num_line];
                                ol->block = block;
                                ol->seq_id = msa->numseq;
                                ol->line[0] = '\n';
                                ol->line[1] = 0;
                                lb->num_line++;
                        }

                        block++;
                        if(f == aln_len){
                                break;
                        }
                }


        }


        qsort(lb->lines , lb->num_line, sizeof(struct out_line *), sort_out_lines);
        for(i = 0; i < lb->num_line;i++){
                ol = lb->lines[i];
                //fprintf(stdout,"%d %d %s\n",ol->seq_id,ol->block,ol->line);
                fprintf(f_ptr, "%s\n", ol->line);
        }
        if(outfile){
                fclose(f_ptr);
        }
        free_line_buffer(lb);
        MFREE(linear_seq);
        return OK;
ERROR:
        return FAIL;
}

int write_msa_clustal(struct msa* msa,char* outfile)
{
        struct line_buffer* lb = NULL;
        struct out_line* ol= NULL;
        struct msa_seq* seq = NULL;

        char* linear_seq = NULL;
        char* ptr;
        FILE* f_ptr = NULL;

        int aln_len;
        int i,j,c,f;
        int block;
        int max_name_len = 0;
        int line_length;

        if(!msa->aligned){
                ERROR_MSG("Sequences appear to be unaligned!");
        }

        if(!outfile){
                f_ptr = stdout;
        }else{
                RUNP(f_ptr = fopen(outfile, "w"));
        }

        for(i = 0; i < msa->numseq;i++){
                max_name_len = MACRO_MAX(max_name_len, strnlen( msa->sequences[i]->name,MSA_NAME_LEN));
        }

        aln_len = 0;
        for (j = 0; j <= msa->sequences[0]->len;j++){
                aln_len+=  msa->sequences[0]->gaps[j];
        }
        aln_len += msa->sequences[0]->len;

        MMALLOC(linear_seq, sizeof(char)* (aln_len+1));


        /* length of write line buffer should be:
           max_name_len + 5 +
           60 (for sequences) +
           1 for newline character
        */
        line_length = max_name_len +5 + 60 + 2;

        RUNP(lb = alloc_line_buffer(line_length));

        /* print header line
           here I will use seq index -1 to make sure this ends up on top.
         */
        ol = lb->lines[lb->num_line];

        snprintf(ol->line, line_length,"Kalign (%s) multiple sequence alignment", PACKAGE_VERSION);
        ol->block = -1;
        ol->seq_id = -2;
        lb->num_line++;

        ol = lb->lines[lb->num_line];
        ol->line[0] = 0;
        ol->block = -1;
        ol->seq_id = -1;
        lb->num_line++;


        /* now the actual sequence lines  */

        for(i = 0; i < msa->numseq;i++){
                block = 0;
                seq = msa->sequences[i];
                RUN(make_linear_sequence(seq,linear_seq));

                f = 0;
                while(1){
                        if(lb->alloc_num_lines == lb->num_line){
                                resize_line_buffer(lb);
                        }
                        ol = lb->lines[lb->num_line];
                        c = strnlen(seq->name, MSA_NAME_LEN);


                        for(j = 0;j < c;j++){
                                ol->line[j] = seq->name[j];
                        }
                        for(j = c;j <  max_name_len+5;j++){
                                ol->line[j] = ' ';
                        }
                        ptr = ol->line + max_name_len+5;
                        for(j = 0;j < 60;j++){
                                if(f == aln_len){

                                        ptr[j] = 0;
                                        break;
                                }
                                ptr[j] = linear_seq[f];
                                f++;
                        }
                        ptr[j] = 0;



                        ol->block = block;
                        ol->seq_id = i;

                        lb->num_line++;

                        if(i == 0){
                                if(lb->alloc_num_lines == lb->num_line){
                                        resize_line_buffer(lb);
                                }

                                ol = lb->lines[lb->num_line];
                                ol->block = block;
                                ol->seq_id = msa->numseq;
                                ol->line[0] = '\n';
                                ol->line[1] = 0;
                                lb->num_line++;
                        }

                        block++;
                        if(f == aln_len){
                                break;
                        }
                }


        }


        qsort(lb->lines , lb->num_line, sizeof(struct out_line *), sort_out_lines);
        for(i = 0; i < lb->num_line;i++){
                ol = lb->lines[i];
                fprintf(f_ptr, "%s\n", ol->line);
                //fprintf(stdout,"%d %d %s\n",ol->seq_id,ol->block,ol->line);

        }
        if(outfile){
                fclose(f_ptr);
        }
        free_line_buffer(lb);
        MFREE(linear_seq);
        return OK;
ERROR:
        return FAIL;
}

struct msa* read_clu(char* infile, struct msa* msa)
{
        //struct msa* msa = NULL;
        struct msa_seq* seq_ptr = NULL;
        FILE* f_ptr = NULL;
        char line[BUFFER_LEN];
        int line_len;
        int i,j;
        char* p;
        int active_seq = 0;

        /* sanity checks  */
        if(!my_file_exists(infile)){
                ERROR_MSG("File: %s does not exist.",infile);
        }
        if(msa == NULL){
                msa = alloc_msa();
        }

        RUNP(f_ptr = fopen(infile, "r"));
        LOG_MSG("GAGA");
        /* scan through first line header  */
        while(fgets(line, BUFFER_LEN, f_ptr)){
                fprintf(stdout,"LINE: %s", line);
                line_len = strnlen(line, BUFFER_LEN);
                line[line_len-1] = 0;
                line_len--;
                break;
        }
        active_seq =0;
        while(fgets(line, BUFFER_LEN, f_ptr)){
                line_len = strnlen(line, BUFFER_LEN);
                line[line_len-1] = 0;
                line_len--;     /* last character is newline  */
                if(!line_len){
                        active_seq = 0;
                }else{
                        if(!isspace(line[0])){
                                if(msa->alloc_numseq == active_seq){
                                        RUN(resize_msa(msa));
                                }
                                seq_ptr = msa->sequences[active_seq];
                                //p = strstr(line,seq_ptr->name);
                                //if(p){
                                //LOG_MSG("Found bitsof seq %s", seq_ptr->name);
                                p = line;
                                j = 0;
                                for(i = 0;i < line_len;i++){
                                        if(isspace((int)p[i])){
                                                j = i;
                                                break;
                                        }
                                        seq_ptr->name[i] = p[i];
                                }
                                seq_ptr->name[j] =0;
                                for(i = j;i < line_len;i++){
                                        msa->letter_freq[(int)p[i]]++;
                                        if(isalpha((int)p[i])){
                                                if(seq_ptr->alloc_len == seq_ptr->len){
                                                        resize_msa_seq(seq_ptr);
                                                }
                                                seq_ptr->seq[seq_ptr->len] = p[i];
                                                seq_ptr->len++;
                                        }
                                        if(ispunct((int)p[i])){
                                                seq_ptr->gaps[seq_ptr->len]++;
                                        }
                                }
                                active_seq++;
                                msa->numseq = MACRO_MAX(msa->numseq, active_seq);

                        }

                }
                //fprintf(stdout,"%d \"%s\"\n",line_len,line);
        }
        RUN(null_terminate_sequences(msa));

        fclose(f_ptr);
        return msa;
ERROR:
        free_msa(msa);
        return NULL;
}



struct msa* read_msf(char* infile,struct msa* msa)
{
        //struct msa* msa = NULL;
        struct msa_seq* seq_ptr = NULL;
        FILE* f_ptr = NULL;
        char line[BUFFER_LEN];
        int line_len;
        int i,j;
        char* p;
        int active_seq = 0;

        /* sanity checks  */
        if(!my_file_exists(infile)){
                ERROR_MSG("File: %s does not exist.",infile);
        }
        if(msa == NULL){
                msa = alloc_msa();
        }

        RUNP(f_ptr = fopen(infile, "r"));


        while(fgets(line, BUFFER_LEN, f_ptr)){
                line_len = strnlen(line, BUFFER_LEN);
                line[line_len-1] = 0;
                line_len--;     /* last character is newline  */
                //fprintf(stdout,"%d \"%s\"\n",line_len,line);
                if(strstr(line, "//")){
                        //      LOG_MSG("Header done");
                        break;
                }
                /* look for name  */
                p = strstr(line,"Name:");
                if(p && strstr(line,"Len:")){
                        if(msa->alloc_numseq == msa->numseq){
                                RUN(resize_msa(msa));
                        }

                        p+= 5;  /* length of name: */
                        while( isspace((int)*p)){
                                p++;
                        }
                        seq_ptr = msa->sequences[active_seq];

                        for(i = 0;i < line_len;i++){
                                if(isspace((int)p[i])){
                                        seq_ptr->name[i] = 0;
                                        break;
                                }
                                seq_ptr->name[i] = p[i];
                        }
                        seq_ptr->seq[0] = 0;
                        msa->numseq++;
                        active_seq++;
                        //LOG_MSG("Got a name %s",p);
                }
        }
        active_seq =0;
        while(fgets(line, BUFFER_LEN, f_ptr)){
                line_len = strnlen(line, BUFFER_LEN);
                line[line_len-1] = 0;
                line_len--;     /* last character is newline  */
                if(!line_len){
                        active_seq = 0;
                }else{
                        if(!isspace(line[0])){
                                seq_ptr = msa->sequences[active_seq];
                                //p = strstr(line,seq_ptr->name);
                                //if(p){
                                //LOG_MSG("Found bitsof seq %s", seq_ptr->name);
                                p = line;
                                j = strnlen(seq_ptr->name, MSA_NAME_LEN);
                                p += j;
                                for(i = 0;i < line_len-j;i++){
                                        msa->letter_freq[(int)p[i]]++;
                                        if(isalpha((int)p[i])){
                                                if(seq_ptr->alloc_len == seq_ptr->len){
                                                        resize_msa_seq(seq_ptr);
                                                }
                                                seq_ptr->seq[seq_ptr->len] = p[i];
                                                seq_ptr->len++;
                                        }
                                        if(ispunct((int)p[i])){
                                                seq_ptr->gaps[seq_ptr->len]++;
                                        }
                                }

                                //}
                                active_seq++;
                        }

                }
                //fprintf(stdout,"%d \"%s\"\n",line_len,line);
        }
        RUN(null_terminate_sequences(msa));

        fclose(f_ptr);
        return msa;
ERROR:
        free_msa(msa);
        return NULL;
}

struct msa* read_fasta(char* infile,struct msa* msa)
{
        //struct msa* msa = NULL;
        struct msa_seq* seq_ptr = NULL;
        FILE* f_ptr = NULL;
        char line[BUFFER_LEN];
        int line_len;
        int i;

        /* sanity checks  */
        if(!my_file_exists(infile)){
                ERROR_MSG("File: %s does not exist.",infile);
        }
        if(msa == NULL){
                msa = alloc_msa();
        }

        RUNP(f_ptr = fopen(infile, "r"));


        while(fgets(line, BUFFER_LEN, f_ptr)){
                line_len = strnlen(line, BUFFER_LEN);
                if(line[0] == '>'){
                        /* alloc seq if buffer is full */
                        if(msa->alloc_numseq == msa->numseq){
                                RUN(resize_msa(msa));
                        }

                        line[line_len-1] = 0;
                        for(i =0 ; i < line_len;i++){
                                if(isspace(line[i])){
                                        line[i] = 0;

                                }

                        }
                        seq_ptr = msa->sequences[msa->numseq];
                        snprintf(seq_ptr->name ,MSA_NAME_LEN ,"%s",line+1);
                        msa->numseq++;

                }else{

                        for(i = 0;i < line_len;i++){
                                msa->letter_freq[(int)line[i]]++;
                                if(isalpha((int)line[i])){

                                        if(seq_ptr->alloc_len == seq_ptr->len){
                                                resize_msa_seq(seq_ptr);
                                        }

                                        seq_ptr->seq[seq_ptr->len] = line[i];
                                        seq_ptr->len++;
                                }
                                if(ispunct((int)line[i])){
                                        seq_ptr->gaps[seq_ptr->len]++;
                                }
                        }
                }
        }
        RUN(null_terminate_sequences(msa));
        fclose(f_ptr);
        return msa;
ERROR:
        free_msa(msa);
        return NULL;
}


int null_terminate_sequences(struct msa* msa)
{
        int i;
        struct msa_seq* seq_ptr = NULL;
        /* 0 terminate sequences  */
        for(i = 0; i < msa->numseq;i++){
                seq_ptr = msa->sequences[i];
                if(seq_ptr->alloc_len == seq_ptr->len){
                        resize_msa_seq(seq_ptr);
                }
                seq_ptr->seq[seq_ptr->len] = 0;

        }
        return OK;
}

int make_linear_sequence(struct msa_seq* seq, char* linear_seq)
{
        int c,j,f;
        f = 0;
        for(j = 0;j < seq->len;j++){
                for(c = 0;c < seq->gaps[j];c++){
                        linear_seq[f] = '-';
                        f++;

                }
                linear_seq[f] = seq->seq[j];
                f++;
        }
        for(c = 0;c < seq->gaps[ seq->len];c++){
                linear_seq[f] = '-';
                f++;
        }
        linear_seq[f] = 0;
        //fprintf(stdout,"LINEAR:%s\n",linear_seq);
        return OK;
}


/* memory functions */
struct msa* alloc_msa(void)
{
        struct msa* msa = NULL;
        int i;
        MMALLOC(msa, sizeof(struct msa));
        msa->sequences = NULL;
        msa->alloc_numseq = 512;
        msa->numseq = 0;
        msa->L = 0;
        msa->aligned = 0;
        msa->plen = NULL;
        msa->sip = NULL;
        msa->nsip = NULL;

        MMALLOC(msa->sequences, sizeof(struct msa_seq*) * msa->alloc_numseq);

        for(i = 0; i < msa->alloc_numseq;i++){
                msa->sequences[i] = NULL;
                RUNP(msa->sequences[i] = alloc_msa_seq());
        }
        for(i = 0; i < 128; i++){
                msa->letter_freq[i] = 0;
        }
        return msa;
ERROR:
        free_msa(msa);
        return NULL;
}

int resize_msa(struct msa* msa)
{
        int i;
        int old_size;

        old_size = msa->alloc_numseq;
        msa->alloc_numseq = msa->alloc_numseq + 512;

        MREALLOC(msa->sequences, sizeof(struct msa_seq*) * msa->alloc_numseq);

        for(i = old_size; i < msa->alloc_numseq;i++){
                msa->sequences[i] = NULL;
                RUNP(msa->sequences[i] = alloc_msa_seq());
        }
        return OK;
ERROR:
        return FAIL;
}

void free_msa(struct msa* msa)
{
        int i;
        if(msa){
                for(i = 0; i < msa->alloc_numseq;i++){
                        free_msa_seq(msa->sequences[i]);
                }

                for (i = msa->num_profiles;i--;){
                        if(msa->sip[i]){
                                MFREE(msa->sip[i]);
                        }
                }
                MFREE(msa->plen);
                MFREE(msa->sip);
                MFREE(msa->nsip);

                MFREE(msa->sequences);
                MFREE(msa);
        }
}


struct msa_seq* alloc_msa_seq(void)
{
        struct msa_seq* seq = NULL;
        int i;
        MMALLOC(seq, sizeof(struct msa_seq));
        seq->name = NULL;
        seq->seq = NULL;
        seq->s = NULL;
        seq->gaps = NULL;
        seq->len = 0;
        seq->alloc_len = 512;

        MMALLOC(seq->name, sizeof(char)* MSA_NAME_LEN);

        MMALLOC(seq->seq, sizeof(char) * seq->alloc_len);
        MMALLOC(seq->s, sizeof(uint8_t) * seq->alloc_len);
        MMALLOC(seq->gaps, sizeof(int) * (seq->alloc_len+1));
        for(i =0;i < seq->alloc_len+1;i++){
                seq->gaps[i] = 0;
        }
        return seq;
ERROR:
        free_msa_seq(seq);
        return NULL;
}


int resize_msa_seq(struct msa_seq* seq)
{
        int old_len;
        int i;
        old_len = seq->alloc_len;
        seq->alloc_len = seq->alloc_len + 512;

        MREALLOC(seq->seq, sizeof(char) * seq->alloc_len);
        MREALLOC(seq->s, sizeof(uint8_t) * seq->alloc_len);
        MREALLOC(seq->gaps, sizeof(int) * (seq->alloc_len+1));

        for(i = old_len;i < seq->alloc_len+1;i++){
                seq->gaps[i] = 0;
        }

        return OK;
ERROR:
        return FAIL;
}

void free_msa_seq(struct msa_seq* seq)
{
        if(seq){

                MFREE(seq->name);
                MFREE(seq->seq);
                MFREE(seq->s);
                MFREE(seq->gaps);
                MFREE(seq);
        }
}



struct line_buffer* alloc_line_buffer(int max_line_len)
{
        struct line_buffer* lb = NULL;
        int i;
        ASSERT(max_line_len > 60, "max_line_len:%d too small", max_line_len);
        MMALLOC(lb, sizeof(struct line_buffer));
        lb->alloc_num_lines = 1024;

        lb->num_line = 0;
        lb->lines = NULL;
        lb->max_line_len = max_line_len;
        MMALLOC(lb->lines, sizeof(struct out_line*) * lb->alloc_num_lines);

        for(i = 0; i < lb->alloc_num_lines;i++){
                lb->lines[i] = NULL;
                MMALLOC(lb->lines[i], sizeof(struct out_line));
                lb->lines[i]->block = 0;
                lb->lines[i]->seq_id =0;
                lb->lines[i]->line = NULL;
                MMALLOC(lb->lines[i]->line,sizeof(char) * lb->max_line_len);

        }

        return lb;
ERROR:
        return NULL;
}


int resize_line_buffer(struct line_buffer* lb)
{
        int old_len = 0;
        int i;
        old_len = lb->alloc_num_lines;
        lb->alloc_num_lines = lb->alloc_num_lines + 1024;

        MREALLOC(lb->lines, sizeof(struct out_line*) * lb->alloc_num_lines);

        for(i = old_len; i < lb->alloc_num_lines;i++){
                lb->lines[i] = NULL;
                MMALLOC(lb->lines[i], sizeof(struct out_line));
                lb->lines[i]->block = 0;
                lb->lines[i]->seq_id =0;
                lb->lines[i]->line = NULL;
                MMALLOC(lb->lines[i]->line,sizeof(char) * lb->max_line_len);

        }


        return OK;
ERROR:
        return FAIL;
}


void free_line_buffer(struct line_buffer* lb)
{
        int i;

        if(lb){
                for(i = 0; i < lb->alloc_num_lines;i++){
                        MFREE(lb->lines[i]->line);
                        MFREE(lb->lines[i]);
                }
                MFREE(lb->lines);
                MFREE(lb);
        }

}


int sort_out_lines(const void *a, const void *b)
{
        struct out_line* const *one = a;
        struct out_line* const *two = b;

        if((*one)->block > (*two)->block){
                return 1;
        }else if((*one)->block == (*two)->block){
                if((*one)->seq_id > (*two)->seq_id){
                        return 1;
                }else if((*one)->seq_id == (*two)->seq_id){
                        return 0;
                }else{
                        return -1;
                }
        }else{
                return -1;
        }
}

/* alignment checksum  */
int GCGMultchecksum(struct msa* msa)
{
        int chk = 0;
        int idx;

        for (idx = 0; idx < msa->numseq; idx++){
                chk = (chk + GCGchecksum(msa->sequences[idx]->seq,  msa->sequences[idx]->len)) % 10000;
        }
        return chk;
}

/* Taken from squid library by Sean Eddy  */
int GCGchecksum(char *seq, int len)
{
        int i;			/* position in sequence */
        int chk = 0;			/* calculated checksum  */

        for (i = 0; i < len; i++){
                chk = (chk + (i % 57 + 1) * (toupper((int) seq[i]))) % 10000;
        }
        return chk;
}

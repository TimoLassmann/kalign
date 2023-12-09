#include "tldevel.h"
#include "tlmisc.h"
#include "esl_stopwatch.h"

#include "version.h"

#include "msa_struct.h"
#include "msa_alloc.h"
#include "msa_misc.h"

/* NEED TO REMOVE  */
#include "msa_op.h"

#include "alphabet.h"

#include <string.h>
#include <ctype.h>

#define  MSA_IO_IMPORT
#include "msa_io.h"



struct in_line{
        char* line;
        int len;
};

struct in_buffer{
        struct in_line** l;
        int n_lines;
        int alloc_lines;
};

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


static int check_for_sequences(struct msa* msa);

static int read_fasta(struct in_buffer* b, struct msa** msa);
static int read_msf(struct in_buffer* b, struct msa** msa);
static int read_clu(struct in_buffer* b, struct msa** msa);

static int write_msa_fasta(struct msa* msa,char* outfile);
static int write_msa_clu(struct msa* msa,char* outfile);
static int write_msa_msf(struct msa* msa,char* outfile);

static struct line_buffer* alloc_line_buffer(int max_line_len);
static int resize_line_buffer(struct line_buffer* lb);
static void free_line_buffer(struct line_buffer* lb);

static int read_file_stdin(struct in_buffer** buffer,char* infile);
static int alloc_in_buffer(struct in_buffer** buffer, int n);
static int resize_in_buffer(struct in_buffer* b);
static void free_in_buffer(struct in_buffer* b);

static int detect_alignment_format(struct in_buffer* b,int* type);

static int null_terminate_sequences(struct msa* msa);
static int sort_out_lines(const void *a, const void *b);

static int parse_format_argument(char* format, int* type);

int kalign_read_input(char* infile, struct msa** msa, int quiet)
{
        struct in_buffer* b = NULL;
        struct msa* m = NULL;
        int type;
        int i,j;
        //ASSERT(infile != NULL,"No input file");
        /* sanity checks  */
        if(infile){
                if(!my_file_exists(infile)){
                        ERROR_MSG("File: %s does not exist.",infile);
                }
        }


        DECLARE_TIMER(timer);
        START_TIMER(timer);


        /* read everything into a in_buffer  */

        RUN(read_file_stdin(&b,infile));

        /* Check if any input was read */
        /* Logic: sum length of lines 1.. up to 5 */
        //LOG_MSG("%s lines: %d", infile, b->n_lines);
        j = 0;
        for(i = 0; i < MACRO_MIN(1, b->n_lines);i++){
                j += b->l[i]->len -1; /* exclude '0' at the end of strings  */
                //fprintf(stdout,"%d %s\n", i, b->l[i]->line);
        }
        //LOG_MSG("Len: %d",j);
        if(j == 0){
                DESTROY_TIMER(timer);
                free_in_buffer(b);
                *msa = NULL;
                return OK;
        }

        RUN(detect_alignment_format(b, &type));

        //LOG_MSG("FORMAT: %d", type);
        if(type == FORMAT_FA){
                //RUNP(msa = read_msf(infile,msa));
                //RUNP(msa = read_clu(infile,msa));
                /* LOG_MSG("reading fasta"); */
                RUN(read_fasta(b,&m));
        }else if(type == FORMAT_MSF){
                RUN(read_msf(b,&m));
        }else if(type == FORMAT_CLU){
                RUN(read_clu(b,&m));
        }else if(type == FORMAT_DETECT_FAIL){
                if(infile){
                        WARNING_MSG("Could not detect input in file: %s", infile);
                }else{
                        WARNING_MSG("Could not detect input in standard input");
                }
                /* clean up allocated structures */
                free_in_buffer(b);
                DESTROY_TIMER(timer);
                *msa = NULL;
                return OK;
        }
        m->quiet = quiet;

        RUN(detect_alphabet(m));
        RUN(detect_aligned(m));
        RUN(set_sip_nsip(m));
        free_in_buffer(b);
        STOP_TIMER(timer);
        if(!quiet){
                if(infile){
                        LOG_MSG("Read %d sequences from %s.", m->numseq, infile);
                }else{
                        LOG_MSG("Read %d sequences from standard input.",m->numseq);
                }
                GET_TIMING(timer);
        }
        DESTROY_TIMER(timer);

        if(*msa != NULL){
                RUN(merge_msa(msa, m));
                kalign_free_msa(m);
        }else{
                *msa = m;
        }
        /* LOG_MSG("%d " , (*msa)->aligned); */
        RUN(check_for_sequences(*msa));
        return OK;
ERROR:
        if(m){
                kalign_free_msa(m);
        }
        return FAIL;
}

int check_for_sequences(struct msa* msa)
{
        if(!msa){
                ERROR_MSG("No sequences were found in the input files or standard input.");
        }
        if(msa->numseq < 2){
                if(msa->numseq == 0){
                        ERROR_MSG("No sequences were found in the input files or standard input.");
                }else if (msa->numseq == 1){
                        ERROR_MSG("Only 1 sequence was found in the input files or standard input");
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int kalign_write_msa(struct msa* msa, char* outfile, char* format)
{
        int type = FORMAT_FA;

        ASSERT(msa!= NULL, "No alignment");
        if(msa->aligned != ALN_STATUS_FINAL){
                ERROR_MSG("Cannot produce msa output: Sequences are not aligned.");
        }
        RUN(parse_format_argument(format, &type));
        /* LOG_MSG("Aligned? %d", msa->aligned); */
        if(type == FORMAT_FA){
                RUN(write_msa_fasta(msa, outfile));
        }else if(type == FORMAT_MSF){
                /* if(msa->aligned == ALN_STATUS_UNALIGNED || msa->aligned == ALN_STATUS_UNKNOWN){ */
                /*         ERROR_MSG("Cannot produce msa output: Sequences are not aligned."); */
                /* } */
                RUN(write_msa_msf(msa, outfile));
        }else if(type == FORMAT_CLU){
                /* if(msa->aligned == ALN_STATUS_UNALIGNED || msa->aligned == ALN_STATUS_UNKNOWN){ */
                /*         ERROR_MSG("Cannot produce clu output: Sequences are not aligned."); */
                /* } */
                RUN(write_msa_clu(msa, outfile));
        }else{
                ERROR_MSG("Output format not recognized.");
        }
        return OK;
ERROR:
        return FAIL;
}

int parse_format_argument(char* format, int* type)
{
        int out_format = -1;
        if(format){
                if(strstr(format,"msf")){
                        out_format = FORMAT_MSF;
                }else if(strstr(format,"clu")){
                        out_format = FORMAT_CLU;
                }else if(strstr(format,"fasta")){
                        out_format  = FORMAT_FA;
                }else if(strstr(format,"fa")){
                        out_format  = FORMAT_FA;
                }else{
                        ERROR_MSG("Format %s not recognized.",format);
                }
        }else{
                out_format  = FORMAT_FA;
        }
        *type = out_format;
        return OK;
ERROR:
        return FAIL;
}
/* detect alignment format  */

int detect_alignment_format(struct in_buffer*b,int* type)
{
        //FILE* f_ptr = NULL;

        char* line = NULL;
        //size_t b_len = 0;
        //ssize_t nread;

        //char line[BUFFER_LEN];
        int hints[3];
        int line_len;
        int line_number;
        int set;
        int i;
        ASSERT(b != NULL,"No input file");
        /* sanity checks  */
        //if(!my_file_exists(infile)){
        //ERROR_MSG("File: %s does not exist.",infile);
        //}

        line_number = 0;
        for(i = 0; i < 3; i++){
                hints[i] =0;
        }
        for(i = 0; i < MACRO_MIN(b->n_lines, 100);i++){
                line = b->l[i]->line;
                line_len = b->l[i]->len;

                //}
        //RUNP(f_ptr = fopen(infile, "r"));

        /* scan through first line header  */
        //while ((nread = getline(&line, &b_len, f_ptr)) != -1){
                //while(fgets(line, BUFFER_LEN, f_ptr)){
                //line_len =  nread;//strnlen(line, BUFFER_LEN);
        //line[line_len-1] = 0;

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

        //fclose(f_ptr);
        //MFREE(line);
        set = 0;

        for(i = 0; i < 3;i++){
                if(hints[i]){
                        set++;
                }
        }
        if(set == 0){
                *type = FORMAT_DETECT_FAIL;
                //ERROR_MSG("Input alignment format could not be detected.");
        }
        if(set > 1){
                *type = FORMAT_DETECT_FAIL;
                //ERROR_MSG("Input format could not be unambiguously detected");
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

int read_file_stdin(struct in_buffer** buffer,char* infile)
{
        struct in_buffer* b = NULL;
        FILE* f_ptr = NULL;
        char* line = NULL;
        char* tmp = NULL;
        size_t b_len = 0;
        ssize_t nread;
        int i;
        //char line[BUFFER_LEN];
        int line_len;

        b = *buffer;

        if(!b){
                RUN(alloc_in_buffer(&b, 1024));
        }
        if(infile){
                if(!my_file_exists(infile)){
                        ERROR_MSG("File: %s does not exist.",infile);
                }
                RUNP(f_ptr = fopen(infile, "r"));
        }else{
                f_ptr = stdin;
        }

        b->n_lines = 0;

        while ((nread = getline(&line, &b_len, f_ptr)) != -1){
                //LOG_MSG("Read %d ", nread);
                //while(fgets(line, BUFFER_LEN, f_ptr)){
                line_len = nread;
                //tmp = b->l[b->n_lines]->line;
                tmp= NULL;
                MMALLOC(tmp, sizeof(char) * (line_len+1));
                for(i = 0; i < line_len;i++){

                        if(iscntrl((int) line[i])){
                                break;
                        }
                        tmp[i] = line[i];
                }
                tmp[i] = 0;
                //memcpy(tmp, line, line_len);
                //tmp[line_len-1] = 0;
                b->l[b->n_lines]->line = tmp;
                b->l[b->n_lines]->len = i;
                b->n_lines++;
                if(b->n_lines == b->alloc_lines){
                        RUN(resize_in_buffer(b));
                }

        }
        fclose(f_ptr);
        MFREE(line);
        *buffer = b;
        return OK;
ERROR:
        return FAIL;
}


int read_fasta( struct in_buffer* b,struct msa** m)
{
        struct msa* msa = NULL;
        struct msa_seq* seq_ptr = NULL;

        char* line = NULL;
        int line_len;
        int i;
        int nl;

        if(msa == NULL){
                RUN(alloc_msa(&msa,512));
                /* msa = alloc_msa(); */
        }

        for(nl = 0; nl < b->n_lines;nl++){
                line = b->l[nl]->line;
                line_len = b->l[nl]->len;
                if(line[0] == '>'){
                        /* alloc seq if buffer is full */
                        if(msa->alloc_numseq == msa->numseq){
                                RUN(resize_msa(msa));
                        }

                        /* line[line_len-1] = 0; */
                        /*for(i = 0; i < line_len;i++){
                                if(isspace(line[i])){
                                        line[i] = 0;
                                }
                                }*/

                        seq_ptr = msa->sequences[msa->numseq];
                        MFREE(seq_ptr->name);
                        MMALLOC(seq_ptr->name, sizeof(char) * line_len);
                        memcpy(seq_ptr->name, line+1, line_len);
                        /* seq_ptr->name_len = line_len; */
                        /* snprintf(seq_ptr->name ,MSA_NAME_LEN ,"%s",line+1); */
                        //LOG_MSG("Name %d: len %d  name : %s ",msa->numseq, line_len, seq_ptr->name);
                        //fprintf(stdout,"%s\n%s\n", line, seq_ptr->name);
                        msa->numseq++;

                }else{
                        for(i = 0;i < line_len;i++){
                                msa->letter_freq[(int)line[i]]++;
                                if(isalpha((int)line[i])){
                                        if(!seq_ptr){
                                                ERROR_MSG("Encountered a sequence before encountering it's name");
                                        }
                                        seq_ptr->seq[seq_ptr->len] = line[i];
                                        seq_ptr->len++;
                                        if(seq_ptr->alloc_len == seq_ptr->len){
                                                resize_msa_seq(seq_ptr);
                                        }
                                }else if(ispunct((int)line[i])){
                                        seq_ptr->gaps[seq_ptr->len]++;
                                }
                        }

                }
        }
        RUN(null_terminate_sequences(msa));


        *m = msa;

        return OK;
ERROR:
        kalign_free_msa(msa);
        return FAIL;
}


int read_clu(struct in_buffer* b , struct msa** m)
{
        struct msa* msa = NULL;
        struct msa_seq* seq_ptr = NULL;

        char* line = NULL;
        int i,j;
        char* p;
        int active_seq = 0;
        int line_len;
        int nl,ni;

        if(msa == NULL){
                RUN(alloc_msa(&msa,512));
                /* msa = alloc_msa(); */
        }

        ni =0;
        for(nl = 0; nl < b->n_lines;nl++){
                line = b->l[nl]->line;
                line_len = b->l[nl]->len;
                ni++;
                break;
        }
        active_seq =0;
        for(nl = ni; nl < b->n_lines;nl++){
                line = b->l[nl]->line;
                line_len = b->l[nl]->len;

                if(!line_len){
                        active_seq = 0;
                }else{
                        if(!isspace(line[0])){
                                if(msa->alloc_numseq == active_seq){
                                        RUN(resize_msa(msa));
                                }
                                seq_ptr = msa->sequences[active_seq];

                                p = line;
                                j = 0;
                                for(i = 0;i < line_len;i++){
                                        if(i == MSA_NAME_LEN-1){
                                                seq_ptr->name[i] = 0;
                                                j = i;
                                                break;
                                        }
                                        if(isspace((int)p[i])){
                                                j = i;
                                                break;
                                        }
                                        seq_ptr->name[i] = p[i];
                                }
                                seq_ptr->name[j] = 0;
                                for(i = j;i < line_len;i++){
                                        msa->letter_freq[(int)p[i]]++;
                                        if(isalpha((int)p[i])){
                                                seq_ptr->seq[seq_ptr->len] = p[i];
                                                seq_ptr->len++;
                                                if(seq_ptr->alloc_len == seq_ptr->len){
                                                        resize_msa_seq(seq_ptr);
                                                }
                                        }else if(ispunct((int)p[i])){
                                                seq_ptr->gaps[seq_ptr->len]++;
                                        }
                                }
                                active_seq++;
                                msa->numseq = MACRO_MAX(msa->numseq, active_seq);
                        }
                }
        }
        RUN(null_terminate_sequences(msa));

        *m = msa;
        return OK;
ERROR:
        kalign_free_msa(msa);
        return FAIL;
}



int read_msf(struct in_buffer* b,struct msa** m)
{
        struct msa* msa = NULL;
        struct msa_seq* seq_ptr = NULL;
        char* line = NULL;
        int line_len;
        int i,j,nl,li;
        char* p;
        int active_seq = 0;
        if(msa == NULL){
                RUN(alloc_msa(&msa,512));
                /* msa = alloc_msa(); */
        }
        li = 0;
        for(nl = 0; nl < b->n_lines;nl++){
                line = b->l[nl]->line;
                line_len = b->l[nl]->len;
                li++;
                /* line_len--;     /\* last character is newline  *\/ */
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
                        p += 5;  /* length of name: */
                        while( isspace((int)*p)){
                                p++;
                        }
                        /* LOG_MSG("Found name: %s len %d", p, strlen(p)); */
                        seq_ptr = msa->sequences[active_seq];
                        for(i = 0;i < line_len;i++){
                                if(i == MSA_NAME_LEN-1){
                                        seq_ptr->name[i] = 0;
                                        break;
                                }
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
        active_seq = 0;
        for(nl = li; nl < b->n_lines;nl++){
                line = b->l[nl]->line;
                line_len = b->l[nl]->len;
                /* line_len--;     /\* last character is newline  *\/ */
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

                                                seq_ptr->seq[seq_ptr->len] = p[i];
                                                seq_ptr->len++;
                                                if(seq_ptr->alloc_len == seq_ptr->len){
                                                        resize_msa_seq(seq_ptr);
                                                }
                                        }else if(ispunct((int)p[i])){
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

        *m = msa;
        //fclose(f_ptr);
        //MFREE(line);
        return OK;
ERROR:
        kalign_free_msa(msa);
        return FAIL;
}



int write_msa_fasta(struct msa* msa,char* outfile)
{
        FILE* f_ptr = NULL;
        int i,j,f;

        if(!outfile){
                f_ptr = stdout;
        }else{
                RUNP(f_ptr = fopen(outfile, "w"));
        }

        for(i = 0; i < msa->numseq;i++){
                fprintf(f_ptr,">%s\n", msa->sequences[i]->name);
                f = 0;
                for(j = 0;j < msa->alnlen;j++){
                        /*for(c = 0;c < msa->sequences[i]->gaps[j];c++){
                                fprintf(f_ptr,"-");
                                f++;
                                if(f == 60){
                                        fprintf(f_ptr, "\n");
                                        f = 0;
                                }
                                }*/
                        fprintf(f_ptr,"%c", msa->sequences[i]->seq[j]);
                        f++;
                        if(f == 60){
                                fprintf(f_ptr, "\n");
                                f = 0;
                        }
                }
                /* for(c = 0;c < msa->sequences[i]->gaps[ msa->sequences[i]->len];c++){ */
                /*         fprintf(f_ptr,"-"); */
                /*         f++; */
                /*         if(f == 60){ */
                /*                 fprintf(f_ptr, "\n"); */
                /*                 f = 0; */
                /*         } */
                /* } */
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

int write_msa_clu(struct msa* msa,char* outfile)
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
                max_name_len = MACRO_MAX(max_name_len, (int)strnlen( msa->sequences[i]->name,MSA_NAME_LEN));
        }


        /* aln_len = 0; */
        /* for (j = 0; j <= msa->sequences[0]->len;j++){ */
        /*         aln_len += msa->sequences[0]->gaps[j]; */
        /* } */
        /* aln_len += msa->sequences[0]->len; */


        /* MMALLOC(linear_seq, sizeof(char)* (aln_len+1)); */


        /* length of write line buffer should be:
           max_name_len + 5 +
           60 (for sequences) +
           1 for newline character
        */
        line_length = max_name_len +5 + 60 + 2;
        line_length = MACRO_MAX(256, line_length );
        RUNP(lb = alloc_line_buffer(line_length));

        /* print header line
           here I will use seq index -1 to make sure this ends up on top.
         */
        ol = lb->lines[lb->num_line];

        snprintf(ol->line, line_length,"Kalign (%s) multiple sequence alignment", KALIGN_PACKAGE_VERSION);
        ol->block = -1;
        ol->seq_id = -2;
        lb->num_line++;

        ol = lb->lines[lb->num_line];
        ol->line[0] = 0;
        ol->block = -1;
        ol->seq_id = -1;
        lb->num_line++;


        aln_len = msa->alnlen;
        /* now the actual sequence lines  */
        for(i = 0; i < msa->numseq;i++){
                block = 0;
                seq = msa->sequences[i];
                /* RUN(make_linear_sequence(seq,linear_seq)); */
                linear_seq = seq->seq;// msa->sequences[i]->seq;
                /* aln_len = seq->len;//msa->sequences[i]->len; */
                f = 0;
                while(1){
                        if(lb->alloc_num_lines == lb->num_line){
                                resize_line_buffer(lb);
                        }
                        ol = lb->lines[lb->num_line];
                        c = (int) strnlen(seq->name, MSA_NAME_LEN);


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
        /* MFREE(linear_seq); */
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
        struct tm local_time;
        char   date[64];		/* today's date in GCG's format "October 3, 1996 15:57" */
        char* linear_seq = NULL;
        char* ptr;
        char* basename = NULL;
        FILE* f_ptr = NULL;

        int aln_len;
        int i,j,c,f;
        int block;
        int max_name_len = 0;
        int line_length;
        int header_index;
        int written;

        if(!msa->aligned){
                ERROR_MSG("Sequences appear to be unaligned!");
        }

        if(!outfile){
                f_ptr = stdout;
        }else{
                RUNP(f_ptr = fopen(outfile, "w"));
        }

        for(i = 0; i < msa->numseq;i++){
                max_name_len = MACRO_MAX(max_name_len, (int)strnlen( msa->sequences[i]->name,MSA_NAME_LEN));
        }

        aln_len = msa->sequences[0]->len;
        /* for (j = 0; j <= msa->sequences[0]->len;j++){ */
        /*         aln_len+=  msa->sequences[0]->gaps[j]; */
        /* } */
        /* aln_len += msa->sequences[0]->len; */

        /* MMALLOC(linear_seq, sizeof(char)* (aln_len+1)); */


        /* length of write line buffer should be:
           max_name_len + 5 +
           60 (for sequences) +
           1 for newline character +
           6 for spaces every 10 letters
        */
        /* BUG: the header lines might be longer!!  */
        line_length = max_name_len +5 + 60 + 2+6 ;
        line_length = MACRO_MAX(256, line_length);
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
        if(msa->L == ALPHA_defPROTEIN){
                snprintf(ol->line, line_length,"!!AA_MULTIPLE_ALIGNMENT 1.0");
        }else if(msa->L == ALPHA_redPROTEIN){
                snprintf(ol->line, line_length,"!!AA_MULTIPLE_ALIGNMENT 1.0");
        }else if(msa->L == ALPHA_defDNA){
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

        if((localtime_r(&now,&local_time)) == NULL){
                ERROR_MSG("could not get local time");
        }

        if (strftime(date, 64, "%B %d, %Y %H:%M", &local_time) == 0){
                ERROR_MSG("time failed???");
        }
        ol = lb->lines[lb->num_line];
        if(outfile){
                RUN(tlfilename(outfile, &basename));
        }

        written = snprintf(ol->line, line_length," %s  MSF: %d  Type: %c  %s  Check: %d  ..", outfile == NULL ? "stdout" :   basename,aln_len, msa->L == ALPHA_defPROTEIN ? 'P' : 'N', date, GCGMultchecksum(msa));

        if(written >= line_length){
                MREALLOC(lb->lines[lb->num_line]->line,sizeof(char) * (written+1));
                ol = lb->lines[lb->num_line];
                written = snprintf(ol->line, written+1," %s  MSF: %d  Type: %c  %s  Check: %d  ..", outfile == NULL ? "stdout" : basename,aln_len, msa->L == ALPHA_defPROTEIN ? 'P' : 'N', date, GCGMultchecksum(msa));

        }

        if(basename){
                MFREE(basename);
        }

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

                written = snprintf(ol->line, line_length," Name: %-*.*s  Len:  %5d  Check: %4d  Weight: %.2f",
                                   max_name_len,max_name_len,
                                   msa->sequences[i]->name ,
                                   aln_len,
                                   GCGchecksum(msa->sequences[i]->seq, msa->sequences[i]->len),
                                   1.0);
                if(written >= line_length){
                        MREALLOC(lb->lines[lb->num_line]->line,sizeof(char) * (written+1));
                        ol = lb->lines[lb->num_line];
                        written = snprintf(ol->line, written+1," Name: %-*.*s  Len:  %5d  Check: %4d  Weight: %.2f",
                                           max_name_len,max_name_len,
                                           msa->sequences[i]->name ,
                                           aln_len,
                                           GCGchecksum(msa->sequences[i]->seq, msa->sequences[i]->len),
                                           1.0);
                }
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

        aln_len = msa->alnlen;
        /* now the actual sequence lines  */
        for(i = 0; i < msa->numseq;i++){
                block = 0;
                seq = msa->sequences[i];
                //RUN(make_linear_sequence(seq,linear_seq));
                linear_seq = seq->seq;
                /* aln_len = seq->len; */
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
        /* MFREE(linear_seq); */
        return OK;
ERROR:
        if(linear_seq){
                MFREE(linear_seq);
        }
        if(basename){
                MFREE(basename);
        }
        return FAIL;
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

int alloc_in_buffer(struct in_buffer** buffer, int n)
{
        struct in_buffer* b = NULL;
        int i;

        MMALLOC(b, sizeof(struct in_buffer));
        b->alloc_lines = n;
        b->l = NULL;
        b->n_lines = 0;

        MMALLOC(b->l, sizeof(struct in_line*) * b->alloc_lines);
        for(i = 0; i < b->alloc_lines;i++){
                b->l[i] = NULL;
                MMALLOC(b->l[i], sizeof(struct in_line));
                b->l[i]->line = NULL;
                b->l[i]->len =0;
        }

        *buffer = b;
        return OK;
ERROR:
        return FAIL;
}

int resize_in_buffer(struct in_buffer* b)
{
        int i,o;
        o = b->alloc_lines;
        b->alloc_lines = b->alloc_lines + b->alloc_lines / 2;
        MREALLOC(b->l, sizeof(struct in_line*) * b->alloc_lines);
        for(i = o; i < b->alloc_lines;i++){
                b->l[i] = NULL;
                MMALLOC(b->l[i], sizeof(struct in_line));
                b->l[i]->line = NULL;
                b->l[i]->len =0;
        }
        return OK;
ERROR:
        return FAIL;
}

void free_in_buffer(struct in_buffer* b)
{
        int i;
        if(b){
                for(i = 0; i < b->n_lines ;i++){
                        MFREE(b->l[i]->line);

                }
                for(i = 0; i < b->alloc_lines;i++){
                        MFREE(b->l[i]);
                }
                MFREE(b->l);
                MFREE(b);
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

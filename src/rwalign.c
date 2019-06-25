#include "global.h"


#define MSA_NAME_LEN 128

struct msa_seq{
        char* name;
        char* seq;
        uint8_t* s;
        int* gaps;
        int len;
        int alloc_len;
};

struct msa{
        struct msa_seq** sequences;
        int letter_freq[128];
        int numseq;
        int alloc_numseq;
        int L;

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


/* rw functions */

struct msa* read_fasta(char* infile);
struct msa* read_msf(char* infile);
struct msa* read_clu(char* infile);

int write_msa_fasta(struct msa* msa,char* outfile);
int write_msa_clustal(struct msa* msa,char* outfile);
/* memory functions  */
struct msa* alloc_msa(void);
int resize_msa(struct msa* msa);
void free_msa(struct msa* msa);

struct msa_seq* alloc_msa_seq(void);
int resize_msa_seq(struct msa_seq* seq);
void free_msa_seq(struct msa_seq* seq);



struct line_buffer* alloc_line_buffer(int max_line_len);
int resize_line_buffer(struct line_buffer* lb);
void free_line_buffer(struct line_buffer* lb);

/* local helper functions  */
static int null_terminate_sequences(struct msa* msa);
static int sort_out_lines(const void *a, const void *b);

int print_msa(struct msa* msa);
int main(int argc, char *argv[])
{
        struct msa* msa = NULL;
        LOG_MSG("Start io tests.");
        RUNP(msa = read_clu(argv[1]));
        print_msa(msa);
        write_msa_clustal(msa,NULL);

        write_msa_fasta(msa, NULL);
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


/* rw functions; I wand fasta, msf and clustal */

int write_msa_fasta(struct msa* msa,char* outfile)
{
        FILE* f_ptr = NULL;
        int i,j,c,f;

        if(!outfile){
                f_ptr = stdout;
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
        int pos;
        int max_name_len = 0;
        int line_length;
        if(!outfile){
                f_ptr = stdout;
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

        snprintf(ol->line, line_length,"\n");
        ol->block = -1;
        ol->seq_id = -1;
        lb->num_line++;


        /* now the actual sequence lines  */

        for(i = 0; i < msa->numseq;i++){
                block = 0;
                seq = msa->sequences[i];
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
                fprintf(stdout,"LINEAR:%s\n",linear_seq);

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
                fprintf(stdout,"%d %d %s\n",ol->seq_id,ol->block,ol->line);

        }

        free_line_buffer(lb);
        MFREE(linear_seq);
        return OK;
ERROR:
        return FAIL;
}

struct msa* read_clu(char* infile)
{
        struct msa* msa = NULL;
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

        msa = alloc_msa();

        RUNP(f_ptr = fopen(infile, "r"));

        /* scan through first line header  */
        while(fgets(line, BUFFER_LEN, f_ptr)){
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
                                for(i = 0;i < line_len;i++){
                                        if(isspace((int)p[i])){
                                                j = i;
                                                break;
                                        }
                                        seq_ptr->name[i] = p[i];
                                }
                                seq_ptr->name[j] =0;
                                for(i = j;i < line_len;i++){
                                        if(isalpha((int)p[i])){
                                                if(seq_ptr->alloc_len == seq_ptr->len){
                                                        resize_msa_seq(seq_ptr);
                                                }
                                                msa->letter_freq[(int)p[i]]++;
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
                fprintf(stdout,"%d \"%s\"\n",line_len,line);
        }
        RUN(null_terminate_sequences(msa));

        fclose(f_ptr);
        return msa;
ERROR:
        free_msa(msa);
        return NULL;
}



struct msa* read_msf(char* infile)
{
        struct msa* msa = NULL;
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

        msa = alloc_msa();

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
                                        if(isalpha((int)p[i])){
                                                if(seq_ptr->alloc_len == seq_ptr->len){
                                                        resize_msa_seq(seq_ptr);
                                                }
                                                msa->letter_freq[(int)p[i]]++;
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

struct msa* read_fasta(char* infile)
{
        struct msa* msa = NULL;
        struct msa_seq* seq_ptr = NULL;
        FILE* f_ptr = NULL;
        char line[BUFFER_LEN];
        int line_len;
        int i;

        /* sanity checks  */
        if(!my_file_exists(infile)){
                ERROR_MSG("File: %s does not exist.",infile);
        }

        msa = alloc_msa();

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
                                if(isalpha((int)line[i])){

                                        if(seq_ptr->alloc_len == seq_ptr->len){
                                                resize_msa_seq(seq_ptr);
                                        }
                                        msa->letter_freq[(int)line[i]]++;
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


#include "global.h"
#include "parameters.h"
#include "align_io.h"

#include "misc.h"


struct infile_buffer{
        char* input;
        unsigned short int input_type;
        unsigned short int input_numseq;

};

struct align_io_buffer{
        struct infile_buffer** in_buf;
        int feature;
        int num_inputs;
        int numseq;
};

struct align_io_buffer* alloc_align_io_buffer(int num_infiles);
int read_all_sequences(struct alignment* aln, struct align_io_buffer* b);
int count_sequences_and_detect(struct align_io_buffer* b);
int check_out_and_errors(struct align_io_buffer* b, struct parameters* param);
void free_align_io_buffer(struct align_io_buffer* b);


struct alignment* aln_alloc(int numseq);
static char* get_input_into_string(char* infile);

int count_sequences_macsim(char* string);
int count_sequences_swissprot(char* string);
int count_sequences_uniprot(char* string);
int count_sequences_stockholm(char* string);
int count_sequences_clustalw(char* string);
int count_sequences_fasta(char* string);

struct alignment* read_alignment(struct alignment* aln,char* string);
struct alignment* read_alignment_from_swissprot(struct alignment* aln,char* string);
struct alignment* read_alignment_macsim_xml(struct alignment* aln,char* string);
struct alignment* read_alignment_uniprot_xml(struct alignment* aln,char* string);
struct alignment* read_alignment_clustal(struct alignment* aln,char* string);
struct alignment* read_alignment_stockholm(struct alignment* aln,char* string);


struct feature* read_ft(struct feature* ft,char* p);

struct alignment* read_sequences(struct alignment* aln,char* string);
struct alignment* read_sequences_macsim_xml(struct alignment* aln,char* string);
struct alignment* read_sequences_from_swissprot(struct alignment* aln,char* string);
struct alignment* read_sequences_uniprot_xml(struct alignment* aln,char* string);
struct alignment* read_sequences_clustal(struct alignment* aln,char* string);
struct alignment* read_sequences_stockholm(struct alignment* aln,char* string);

struct alignment* detect_and_read_sequences(struct parameters* param)
{
        struct alignment* aln = NULL;
        struct align_io_buffer* b = NULL;


        int i = 0;



        ASSERT(param!= NULL, "No parameters");
        /* if profile do something else... */
        /* allocate input buffers  */

        RUNP(b = alloc_align_io_buffer(param->num_infiles)); /* this will allocate and extra buffer for stdin (param->num_profiles) */
        /* read from stdin */
        b->in_buf[param->num_infiles]->input = get_input_into_string(NULL);
        /* try to read in as much as possible  */
        for(i = 0; i < param->num_infiles;i++){
                if(my_file_exists(param->infile[i])){
                RUNP(b->in_buf[i]->input = get_input_into_string(param->infile[i]));
                }
        }
        /* count  */
        RUN(count_sequences_and_detect(b));
        /* detect errors  */
        RUN(check_out_and_errors(b,param));

        /* allocate alignment structure */
        RUNP(aln = aln_alloc(b->numseq));
        /* read sequences */
        RUN(read_all_sequences(aln,b));

        if(aln->numseq < 2){
                ERROR_MSG("No sequences could be read.");
        }
        if(!param->format && param->outfile){
                if (byg_start("msf",param->outfile) != -1){
                        param->format = "msf";
                }else if (byg_start("clustal",param->outfile) != -1){
                        param->format = "clustal";
                }else if (byg_start("aln",param->outfile) != -1){
                        param->format = "clustal";
                }else if (byg_start("macsim",param->outfile) != -1){
                        param->format = "macsim";
                }
                fprintf(stderr,"Output file: %s, in %s format.\n",param->outfile,param->format);
        }
        free_align_io_buffer(b);
        return aln;
ERROR:
        free_align_io_buffer(b);
        free_aln(aln);
        return NULL;
}



int count_sequences_macsim(char* string)
{
        int n = 0;
        n = byg_count("<seq-name>",string);
        if(!n){
                return -1;
        }
        return n;
}

int count_sequences_swissprot(char* string)
{
        int n = 0;
        n = byg_count("ID   ",string);
        if(!n){
                return 0;
        }
        return n;
}

int count_sequences_uniprot(char* string)
{
        int n = 0;
        n = byg_count("<entry",string);
        if(!n){
                return 0;
        }
        return n;
}

int count_sequences_stockholm(char* string)
{
        char* p1 = string;
        int i = 0;
        int j = 0;
        int n = 0;
        while((i = byg_end("\n",p1))!=-1){
                p1+=i;
                if (!(byg_start("//",p1))){
                        break;
                }
                j = byg_end("#",p1);
                if(j != 1){
                        n++;
                }
        }
        if(!n){
                return 0;
        }
        return n;
}

int count_sequences_clustalw(char* string)
{
        char* p1 = string;
        int i = 0;
        int j = 0;
        int c = 0;
        int n = 0;
        int f = 0;


        while((i = byg_end("\n",p1))!=-1){
                p1+=i;
                j = byg_end(" ",p1);
                f = byg_end("\n",p1);
                if(f > 2 && f>j && j!= 1){
                        if(c ==0){
                                i = j;
                                while(p1[i] != '\n'){
                                        //if (!isspace((int)p1[i])){
                                        //	len++;
                                        //}
                                        i++;
                                }
                        }
                        c++;
                }else{
                        if (c){
                                if(c > n){
                                        n = c;
                                }
                                c =0;
                        }
                }
        }
        if(!n){
                return 0;
        }
        return n;
}

int count_sequences_fasta(char* string)
{
        int nbytes;
        int i;
        int n = 0;
        int stop = 0;
        nbytes = strlen(string);
        for (i =0;i < nbytes;i++){
                if (string[i] == '>'&& stop == 0){
                        stop = 1;
                        n++;
                }
                if (string[i] == '\n'){
                        stop = 0;
                }
        }
        if(!n){
                return 0;
        }
        return n;
}

char* get_input_into_string(char* infile)
{
        int i = 0;
        int string_length = 2;
        char c = 0;
        FILE *file = NULL;
        char* string = NULL;

        if(infile){
                if (!(file = fopen( infile, "r" ))){
                        return NULL;
                        fprintf(stderr,"Cannot open file '%s'\n", infile);
                        exit(1);
                }
                if (fseek(file,0,SEEK_END) != 0){
                        (void)fprintf(stderr, "ERROR: fseek failed\n");
                        (void)exit(EXIT_FAILURE);
                }
                i= ftell (file);
                if (fseek(file,0,SEEK_START) != 0){
                        (void)fprintf(stderr, "ERROR: fseek failed\n");
                        (void)exit(EXIT_FAILURE);
                }
                MMALLOC(string, sizeof(char)* (i+1));

                fread(string,sizeof(char), i, file);
                string[i] = 0;
                fclose(file);
        }else{
                if (!isatty(0)){
                        MMALLOC(string, sizeof(char)* (string_length));
                        while (!feof (stdin)){
                                c = getc(stdin);
                                if (i == string_length){

                                        string_length = string_length << 1;
                                        MREALLOC(string, sizeof(char)* (string_length));
                                }
                                string[i] = c;
                                i++;
                        }
                        string[i-1] = 0;
                }else{
                        return NULL;
                }
        }
        return string;
ERROR:
        return NULL;
}

struct alignment* read_sequences_from_swissprot(struct alignment* aln,char* string)
{
        //int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};
        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        int i,j,c,n;
        char* p = 0;
        p = string;
        /*numseq = byg_count("ID   ",p);
          if(!numseq){
          fprintf(stderr,"No sequences found!\n");
          exit(1);
          }
          aln = (struct alignment *) malloc(sizeof(struct alignment));
          numprofiles = (numseq << 1) - 1;
          aln->ft = 0;
          aln->si = 0;
          aln->s = malloc(sizeof(int*) * (numseq ));
          aln->seq = malloc(sizeof(char*) * (numseq ));
          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);
          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);
          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }

          for (i = numseq;i--;){
          aln->sip[i] = malloc(sizeof(int)*1);
          aln->nsip[i] = 1;
          aln->sip[i][0] = i;
          }*/
        c = 0;
        while(aln->sl[c]){
                c++;
        }


        while ((i = byg_end("ID   ",p)) != -1){
                p+=i;
                j = byg_start(" ",p);
                aln->lsn[c] = j;
                aln->sn[c] = malloc(sizeof(char)*(j+1));
                for (i = 0;i < j;i++){
                        aln->sn[c][i] = p[i];
                }
                aln->sn[c][j] = 0;
                p+= j;
                j = byg_end("SQ   ",p);
                p+= j;
                j = byg_end("\n",p);
                p+= j;
                j = byg_start("//",p);

                aln->s[c] = malloc(sizeof(int)*(j+1));
                aln->seq[c] = malloc(sizeof(char)*(j+1));
                n = 0;
                for (i = 0;i < j;i++){
                        if(isalpha((int)p[i])){
                                aln->s[c][n] = aacode[toupper(p[i])-65];
                                aln->seq[c][n] = p[i];
                                n++;
                        }
                }
                aln->s[c][n] = 0;
                aln->seq[c][n] = 0;
                aln->sl[c] = n;
                c++;
        }
        return aln;
}


struct alignment* read_alignment_from_swissprot(struct alignment* aln,char* string)
{
        //int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};
        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        int i,j,c,n;
        char* p = 0;
        p = string;
        /*numseq = byg_count("ID   ",p);
          if(!numseq){
          fprintf(stderr,"No sequences found!\n");
          exit(1);
          }
          aln = (struct alignment *) malloc(sizeof(struct alignment));
          numprofiles = (numseq << 1) - 1;
          aln->ft = 0;
          aln->si = 0;
          aln->s = malloc(sizeof(int*) * (numseq ));
          aln->seq = malloc(sizeof(char*) * (numseq ));
          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);
          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);
          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }

          for (i = numseq;i--;){
          aln->sip[i] = malloc(sizeof(int)*1);
          aln->nsip[i] = 1;
          aln->sip[i][0] = i;
          }*/
        c = 0;
        while(aln->sl[c]){
                c++;
        }

        fprintf(stderr,"found sequence:\n");
        while ((i = byg_end("ID   ",p)) != -1){
                p+=i;
                j = byg_start(" ",p);
                aln->lsn[c] = j;
                MMALLOC(aln->sn[c],sizeof(char)*(j+1));
                for (i = 0;i < j;i++){
                        aln->sn[c][i] = p[i];
                }
                aln->sn[c][j] = 0;
                p+= j;
                j = byg_end("SQ   ",p);
                p+= j;
                j = byg_end("\n",p);
                p+= j;
                j = byg_start("//",p);
                fprintf(stderr,"found sequence:\n");
                MMALLOC(aln->s[c],sizeof(int)*(j+1));
                MMALLOC(aln->seq[c],sizeof(char)*(j+1));
                n = 0;
                for (i = 0;i < j;i++){
                        if((int)p[i] > 32){
                                if(isalpha((int)p[i])){
                                        aln->s[c][n] = aacode[toupper(p[i])-65];
                                }else{
                                        aln->s[c][n] = -1;
                                }
                                fprintf(stderr,"%c",p[i]);
                                aln->seq[c][n] = p[i];
                                n++;
                        }
                }

                fprintf(stderr,"\n\n");
                aln->s[c][n] = 0;
                aln->seq[c][n] = 0;
                aln->sl[c] = n;
                c++;
        }
        return aln;
ERROR:
        free_aln(aln);
        return NULL;
}

struct alignment* read_sequences_macsim_xml(struct alignment* aln,char* string)
{
        int c = 0;
        int n = 0;
        int i = 0;
        int j = 0;
        char *p = 0;
        int max = 0;

        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        //int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};

        /*aln = (struct alignment*) malloc(sizeof(struct alignment));
          numseq = byg_count("<seq-name>",string);
          if(!numseq){
          fprintf(stderr,"No sequences found!\n");
          exit(1);
          }

          numprofiles = (numseq << 1) - 1;
          aln->s = malloc(sizeof(int*) * (numseq ));
          aln->seq = malloc(sizeof(char*) * (numseq ));
          aln->ft =  malloc(sizeof(struct feature* ) * (numseq));
          aln->si  =  malloc(sizeof(struct sequence_information* ) * (numseq));

          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);
          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);

          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }
          for(i =0;i < numseq;i++){
          aln->ft[i] = 0;
          aln->si[i] = 0;
          aln->sip[i] = malloc(sizeof(int)*1);
          aln->nsip[i] = 1;
          aln->sip[i][0] = i;
          }*/

        p = string;

        if(byg_count("<g>",p)){
                while((i = byg_start("<g>",p))!=-1){
                        p+=i;
                        j = byg_end("<r>",p);
                        for(i = 0; i< j;i++){
                                p[i] = ' ';
                        }
                        i = byg_start("</r>",p);
                        p+=i;

                        j = byg_end("</g>",p);
                        for(i = 0; i< j;i++){
                                p[i] = ' ';
                        }

                }
        }
        p = string;

        c = 0;
        while(aln->sl[c]){
                c++;
        }



        while((i = byg_end("<sequence",p))!=-1){
                p+=i;// p1 is at start of entry;
                max = byg_end("</sequence>",p);

                i = byg_end("<seq-name>",p);
                if(i < max){
                        p +=i; //p1 is at the end of the sequence name tag
                        j = byg_start("</seq-name>",p);

                        aln->lsn[c] = j;
                        aln->sn[c] = malloc(sizeof(char)*(j+1));
                        for (i = 0;i < j;i++){
                                aln->sn[c][i] = p[i];
                        }
                        aln->sn[c][j] = 0;

                }
                i = byg_end("<ftable>",p);
                if(i < max){
                        aln->ft[c] = read_ft(aln->ft[c],p);
                }
                i = byg_end("<seq-data>",p);
                if(i < max){
                        p+= i;
                        j = byg_start("</seq-data>",p);
                        aln->s[c] = malloc(sizeof(int)*(j+1));
                        aln->seq[c] = malloc(sizeof(char)*(j+1));
                        n = 0;
                        for (i = 0;i < j;i++){
                                if(isalpha((int)p[i])){
                                        aln->s[c][n] = aacode[toupper(p[i])-65];
                                        aln->seq[c][n] = p[i];
                                        n++;
                                }
                        }
                        aln->s[c][n] = 0;
                        aln->seq[c][n] = 0;
                        aln->sl[c] = n;
                }

                c++;
        }
        return aln;
}


struct alignment* read_alignment_macsim_xml(struct alignment* aln,char* string)
{
        int c = 0;
        int n = 0;
        int i = 0;
        int j = 0;
        char *p = 0;
        int max = 0;

        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        //int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};

        /*aln = (struct alignment*) malloc(sizeof(struct alignment));
          numseq = byg_count("<seq-name>",string);
          if(!numseq){
          fprintf(stderr,"No sequences found!\n");
          exit(1);
          }

          numprofiles = (numseq << 1) - 1;
          aln->s = malloc(sizeof(int*) * (numseq ));
          aln->seq = malloc(sizeof(char*) * (numseq ));
          aln->ft =  malloc(sizeof(struct feature* ) * (numseq));
          aln->si  =  malloc(sizeof(struct sequence_information* ) * (numseq));

          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);
          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);

          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }
          for(i =0;i < numseq;i++){
          aln->ft[i] = 0;
          aln->si[i] = 0;
          aln->sip[i] = malloc(sizeof(int)*1);
          aln->nsip[i] = 1;
          aln->sip[i][0] = i;
          }*/

        p = string;

        if(byg_count("<g>",p)){
                while((i = byg_start("<g>",p))!=-1){
                        p+=i;
                        j = byg_end("<r>",p);
                        for(i = 0; i< j;i++){
                                p[i] = ' ';
                        }
                        i = byg_start("</r>",p);
                        p+=i;

                        j = byg_end("</g>",p);
                        for(i = 0; i< j;i++){
                                p[i] = ' ';
                        }

                }
        }
        p = string;

        c = 0;
        while(aln->sl[c]){
                c++;
        }



        while((i = byg_end("<sequence",p))!=-1){
                p+=i;// p1 is at start of entry;
                max = byg_end("</sequence>",p);

                i = byg_end("<seq-name>",p);
                if(i < max){
                        p +=i; //p1 is at the end of the sequence name tag
                        j = byg_start("</seq-name>",p);

                        aln->lsn[c] = j;
                        MMALLOC(aln->sn[c],sizeof(char)*(j+1));
                        for (i = 0;i < j;i++){
                                aln->sn[c][i] = p[i];
                        }
                        aln->sn[c][j] = 0;

                }
                i = byg_end("<ftable>",p);
                if(i < max){
                        aln->ft[c] = read_ft(aln->ft[c],p);
                }
                i = byg_end("<seq-data>",p);
                if(i < max){
                        p+= i;
                        j = byg_start("</seq-data>",p);

                        MMALLOC(aln->s[c],sizeof(int)*(j+1));
                        MMALLOC(aln->seq[c],sizeof(char)*(j+1));
                        n = 0;
                        for (i = 0;i < j;i++){
                                if((int)p[i]>32){
                                        if(isalpha((int)p[i])){
                                                aln->s[c][n] = aacode[toupper(p[i])-65];
                                        }else{
                                                aln->s[c][n] = -1;
                                        }
                                        aln->seq[c][n] = p[i];
                                        n++;
                                }
                        }
                        aln->s[c][n] = 0;
                        aln->seq[c][n] = 0;
                        aln->sl[c] = n;
                }

                c++;
        }
        return aln;
ERROR:
        free_aln(aln);
        return NULL;
}


struct feature* read_ft(struct feature* ft,char* p)
{

        int i,j;
        struct feature *n = NULL;
        struct feature *old_n = NULL;
        char tmp[10];
        char* p1 = 0;
        p1 = p;
        while((j = byg_end("<fitem>",p1))!= -1){
                i = byg_end("</seq-info>",p1);

                if(j >i){
                        break;
                }

                n = NULL;
                MMALLOC(n,sizeof(struct feature));
                n->next = 0;
                n->color = -1;
                n->type = NULL;
                n->note = NULL;
                n->next = NULL;

                p1+=j;// p1 is at start of entry;
                i = byg_end("<ftype>",p1);
                p1 +=i; //p1 is at the end of the sequence name tag
                j = byg_start("</ftype>",p1);

                MMALLOC(n->type,sizeof(char*)*(j+1));
                for (i = 0; i < j;i++){
                        n->type[i] = p1[i];
                }
                n->type[j] = 0;

                i = byg_end("<fstart>",p1);
                p1+= i;
                j = byg_start("</fstart>",p1);

                for (i = 0; i < j;i++){
                        tmp[i] = p1[i];
                }
                tmp[j] = 0;
                n->start = atoi(tmp);
                i = byg_end("<fstop>",p1);
                p1+= i;
                j = byg_start("</fstop>",p1);
                for (i = 0; i < j;i++){
                        tmp[i] = p1[i];
                }
                tmp[j] = 0;
                n->end = atoi(tmp);

                i = byg_end("<fnote>",p1);
                p1+= i;
                j = byg_start("</fnote>",p1);
                MMALLOC(n->note,sizeof(char*)*(j+1));
                for (i = 0; i < j;i++){
                        n->note[i] = p1[i];
                }

                n->note[j] = 0;


                if((old_n = ft)!= 0){
                        while(old_n->next!=0){
                                old_n = old_n->next;
                        }
                        old_n->next = n;
                }else{
                        ft = n;
                }
                n = NULL;
        }
        return ft;
ERROR:
        return NULL;
}

struct alignment* read_sequences_uniprot_xml(struct alignment* aln,char* string)
{
        int c = 0;
        int n = 0;
        int i = 0;
        int j = 0;
        char *p1 = 0;

        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        //int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};

        /*aln = (struct alignment *) malloc(sizeof(struct alignment));
          numseq = byg_count("<entry",string);
          if(!numseq){
          fprintf(stderr,"No sequences found!\n");
          exit(1);
          }

          numprofiles = (numseq << 1) - 1;
          aln->s = malloc(sizeof(int*) * (numseq ));
          aln->seq = malloc(sizeof(char*) * (numseq ));
          aln->si = 0;
          aln->ft = 0;
          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);
          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);
          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }
          for(i =0;i < numseq;i++){
          aln->sip[i] = malloc(sizeof(int)*1);
          aln->nsip[i] = 1;
          aln->sip[i][0] = i;
          }*/

        p1 = string;


        c = 0;
        while(aln->sl[c]){
                c++;
        }

        while((i = byg_end("<entry",p1))!=-1){

                p1+=i;// p1 is at start of entry;
                i = byg_end("<name>",p1);
                p1 +=i; //p1 is at the end of the sequence name tag
                j = byg_start("</name>",p1);
                aln->lsn[c] = j;
                aln->sn[c] = malloc(sizeof(char)*(j+1));
                for (i = 0;i < j;i++){
                        aln->sn[c][i] = p1[i];
                }
                aln->sn[c][j] = 0;

                while((i = byg_end("<sequence",p1))!= -1 ){
                        i = byg_end("<sequence",p1);
                        p1+= i;
                        i = byg_end(">",p1);
                        p1 +=i;
                }

                j = byg_start("</sequence>",p1);

                aln->s[c] = malloc(sizeof(int)*(j+1));
                aln->seq[c] = malloc(sizeof(char)*(j+1));
                n = 0;
                for (i = 0;i < j;i++){
                        if(isalpha((int)p1[i])){
                                aln->s[c][n] = aacode[toupper(p1[i])-65];
                                aln->seq[c][n] = p1[i];
                                n++;
                        }
                }
                aln->s[c][n] = 0;
                aln->seq[c][n] = 0;
                aln->sl[c] = n;
                c++;
        }
        return aln;
}



struct alignment* read_alignment_uniprot_xml(struct alignment* aln,char* string)
{
        int c = 0;
        int n = 0;
        int i = 0;
        int j = 0;
        char *p1 = 0;

        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        //int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};

        /*aln = (struct alignment *) malloc(sizeof(struct alignment));
          numseq = byg_count("<entry",string);
          if(!numseq){
          fprintf(stderr,"No sequences found!\n");
          exit(1);
          }

          numprofiles = (numseq << 1) - 1;
          aln->s = malloc(sizeof(int*) * (numseq ));
          aln->seq = malloc(sizeof(char*) * (numseq ));
          aln->si = 0;
          aln->ft = 0;
          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);
          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);
          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }
          for(i =0;i < numseq;i++){
          aln->sip[i] = malloc(sizeof(int)*1);
          aln->nsip[i] = 1;
          aln->sip[i][0] = i;
          }*/

        p1 = string;


        c = 0;
        while(aln->sl[c]){
                c++;
        }

        while((i = byg_end("<entry",p1))!=-1){
                p1+=i;// p1 is at start of entry;
                i = byg_end("<name>",p1);
                p1 +=i; //p1 is at the end of the sequence name tag
                j = byg_start("</name>",p1);
                aln->lsn[c] = j;
                MMALLOC(aln->sn[c],sizeof(char)*(j+1));
                for (i = 0;i < j;i++){
                        aln->sn[c][i] = p1[i];
                }
                aln->sn[c][j] = 0;
                i = byg_end("<sequence",p1);
                p1+= i;
                i = byg_end(">",p1);
                p1 +=i;
                j = byg_start("</sequence>",p1);
                MMALLOC( aln->s[c],sizeof(int)*(j+1));
                MMALLOC(aln->seq[c],sizeof(char)*(j+1));
                n = 0;
                for (i = 0;i < j;i++){
                        if((int)p1[i] > 32){
                                if(isalpha((int)p1[i])){
                                        aln->s[c][n] = aacode[toupper(p1[i])-65];
                                }else{
                                        aln->s[c][n] = -1;
                                }
                                aln->seq[c][n] = p1[i];
                                n++;
                        }
                }
                aln->s[c][n] = 0;
                aln->seq[c][n] = 0;
                aln->sl[c] = n;
                c++;
        }
        return aln;
ERROR:
        free_aln(aln);
        return NULL;
}

struct alignment* read_sequences_stockholm(struct alignment* aln,char* string)
{
        int c = 0;
        int n = 0;
        int i = 0;
        int j = 0;
        char *p1 = 0;

        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        //int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};

        /*aln = (struct alignment*) malloc(sizeof(struct alignment));
          p1 = string;
          while((i = byg_end("\n",p1))!=-1){
          p1+=i;
          if (!(byg_start("//",p1))){
          break;
          }
          j = byg_end("#",p1);
          if(j != 1){
          numseq++;
          }
          }

          numprofiles = (numseq << 1) - 1;
          aln->s = malloc(sizeof(int*) * (numseq ));
          aln->seq = malloc(sizeof(char*) * (numseq ));
          aln->ft = 0;
          aln->si = 0;
          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);

          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);
          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }
          for(i =0;i < numseq;i++){
          aln->sip[i] = malloc(sizeof(int)*1);
          aln->nsip[i] = 1;
          aln->sip[i][0] = i;
          }*/

        c = 0;
        while(aln->sl[c]){
                c++;
        }

        p1 = string;
        while((i = byg_end("\n",p1))!=-1){
                p1+=i;
                if (!(byg_start("//",p1))){
                        break;
                }
                j = byg_end("#",p1);
                if(j != 1){
                        j = byg_start(" ",p1);
                        aln->lsn[c] = j;
                        aln->sn[c] = malloc(sizeof(char)*(j+1));
                        for (i = 0;i < j;i++){
                                aln->sn[c][i] = p1[i];
                        }
                        aln->sn[c][j] = 0;


                        p1+=j;
                        j = byg_start("\n",p1);

                        aln->s[c] = malloc(sizeof(int)*(j+1));
                        aln->seq[c] = malloc(sizeof(char)*(j+1));
                        n = 0;
                        for (i = 0;i < j;i++){
                                if(isalpha((int)p1[i])){
                                        aln->s[c][n] = aacode[toupper(p1[i])-65];
                                        aln->seq[c][n] = p1[i];
                                        n++;
                                }
                        }
                        aln->s[c][n] = 0;
                        aln->seq[c][n] = 0;
                        aln->sl[c] = n;
                        c++;
                }
        }
        return aln;
}

struct alignment* read_alignment_stockholm(struct alignment* aln,char* string)
{
        int c = 0;
        int n = 0;
        int i = 0;
        int j = 0;
        char *p1 = 0;

        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        //int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};

        /*aln = (struct alignment*) malloc(sizeof(struct alignment));
          p1 = string;
          while((i = byg_end("\n",p1))!=-1){
          p1+=i;
          if (!(byg_start("//",p1))){
          break;
          }
          j = byg_end("#",p1);
          if(j != 1){
          numseq++;
          }
          }

          numprofiles = (numseq << 1) - 1;
          aln->s = malloc(sizeof(int*) * (numseq ));
          aln->seq = malloc(sizeof(char*) * (numseq ));
          aln->ft = 0;
          aln->si = 0;
          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);

          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);
          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }
          for(i =0;i < numseq;i++){
          aln->sip[i] = malloc(sizeof(int)*1);
          aln->nsip[i] = 1;
          aln->sip[i][0] = i;
          }*/

        c = 0;
        while(aln->sl[c]){
                c++;
        }

        p1 = string;
        while((i = byg_end("\n",p1))!=-1){
                p1+=i;
                if (!(byg_start("//",p1))){
                        break;
                }
                j = byg_end("#",p1);
                if(j != 1){
                        j = byg_start(" ",p1);
                        aln->lsn[c] = j;
                        MMALLOC(aln->sn[c],sizeof(char)*(j+1));
                        for (i = 0;i < j;i++){
                                aln->sn[c][i] = p1[i];
                        }
                        aln->sn[c][j] = 0;


                        p1+=j;
                        j = byg_start("\n",p1);

                        MMALLOC(aln->s[c],sizeof(int)*(j+1));
                        MMALLOC(aln->seq[c],sizeof(char)*(j+1));
                        n = 0;
                        for (i = 0;i < j;i++){
                                if((int)p1[i] > 32){
                                        if(isalpha((int)p1[i])){
                                                aln->s[c][n] = aacode[toupper(p1[i])-65];
                                        }else{
                                                aln->s[c][n] = -1;
                                        }
                                        aln->seq[c][n] = p1[i];
                                        n++;
                                }
                        }
                        aln->s[c][n] = 0;
                        aln->seq[c][n] = 0;
                        aln->sl[c] = n;
                        c++;
                }
        }
        return aln;
ERROR:
        free_aln(aln);
        return NULL;
}


struct alignment* read_sequences_clustal(struct alignment* aln,char* string)
{
        int c = 0;
        int n = 0;
        int len = 0;
        int i = 0;
        int j = 0;
        int start = 0;
        char *p1 = 0;
        int local_numseq = 0;

        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        //int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};


        //aln = (struct alignment*) malloc(sizeof(struct alignment));
        p1 = string;

        while((i = byg_end("\n",p1))!=-1){
                p1+=i;
                j = byg_end(" ",p1);
                n = byg_end("\n",p1);
                if(n > 2 && n>j && j!= 1){
                        if(c ==0){
                                i = j;
                                while(p1[i] != '\n'){
                                        if (!isspace((int)p1[i])){
                                                len++;
                                        }
                                        i++;
                                }
                        }
                        c++;
                }else{
                        if (c){
                                if(c > local_numseq){
                                        local_numseq = c;
                                }
                                c =0;
                        }
                }
        }

        /*numprofiles = (numseq << 1) - 1;
          aln->s = malloc(sizeof(int*) * (numseq ));
          aln->seq = malloc(sizeof(char*) * (numseq ));
          aln->ft = 0;
          aln->si = 0;
          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);
          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);

          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }

          for(i =0;i < numseq;i++){
          aln->lsn[i] = 0;
          aln->sip[i] = malloc(sizeof(int)*1);
          aln->nsip[i] = 1;
          aln->sip[i][0] = i;
          aln->sl[i] = 0;*/
        start = 0;
        while(aln->sl[start]){
                start++;
        }

        for(i =start;i < local_numseq+start;i++){
                aln->s[i] = malloc(sizeof(int)*(len+1));
                aln->seq[i] = malloc(sizeof(char)*(len+1));
        }

        p1 = string;
        c = start;
        while((i = byg_end("\n",p1))!=-1){
                p1+=i;
                j = byg_end(" ",p1);
                n = byg_end("\n",p1);
                if(n > 2 && n>j && j!= 1){
                        if(aln->lsn[c] == 0){
                                aln->lsn[c] = j;
                                aln->sn[c] = malloc(sizeof(char)*(j+1));
                                for (i = 0;i < j;i++){
                                        aln->sn[c][i] = p1[i];
                                }
                                aln->sn[c][j] = 0;
                        }
                        for (i = j;i < n;i++){
                                if(isalpha((int)p1[i])){
                                        aln->s[c][aln->sl[c]] = aacode[toupper(p1[i])-65];
                                        aln->seq[c][aln->sl[c]] = p1[i];
                                        aln->sl[c]++;
                                }
                        }
                        c++;
                }else{
                        if (c != start){
                                //c =0;
                                c = start;
                        }
                }
        }
        for (i = start; i < local_numseq+start;i++){
                aln->s[i][aln->sl[i]] = 0;
        }
        return aln;
}


struct alignment* read_alignment_clustal(struct alignment* aln,char* string)
{
        int c = 0;
        int n = 0;
        int len = 0;
        int i = 0;
        int j = 0;
        int start = 0;
        char *p1 = 0;
        int local_numseq = 0;

        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        //int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};


        //aln = (struct alignment*) malloc(sizeof(struct alignment));
        p1 = string;

        while((i = byg_end("\n",p1))!=-1){
                p1+=i;
                j = byg_end(" ",p1);
                n = byg_end("\n",p1);
                if(n > 2 && n>j && j!= 1){
                        if(c ==0){
                                i = j;
                                while(p1[i] != '\n'){
                                        if ((int)p1[i] > 32){
                                                len++;
                                        }
                                        i++;
                                }
                        }
                        c++;
                }else{
                        if (c){
                                if(c > local_numseq){
                                        local_numseq = c;
                                }
                                c =0;
                        }
                }
        }

        /*numprofiles = (numseq << 1) - 1;
          aln->s = malloc(sizeof(int*) * (numseq ));
          aln->seq = malloc(sizeof(char*) * (numseq ));
          aln->ft = 0;
          aln->si = 0;
          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);
          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);

          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }

          for(i =0;i < numseq;i++){
          aln->lsn[i] = 0;
          aln->sip[i] = malloc(sizeof(int)*1);
          aln->nsip[i] = 1;
          aln->sip[i][0] = i;
          aln->sl[i] = 0;*/
        start = 0;
        while(aln->sl[start]){
                start++;
        }

        for(i =start;i < local_numseq+start;i++){
                MMALLOC(aln->s[i],sizeof(int)*(len+1));
                MMALLOC(aln->seq[i],sizeof(char)*(len+1));
        }

        p1 = string;
        c = start;
        while((i = byg_end("\n",p1))!=-1){
                p1+=i;
                j = byg_end(" ",p1);
                n = byg_end("\n",p1);
                if(n > 2 && n>j && j!= 1){
                        if(aln->lsn[c] == 0){
                                aln->lsn[c] = j;
                                MMALLOC(aln->sn[c],sizeof(char)*(j+1));
                                for (i = 0;i < j;i++){
                                        aln->sn[c][i] = p1[i];
                                }
                                aln->sn[c][j] = 0;
                        }
                        for (i = j;i < n;i++){
                                if((int)p1[i] > 32){
                                        if(isalpha((int)p1[i])){
                                                aln->s[c][aln->sl[c]] = aacode[toupper(p1[i])-65];
                                        }else{
                                                aln->s[c][aln->sl[c]] = -1;
                                        }
                                        aln->seq[c][aln->sl[c]] = p1[i];
                                        aln->sl[c]++;
                                }
                        }
                        c++;
                }else{
                        if (c != start){
                                //c =0;
                                c = start;
                        }
                }
        }
        for (i = start; i < local_numseq+start;i++){
                aln->s[i][aln->sl[i]] = 0;
                aln->seq[i][aln->sl[i]] = 0;
        }
        return aln;
ERROR:
        free_aln(aln);
        return NULL;
}

struct alignment* read_sequences(struct alignment* aln,char* string)
{
        int c = 0;
        int n = 0;
        int i = 0;
        int j = 0;
        int stop = 0;
        int start = 0;
        int nbytes;
        int local_numseq = 0;				// O	12				//U17
        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        nbytes = strlen(string);

        //aln = (struct alignment*) malloc(sizeof(struct alignment));
        for (i =0;i < nbytes;i++){
                if (string[i] == '>'&& stop == 0){
                        stop = 1;
                        local_numseq++;
                }
                if (string[i] == '\n'){
                        stop = 0;
                }
        }
        /*
          if(!numseq){
          fprintf(stderr,"No sequences found!\n");
          exit(1);
          }
          numprofiles = (numseq << 1) - 1;
          aln->s = malloc(sizeof(int*) * (numseq));
          aln->seq = malloc(sizeof(char*) * (numseq));
          aln->ft = 0;
          aln->si = 0;
          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);
          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);

          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }*/
        start = 0;
        while(aln->sl[start]){
                start++;
        }
        j = start;

        for (i =0;i < nbytes;i++){
                if (string[i] == '>' && stop == 0){
                        stop = 1;
                        aln->sl[j] =c;
                        j++;
                        c = 0;
                }
                if (string[i] == '\n'){
                        if(stop == 1){
                                aln->lsn[j-1] = n;
                                n = 0;
                        }
                        stop = 0;
                }
                if (stop == 1 && string[i] != '\n' && string[i] != 0 && string[i] != '>' ){
                        n++;
                }
                if (stop == 0 && string[i] != '\n' && string[i] != 0 ){
                        if (isalpha((int)string[i])){
                                c++;
                        }
                }
        }
        aln->sl[j] = c;

        for (i =1+start;i < local_numseq+1+start;i++){
                if(!aln->sl[i]){
                        fprintf(stderr,"Sequence %d has a length of 0!!\n",i-1);
                        exit(1);
                }
                aln->sl[i-1] = aln->sl[i];
        }
        aln->sl[start+local_numseq] = 0;

        //for (i = numseq;i--;){
        for (i = start; i < local_numseq+start;i++){
                aln->s[i] = malloc(sizeof(int)*(aln->sl[i]+1));
                aln->seq[i] = malloc(sizeof(char)*(aln->sl[i]+1));
                aln->sn[i] = malloc(sizeof(char)*(aln->lsn[i]+1));
                //aln->sip[i] = malloc(sizeof(int)*1);
                //aln->nsip[i] = 1;
                //aln->sip[i][0] = i;
        }

        stop = 0;
        j = start;
        for (i =0;i < nbytes;i++){
                if (string[i] == '>' && stop == 0 ){
                        stop = 1;
                        j++;
                        c = 0;
                }
                if (string[i] == '\n'){
                        if(stop == 1){
                                n = 0;
                        }
                        stop = 0;
                }
                if (stop == 1 &&string[i] != '\n' && string[i] != 0 && string[i] != '>' ){
                        aln->sn[j-1][n] = string[i];
                        n++;
                }
                if (stop == 0 && string[i] != '\n' && string[i] != 0 ){
                        if(isalpha((int)string[i])){
                                aln->s[j-1][c] = aacode[toupper(string[i])-65];
                                aln->seq[j-1][c] = string[i];
                                c++;
                        }
                }
        }

        for (i = start;i< local_numseq+start;i++){
                aln->s[i][aln->sl[i]] = 0;
                aln->seq[i][aln->sl[i]] = 0;
                aln->sn[i][aln->lsn[i]] = 0;
        }
        return aln;
}


struct alignment* read_alignment(struct alignment* aln,char* string)
{
        int c = 0;
        int n = 0;
        int i = 0;
        int j = 0;
        int stop = 0;
        int start = 0;
        int nbytes;
        int local_numseq = 0;				// O	12				//U17
        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        nbytes = strlen(string);

        //aln = (struct alignment*) malloc(sizeof(struct alignment));
        for (i =0;i < nbytes;i++){
                if (string[i] == '>'&& stop == 0){
                        stop = 1;
                        local_numseq++;
                }
                if (string[i] == '\n'){
                        stop = 0;
                }
        }
        /*
          if(!numseq){
          fprintf(stderr,"No sequences found!\n");
          exit(1);
          }
          numprofiles = (numseq << 1) - 1;
          aln->s = malloc(sizeof(int*) * (numseq));
          aln->seq = malloc(sizeof(char*) * (numseq));
          aln->ft = 0;
          aln->si = 0;
          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);
          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);

          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }*/
        start = 0;
        while(aln->sl[start]){
                start++;
        }
        j = start;

        for (i =0;i < nbytes;i++){
                if (string[i] == '>' && stop == 0){
                        stop = 1;
                        aln->sl[j] =c;
                        j++;
                        c = 0;
                }
                if (string[i] == '\n'){
                        if(stop == 1){
                                aln->lsn[j-1] = n;
                                n = 0;
                        }
                        stop = 0;
                }
                if (stop == 1 && string[i] != '\n' && string[i] != 0 && string[i] != '>' ){
                        n++;
                }
                if (stop == 0 && string[i] != '\n' && string[i] != 0 ){
                        if ((int)string[i] > 32){
                                c++;
                        }
                }
        }
        aln->sl[j] = c;

        for (i =1+start;i < local_numseq+1+start;i++){
                if(!aln->sl[i]){
                        fprintf(stderr,"Sequence %d has a length of 0!!\n",i-1);
                        exit(1);
                }
                aln->sl[i-1] = aln->sl[i];
        }
        aln->sl[start+local_numseq] = 0;
        //fprintf(stderr,"set to 0 : %d\n",start+local_numseq);
        //for (i = numseq;i--;){
        for (i = start; i < local_numseq+start;i++){
                //	fprintf(stderr,"len:%d %d\n",i,aln->sl[i]);
                MMALLOC(aln->s[i],sizeof(int)*(aln->sl[i]+1));
                MMALLOC(aln->seq[i],sizeof(char)*(aln->sl[i]+1));
                MMALLOC(aln->sn[i],sizeof(char)*(aln->lsn[i]+1));
                //aln->sip[i] = malloc(sizeof(int)*1);
                //aln->nsip[i] = 1;
                //aln->sip[i][0] = i;
        }

        stop = 0;
        j = start;
        for (i =0;i < nbytes;i++){
                if (string[i] == '>' && stop == 0 ){
                        stop = 1;
                        j++;
                        c = 0;
                }
                if (string[i] == '\n'){
                        if(stop == 1){
                                n = 0;
                        }
                        stop = 0;
                }
                if (stop == 1 &&string[i] != '\n' && string[i] != 0 && string[i] != '>' ){
                        aln->sn[j-1][n] = string[i];
                        n++;
                }
                if (stop == 0 && string[i] != '\n' && string[i] != 0 ){
                        if((int) string[i] > 32 ){
                                if(isalpha((int)string[i])){
                                        aln->s[j-1][c] = aacode[toupper(string[i])-65];
                                }else{
                                        aln->s[j-1][c] = -1;
                                }
                                aln->seq[j-1][c] = string[i];
                                c++;
                        }
                }
        }

        for (i = start;i< local_numseq+start;i++){
                aln->s[i][aln->sl[i]] = 0;
                aln->seq[i][aln->sl[i]] = 0;
                aln->sn[i][aln->lsn[i]] = 0;
        }
        return aln;
ERROR:

        free_aln(aln);
        return NULL;
}


struct alignment* aln_alloc(int numseq)
{
        int i;
        struct alignment* aln = NULL;

        MMALLOC(aln, sizeof(struct alignment));
        aln->s = NULL;
        aln->seq = NULL;
        aln->ft = NULL;
        aln->si = NULL;
        aln->sl = NULL;
        aln->sip = NULL;
        aln->nsip = NULL;
        aln->sn = NULL;
        aln->lsn = NULL;
        aln->numseq = numseq;
        aln->num_profiles = (numseq << 1) - 1;



        MMALLOC(aln->s,sizeof(int*) * numseq);
        MMALLOC(aln->seq,sizeof(char*) * numseq);
        MMALLOC(aln->ft,sizeof(struct feature* ) * numseq);
        MMALLOC(aln->si,sizeof(struct sequence_information* ) * numseq);
        MMALLOC(aln->sn,sizeof(char*) * numseq);
        MMALLOC(aln->lsn,sizeof(unsigned int) * numseq);


        MMALLOC(aln->sl,sizeof(unsigned int) * aln->num_profiles);
        MMALLOC(aln->sip,sizeof(unsigned int*)* aln->num_profiles);
        MMALLOC(aln->nsip,sizeof(unsigned int)* aln->num_profiles);


        for (i =0;i < aln->num_profiles;i++){
                aln->sip[i] = NULL;
                aln->nsip[i] = 0;
                aln->sl[i] = 0;
        }

        for(i =0;i < numseq;i++){
                aln->s[i] = NULL;
                aln->seq[i] = NULL;
                aln->ft[i] = NULL;
                aln->si[i] = NULL;
                aln->sn[i] = NULL;
                aln->lsn[i] = 0;

                MMALLOC(aln->sip[i],sizeof(int));
                aln->nsip[i] = 1;
                aln->sip[i][0] = i;
        }
        return aln;
ERROR:
        free_aln(aln);
        return NULL;
}


void free_aln(struct alignment* aln)
{
        int i;
        struct feature* tmp = NULL;
        struct feature* next = NULL;
        if(aln){
                for (i = aln->numseq;i--;){

                        MFREE(aln->s[i]);
                        MFREE(aln->seq[i]);
                        MFREE(aln->sn[i]);
                }

                if(aln->ft){
                        for(i = aln->numseq;i--;){
                                if(aln->ft[i]){
                                        tmp = aln->ft[i];
                                        while(tmp){
                                                next = tmp->next;
                                                MFREE(tmp->type);
                                                MFREE(tmp->note);
                                                MFREE(tmp);
                                                tmp = next;
                                        }
                                }
                        }
                        MFREE(aln->ft);
                }
                if(aln->si){
                        MFREE(aln->si);
                }

                for (i = aln->num_profiles;i--;){
                        if(aln->sip[i]){
                                MFREE(aln->sip[i]);
                        }
                }
                MFREE(aln->seq);
                MFREE(aln->s);
                MFREE(aln->sn);
                MFREE(aln->sl);
                MFREE(aln->lsn);
                MFREE(aln->sip);
                MFREE(aln->nsip);
                MFREE(aln);


        }
}


int read_all_sequences(struct alignment* aln, struct align_io_buffer* b)
{
        struct infile_buffer* in = NULL;
        int i;
        ASSERT(aln != NULL, "No alignment");
        ASSERT(b != NULL, "No IO buffer");
        for(i = 0; i < b->num_inputs;i++){
                in = b->in_buf[i];
                if(in->input){
                        switch(in->input_type){
                        case 0:
                                aln = read_sequences(aln,in->input);
                                break;
                        case 1:
                                aln = read_sequences_macsim_xml(aln,in->input);
                                break;
                        case 2:
                                aln = read_sequences_uniprot_xml(aln,in->input);
                                break;
                        case 3:
                                aln = read_sequences_from_swissprot(aln,in->input);
                                break;
                        case 4:
                                aln = read_sequences_clustal(aln,in->input);
                                break;
                        case 5:
                                aln = read_sequences_stockholm(aln,in->input);
                                break;
                        default:
                                aln = read_sequences(aln,in->input);
                                break;
                        }
                }
        }
        return OK;
ERROR:
        return FAIL;

}


int check_out_and_errors(struct align_io_buffer* b, struct parameters* param)

{
        int i;
        int c;

        int empty_file;
        c = 0;
        empty_file = -1;
        for(i = 0; i < param->num_infiles;i++){
                if(b->in_buf[i]->input_numseq == 0){
                        c++;
                        empty_file = i;
                }

        }
        if(c > 1){
                ERROR_MSG("Multiple input files have no sequences!");
        }
        if(c == 1 && !param->outfile){
                param->outfile = param->infile[empty_file];
                if(!param->format){
                        if (byg_start("msf",param->outfile) != -1){
                                param->format = "msf";
                        }else if (byg_start("clustal",param->outfile) != -1){
                                param->format = "clustal";
                        }else if (byg_start("aln",param->outfile) != -1){
                                param->format = "clustal";
                        }else if (byg_start("macsim",param->outfile) != -1){
                                param->format = "macsim";
                        }else{
                                param->format = "fasta";
                        }
                        if(param->reformat){
                                LOG_MSG("unaligned fasta format.");
                        }else if(param->format){
                                LOG_MSG("%s format.",param->format);
                        }else{
                                LOG_MSG("fasta format.");
                        }
                }
        }

        if(b->numseq < 2){
                if(b->numseq){
                        ERROR_MSG("No sequences found.");
                }else{
                        ERROR_MSG("Only one sequence found.");
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int count_sequences_and_detect(struct align_io_buffer* b)
{
        int i;
        struct infile_buffer* in = NULL;
        ASSERT(b != NULL, "No buffer.");

        for(i = 0; i < b->num_inputs;i++){
                if(b->in_buf[i]->input){
                        in = b->in_buf[i];
                        //fprintf(stdout,"%s",in->input);
                        if (byg_start("<macsim>",in->input) != -1){
                                in->input_numseq = count_sequences_macsim(in->input);
                                b->feature =1;
                                in->input_type = 1;
                        }else if (byg_start("<uniprot",in->input) != -1){
                                in->input_numseq = count_sequences_uniprot(in->input);
                                in->input_type = 2;
                        }else if(byg_start("This SWISS-PROT",in->input) != -1){
                                in->input_numseq = count_sequences_swissprot(in->input);
                                in->input_type = 3;
                        }else if (byg_start("This Swiss-Prot",in->input) != -1){
                                in->input_numseq = count_sequences_swissprot(in->input);
                                in->input_type = 3;
                        }else if (byg_start("CLUSTAL W",in->input) != -1){
                                in->input_numseq = count_sequences_clustalw(in->input);
                                in->input_type = 4;
                        }else if (byg_start("PileUp",in->input) != -1){
                                in->input_numseq = count_sequences_clustalw(in->input);
                                in->input_type = 4;
                        }else if (byg_start("MSF:",in->input) != -1){
                                in->input_numseq = count_sequences_clustalw(in->input);
                                in->input_type = 4;
                        }else if (byg_start("STOCKHOLM",in->input) != -1){
                                in->input_numseq = count_sequences_stockholm(in->input);
                                in->input_type = 5;
                        }else{
                                //fprintf(stdout,"Fasta\n");
                                in->input_numseq  = count_sequences_fasta(in->input);
                                in->input_type = 0;
                                //fprintf(stdout,"Detected %d\n", in->input_numseq);
                        }

                        if(in->input_numseq < 1){
                                MFREE(in->input);
                                in->input = NULL;
                        }else{
                                b->numseq += in->input_numseq;
                        }
                        fprintf(stdout,"%d\n",b->numseq);
                }
        }
        return OK;
ERROR:
        return FAIL;
}

struct align_io_buffer* alloc_align_io_buffer(int num_infiles)
{
        struct align_io_buffer* b = NULL;
        int i;
        MMALLOC(b, sizeof(struct align_io_buffer));
        b->num_inputs = num_infiles+1;
        b->in_buf = NULL;
        b->feature = 0;
        b->numseq = 0;
        MMALLOC(b->in_buf, sizeof(struct infile_buffer*)* b->num_inputs);
        for(i = 0; i < b->num_inputs;i++){
                b->in_buf[i] = NULL;
                MMALLOC(b->in_buf[i], sizeof(struct infile_buffer));
                b->in_buf[i]->input = NULL;
                b->in_buf[i]->input_numseq = 0;
                b->in_buf[i]->input_type = 0;
        }
        return b;
ERROR:
        free_align_io_buffer(b);
        return NULL;
}

void free_align_io_buffer(struct align_io_buffer* b)
{
        int i;
        if(b){
                if(b->in_buf){
                        for(i = 0;i< b->num_inputs;i++){
                                if(b->in_buf[i]){
                                        if(b->in_buf[i]->input){
                                                MFREE(b->in_buf[i]->input);
                                        }
                                        MFREE(b->in_buf[i]);
                                }
                        }
                        MFREE(b->in_buf);
                }
                MFREE(b);
        }
}

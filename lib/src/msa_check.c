#include "tldevel.h"

#include <string.h>
#include <ctype.h>

#include "msa_struct.h"

#define  MSA_CHECK_IMPORT
#include "msa_check.h"

struct sort_struct_name_chksum{
        struct msa_seq* seq;
        char** name;
        int chksum;
        int action;
};

static int GCGchecksum(char *seq, int len);
static int sort_by_name(const void *a, const void *b);
static int sort_by_chksum(const void *a, const void *b);
static int sort_by_both(const void *a, const void *b);

int sort_seq_by_len(const void *a, const void *b);

int kalign_sort_msa(struct msa *msa)
{
        struct sort_struct_name_chksum** a = NULL;

        MMALLOC(a, sizeof(struct sort_struct_name_chksum *) * msa->numseq);

        for(int i = 0; i < msa->numseq;i++){
                a[i] = NULL;
                MMALLOC(a[i], sizeof(struct sort_struct_name_chksum));
                a[i]->seq = msa->sequences[i];
                a[i]->name = &msa->sequences[i]->name;
                a[i]->chksum = GCGchecksum(msa->sequences[i]->seq, msa->sequences[i]->len);
                a[i]->action = 0;
        }

        qsort(a, msa->numseq, sizeof(struct sort_struct*),sort_by_both);

        for(int i = 0; i < msa->numseq;i++){
                msa->sequences[i] = a[i]->seq;
        }

        for(int i = 0; i < msa->numseq;i++){
                MFREE(a[i]);
        }
        MFREE(a);
        return OK;
ERROR:
        if(a){
                for(int i = 0; i < msa->numseq;i++){
                        MFREE(a[i]);
                }
                MFREE(a);
        }
        return FAIL;
}


int kalign_essential_input_check(struct msa *msa, int exit_on_error)
{
        int problem_len0 = 0;
        ASSERT(msa != NULL, "No alignment");

        ASSERT(msa->numseq > 1,"only %d sequences found.", msa->numseq);
        for(int i = 0; i < msa->numseq;i++){
                if(msa->sequences[i]->len == 0){
                        if(!msa->quiet){
                                WARNING_MSG("No sequence found for sequence %s ",msa->sequences[i]->name);
                        }
                        problem_len0++;
                }
                msa->sequences[i]->rank = i;
        }

        if(!exit_on_error){
                /* Here we attempt to fix the zero length problem  */
                if(problem_len0){


                        if(problem_len0 == 1){
                                if(!msa->quiet){
                                        LOG_MSG("Removing %d sequence with a length of 0.", problem_len0);
                                }
                        }else{
                                if(!msa->quiet){
                                        LOG_MSG("Removing %d sequences with a length of 0.",problem_len0);
                                }
                        }

                        struct msa_seq** tmp = NULL;
                        MMALLOC(tmp, sizeof(struct msa_seq* )  * msa->alloc_numseq);
                        int c = 0;
                        int e = msa->numseq-1;

                        for(int i = 0 ; i < msa->numseq;i++){
                                if(msa->sequences[i]->len){
                                        tmp[c] = msa->sequences[i];
                                        c++;
                                }else{
                                        tmp[e] = msa->sequences[i];
                                        e--;
                                }
                        }
                        for(int i = msa->numseq; i < msa->alloc_numseq;i++){
                                 tmp[i] = NULL;
                        }

                        MFREE(msa->sequences);
                        msa->sequences = tmp;
                        /* for(int i = msa->numseq-500; i < msa->numseq;i++){ */
                                  /* LOG_MSG("%d\t%s", msa->sequences[i]->len,msa->sequences[i]->name); */
                        /* } */
                        /* LOG_MSG("%d %d %d ", msa->numseq, msa->numseq -c , problem_len0); */
                        /* qsort(msa->sequences, msa->numseq, sizeof(struct msa_seq*),sort_seq_by_len); */
                        /* int c = 0; */
                        /* for(int i = msa->numseq-1;i >= 0;i--){ */
                        /*         if(msa->sequences[i]->len != 0){ */
                        /*                 c = i; */
                        /*                 break; */
                        /*         } */
                        /* } */
                        /* c++; */
                        msa->numseq = c;
                        ASSERT(msa->numseq > 1,"only %d sequences found.", msa->numseq);
                        /* exit(0); */
                }
        }else{
                ERROR_MSG("%d sequences found with length 0.", problem_len0);
        }

        return OK;
ERROR:
        return FAIL;
}

int kalign_check_msa(struct msa* msa, int exit_on_error)
{
        char* tmp_name = NULL;
        /* char* tmp_ptr; */
        struct sort_struct_name_chksum** a = NULL;
        int i;
        /* int j; */
        int c;
        int l;

        MMALLOC(a, sizeof(struct sort_struct_name_chksum *) * msa->numseq);

        for(i = 0; i < msa->numseq;i++){
                a[i] = NULL;
                MMALLOC(a[i], sizeof(struct sort_struct_name_chksum));
                a[i]->seq = msa->sequences[i];
                a[i]->name = &msa->sequences[i]->name;
                a[i]->chksum = GCGchecksum(msa->sequences[i]->seq, msa->sequences[i]->len);
                a[i]->action = 0;
        }

        qsort(a, msa->numseq, sizeof(struct sort_struct*),sort_by_name);

        for(i = 1; i < msa->numseq;i++){
                if(strncmp(*a[i-1]->name, *a[i]->name,MSA_NAME_LEN) == 0){
                        /* WARNING_MSG("Name: %s is duplicated", a[i]->name); */
                        if(a[i-1]->chksum == a[i]->chksum){
                                if(!msa->quiet){
                                LOG_MSG("Found duplicated sequence:\n%s checksum: %d\n%s checksum: %d\n",
                                        *a[i-1]->name,
                                        a[i-1]->chksum,
                                        *a[i]->name,
                                        a[i]->chksum);
                                }
                                a[i-1]->action = 1;
                                a[i]->action = 1;
                                if(exit_on_error){
                                        ERROR_MSG("Same seq with same name!");
                                }
                        }else{
                                if(!msa->quiet){
                                WARNING_MSG("Found sequence pair with same name but different sequence:\n%s checksum: %d\n%s checksum: %d\n",
                                            *a[i-1]->name,
                                            a[i-1]->chksum,
                                            *a[i]->name,
                                            a[i]->chksum);
                                }

                                a[i-1]->action = 1;
                                a[i]->action = 1;
                                if(exit_on_error){
                                        ERROR_MSG("This is a problem if we want to compare alignments down the track. Will append \"_X\" to the sequence name.");
                                }else{
                                        WARNING_MSG("This is a problem if we want to compare alignments down the track. Will append \"_X\" to the sequence name.");
                                }
                        }
                }
        }

        c = 1;
        for(i = 0; i < msa->numseq;i++){
                if(a[i]->action){

                        tmp_name = NULL;
                        l = strnlen(*a[i]->name, MSA_NAME_LEN)+ 16;
                        MMALLOC(tmp_name, sizeof(char) * l);

                        snprintf(tmp_name, l, "%s_%d",*a[i]->name,c);

                        MFREE(*a[i]->name);
                        *a[i]->name = tmp_name;
                        c++;
                }

        }


        qsort(a, msa->numseq, sizeof(struct sort_struct*),sort_by_chksum);
        for(i = 1; i < msa->numseq;i++){
                if(a[i-1]->chksum == a[i]->chksum){
                        if(!msa->quiet){
                        WARNING_MSG("Found identical sequences:\n%s checksum: %d\n%s checksum: %d\n",
                                    *a[i-1]->name,
                                    a[i-1]->chksum,
                                    *a[i]->name,
                                    a[i]->chksum);
                        }
                }
        }
        for(i = 0; i < msa->numseq;i++){
                MFREE(a[i]);
        }
        MFREE(a);
        return OK;
ERROR:
        if(a){
                for(i = 0; i < msa->numseq;i++){
                        MFREE(a[i]);
                }
                MFREE(a);
        }
        return FAIL;
}

int sort_by_both(const void *a, const void *b)
{
        struct sort_struct_name_chksum* const *one = a;
        struct sort_struct_name_chksum* const *two = b;

        if(strncmp(*(*one)->name, *(*two)->name,MSA_NAME_LEN) < 0){
                return -1;
        }else if(strncmp(*(*one)->name, *(*two)->name,MSA_NAME_LEN) == 0 ){
                if((*one)->chksum > (*two)->chksum){
                        return -1;
                }else{
                        return 1;
                }
        }else{
                return 1;
        }
}
int sort_by_name(const void *a, const void *b)
{
        struct sort_struct_name_chksum* const *one = a;
        struct sort_struct_name_chksum* const *two = b;

        if(strncmp(*(*one)->name, *(*two)->name,MSA_NAME_LEN) < 0){
                return -1;
        }else{
                return 1;
        }
}


int sort_by_chksum(const void *a, const void *b)
{
        struct sort_struct_name_chksum* const *one = a;
        struct sort_struct_name_chksum* const *two = b;

        if((*one)->chksum > (*two)->chksum){
                return -1;
        }else{
                return 1;
        }
}

int sort_seq_by_len(const void *a, const void *b)
{
        struct msa_seq* const *one = a;
        struct msa_seq* const *two = b;
        if((*one)->len > (*two)->len){
                return -1;
        }else{
                return 1;
        }
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

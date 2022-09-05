#include "tldevel.h"

#include <string.h>
#include <ctype.h>

#include "msa_struct.h"

#define  MSA_CHECK_IMPORT
#include "msa_check.h"

struct sort_struct_name_chksum{
        char** name;
        int chksum;
        int action;
};

static int GCGchecksum(char *seq, int len);
static int sort_by_name(const void *a, const void *b);
static int sort_by_chksum(const void *a, const void *b);
int check_msa(struct msa* msa)
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
                a[i]->name = &msa->sequences[i]->name;
                a[i]->chksum = GCGchecksum(msa->sequences[i]->seq, msa->sequences[i]->len);
                a[i]->action = 0;
        }

        qsort(a, msa->numseq, sizeof(struct sort_struct*),sort_by_name);

        for(i = 1; i < msa->numseq;i++){
                if(strncmp(*a[i-1]->name, *a[i]->name,MSA_NAME_LEN) == 0){
                        /* WARNING_MSG("Name: %s is duplicated", a[i]->name); */
                        if(a[i-1]->chksum == a[i]->chksum){
                                LOG_MSG("Found duplicated sequence:\n%s checksum: %d\n%s checksum: %d\n",
                                        *a[i-1]->name,
                                        a[i-1]->chksum,
                                        *a[i]->name,
                                        a[i]->chksum);
                                a[i-1]->action = 1;
                                a[i]->action = 1;
                        }else{
                                WARNING_MSG("Found sequence pair with same name but different sequence:\n%s checksum: %d\n%s checksum: %d\n",
                                            *a[i-1]->name,
                                            a[i-1]->chksum,
                                            *a[i]->name,
                                            a[i]->chksum);

                                a[i-1]->action = 1;
                                a[i]->action = 1;
                                WARNING_MSG("This is a problem if we want to compare alignments down the track. Will append \"_X\" to the sequence name.");
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
                        WARNING_MSG("Found identical sequences:\n%s checksum: %d\n%s checksum: %d\n",
                                    *a[i-1]->name,
                                    a[i-1]->chksum,
                                    *a[i]->name,
                                    a[i]->chksum);
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

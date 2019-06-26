#ifndef MSA_H
#define MSA_H


#define MSA_NAME_LEN 128
#define FORMAT_FA 1
#define FORMAT_MSF 2
#define FORMAT_CLU 3

#include <stdint.h>

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
        int** sip;
        int* nsip;
        int* plen;
        int numseq;
        int num_profiles;
        int alloc_numseq;
        int aligned;
        int letter_freq[128];
        int L;



};

/* dealign */
int dealign_msa(struct msa* msa);


/* convert */
int convert_msa_to_internal(struct msa* msa, int type);
/* rw functions */

struct msa* read_input(char* infile,struct msa* msa);
int write_msa(struct msa* msa, char* outfile, int type);
void free_msa(struct msa* msa);

#endif

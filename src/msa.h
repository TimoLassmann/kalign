#ifndef MSA_H
#define MSA_H


#define FORMAT_FA 1
#define FORMAT_MSF 2
#define FORMAT_CLU 3


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

        int numseq;
        int num_profiles;
        int alloc_numseq;
        int aligned;
        int letter_freq[128];
        int L;



};

/* dealign */
int dealign_msa(struct msa* msa);
/* rw functions */

struct msa* read_input(char* infile);
int write_msa(struct msa* msa, char* outfile, int type);
void free_msa(struct msa* msa);

#endif

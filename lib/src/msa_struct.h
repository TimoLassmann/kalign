#ifndef MSA_STRUCT_H
#define MSA_STRUCT_H

#include <stdint.h>

#ifdef MSA_STRUCT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define MSA_NAME_LEN 256

#define ALN_STATUS_UNALIGNED 1   /* no gaps sequences may or may not have equal lengths  */
#define ALN_STATUS_ALIGNED 2   /* sequences have equal lengths and may or may not contain gaps*/
#define ALN_STATUS_FINAL 3   /* sequences have equal lengths and may or may not contain gaps*/
#define ALN_STATUS_UNKNOWN 3     /* sequences have un-equal length and contain gaps  */

#define ALN_BIOTYPE_PROTEIN 0
#define ALN_BIOTYPE_DNA 1
#define ALN_BIOTYPE_UNDEF 2

struct msa_seq{
        char* name;
        char* seq;
        uint8_t* s;
        int* gaps;
        int rank;
        int len;
        int alloc_len;
};

struct msa{
        struct msa_seq** sequences;
        int** sip;
        int* nsip;
        int* plen;
        uint8_t run_parallel;
        int numseq;
        int num_profiles;
        int alloc_numseq;
        int aligned;
        int alnlen;
        int letter_freq[128];
        uint8_t L;
        uint8_t biotype;
        int quiet;
};

#undef MSA_STRUCT_IMPORT
#undef EXTERN

#endif

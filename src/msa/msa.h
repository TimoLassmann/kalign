#ifndef MSA_H
#define MSA_H

#include <stdint.h>

#ifdef MSA_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define MSA_NAME_LEN 128

#define ALN_STATUS_UNALIGNED 1   /* no gaps sequences may or may not have equal lengths  */
#define ALN_STATUS_ALIGNED 2   /* sequences have equal lengths and may or may not contain gaps*/
#define ALN_STATUS_UNKNOWN 3     /* sequences have un-equal length and contain gaps  */

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
        int quiet;
};

EXTERN int alloc_msa(struct msa** msa);
EXTERN int resize_msa(struct msa* msa);
EXTERN void free_msa(struct msa* msa);

EXTERN int alloc_msa_seq(struct msa_seq** s);
EXTERN int resize_msa_seq(struct msa_seq* seq);
EXTERN void free_msa_seq(struct msa_seq* seq);

EXTERN int convert_msa_to_internal(struct msa* msa, int type);

EXTERN int dealign_msa(struct msa *msa);
EXTERN int detect_alphabet(struct msa *msa);
EXTERN int detect_aligned(struct msa *msa);
EXTERN int set_sip_nsip(struct msa *msa);

#undef MSA_IMPORT
#undef EXTERN


#endif

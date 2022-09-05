#include <ctype.h>

#include "msa_struct.h"

#define MSA_MISC_IMPORT
#include "msa_misc.h"




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

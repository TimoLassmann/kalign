#include "tldevel.h"

#include "alphabet.h"

#include <ctype.h>
#define MSA_IMPORT
#include "msa.h"

static int aln_unknown_warning_message_gaps_but_len_diff(struct msa *msa);
static int aln_unknown_warning_message_same_len_no_gaps(void);

int alloc_msa(struct msa** msa)
{
        struct msa* m = NULL;
        int i;
        MMALLOC(m, sizeof(struct msa));
        m->sequences = NULL;
        m->alloc_numseq = 512;
        m->numseq = 0;
        m->num_profiles = 0;
        m->L = ALPHA_UNDEFINED;
        m->aligned = 0;
        m->quiet = 0;
        m->plen = NULL;
        m->sip = NULL;
        m->nsip = NULL;


        MMALLOC(m->sequences, sizeof(struct msa_seq*) * m->alloc_numseq);

        for(i = 0; i < m->alloc_numseq;i++){
                m->sequences[i] = NULL;
                RUN(alloc_msa_seq(&m->sequences[i]));
        }
        for(i = 0; i < 128; i++){
                m->letter_freq[i] = 0;
        }
        *msa = m;
        return OK;
ERROR:
        free_msa(m);
        return FAIL;
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
                RUN(alloc_msa_seq(&msa->sequences[i]));
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

                for (i = msa->num_profiles;i--;){
                        if(msa->sip[i]){
                                MFREE(msa->sip[i]);
                        }
                }
                if(msa->plen){
                        MFREE(msa->plen);
                }
                if(msa->sip){
                        MFREE(msa->sip);
                }
                if(msa->nsip){
                        MFREE(msa->nsip);
                }

                MFREE(msa->sequences);
                MFREE(msa);
        }
}


int alloc_msa_seq(struct msa_seq** s)
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
        *s = seq;
        return OK;

ERROR:
        free_msa_seq(seq);
        return FAIL;
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

        for(i = old_len+1;i < seq->alloc_len+1;i++){
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

int dealign_msa(struct msa* msa)
{
        struct msa_seq* seq = NULL;
        int i;
        int j;

        for(i = 0; i < msa->numseq;i++){
                seq = msa->sequences[i];
                for(j = 0; j <=  seq->len;j++){
                        seq->gaps[j] = 0;
                }
        }
        msa->aligned = ALN_STATUS_UNALIGNED;
        return OK;
}

int detect_alphabet(struct msa* msa)
{
        int i;

        double DNA[128];
        double protein[128];
        char DNA_letters[12]= "acgtunACGTUN";
        char protein_letters[40] = "acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY";

        double dna_prob;
        double prot_prob;

        ASSERT(msa != NULL, "No alignment");

        for(i = 0; i < 128;i++){
                DNA[i] = log(0.0001 * 1.0 / 116.0);
                protein[i] = log(0.0001 * 1.0 / 88.0);
        }

        for(i = 0 ; i < 12;i++){
                DNA[(int) DNA_letters[i]] = log(0.9999 * 1.0 / 12.0);
        }

        for(i = 0 ; i < 40;i++){
                protein[(int) protein_letters[i]] = log(0.9999 * 1.0 / 40.0);
        }
        /* dna_prob = 0.0; */
        /* prot_prob = 0.0; */
        /* for(i = 0; i <128;i++){ */
        /*         dna_prob += exp(DNA[i]); */

        /*         prot_prob += exp(protein[i]); */
        /* } */
        /* LOG_MSG("DNA: %f PROT: %f",dna_prob,prot_prob); */

        dna_prob = 0.0;
        prot_prob = 0.0;
        for(i = 0; i < 128;i++){
                if(msa->letter_freq[i]){
                        dna_prob += DNA[i] * (double) msa->letter_freq[i];
                        prot_prob += protein[i]* (double) msa->letter_freq[i];
                }
        }

        /* LOG_MSG("DNA: %f PROT: %f",dna_prob,prot_prob); */
        /* exit(0); */
        if( dna_prob == prot_prob){
                WARNING_MSG("Could not determine whether we have a DNA or Protein alignment");
                msa->L = ALPHA_UNKNOWN;
        }else{
                if(dna_prob > prot_prob){
                        if(!msa->quiet){
                                LOG_MSG("Detected DNA sequences.");
                        }
                        msa->L = ALPHA_defDNA;
                        RUN(convert_msa_to_internal(msa, ALPHA_defDNA));
                }else if(prot_prob > dna_prob){
                        if(!msa->quiet){
                                LOG_MSG("Detected protein sequences.");
                        }
                        msa->L = ALPHA_redPROTEIN;
                        RUN(convert_msa_to_internal(msa, ALPHA_redPROTEIN));
                }else{
                        ERROR_MSG("Alphabet not recognized.");
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int detect_aligned(struct msa* msa)
{
        int min_len;
        int max_len;
        int l;
        int i;
        int j;
        int n;
        int gaps = 0;
        /* assume that sequences are not aligned */
        msa->aligned = 0;

        /* Improved logic:
           Lets sum up the number of gaps plus sequence length of the first
           X sequences. if min == max it is aligned.
        */
        min_len = INT32_MAX;
        max_len = 0;
        gaps = 0;
        /* n = MACRO_MIN(50, msa->numseq); */
        n = msa->numseq;
        for(i = 0; i < n;i++){
                l = 0;
                for (j = 0; j <= msa->sequences[i]->len;j++){
                        l += msa->sequences[i]->gaps[j];
                }
                gaps += l;

                l += msa->sequences[i]->len;
                min_len = MACRO_MIN(min_len, l);
                max_len = MACRO_MAX(max_len, l);
                /* LOG_MSG("%d %d", max_len, min_len); */
        }
        if(gaps){
                if(min_len == max_len){ /* sequences have gaps and total length is identical - clearly aligned  */
                        msa->aligned = ALN_STATUS_ALIGNED;
                }else{          /* odd there are gaps but total length differs - unknown status  */
                        if(!msa->quiet){
                                aln_unknown_warning_message_gaps_but_len_diff(msa);
                        }
                        msa->aligned = ALN_STATUS_UNKNOWN;
                }
        }else{
                if(min_len == max_len){ /* no gaps and sequences have same length. Can' tell if they are aligned  */
                        if(!msa->quiet){
                                aln_unknown_warning_message_same_len_no_gaps();
                        }
                        msa->aligned = ALN_STATUS_UNKNOWN;
                }else{          /* No gaps and sequences have different lengths - unaligned */
                        msa->aligned = ALN_STATUS_UNALIGNED;
                }
        }
        //LOG_MSG("Aligned: %d gaps: %d",msa->aligned,gaps);
        return OK;
}

int set_sip_nsip(struct msa* msa)
{
        int i;
        ASSERT(msa!= NULL, "No msa");
        if(msa->plen){
                for (i = msa->num_profiles;i--;){
                        if(msa->sip[i]){
                                MFREE(msa->sip[i]);
                        }
                }
                if(msa->plen){
                        MFREE(msa->plen);
                }
                if(msa->sip){
                        MFREE(msa->sip);
                }
                if(msa->nsip){
                        MFREE(msa->nsip);
                }
                msa->plen = NULL;
                msa->sip = NULL;
                msa->nsip = NULL;
        }

        msa->num_profiles = (msa->numseq << 1 )-1;

        MMALLOC(msa->sip,sizeof(int*)* msa->num_profiles);
        MMALLOC(msa->nsip,sizeof(int)* msa->num_profiles);
        MMALLOC(msa->plen,sizeof(int)* msa->num_profiles);


        for (i =0;i < msa->num_profiles;i++){
                msa->sip[i] = NULL;
                msa->nsip[i] = 0;

        }

        for(i = 0;i < msa->numseq;i++){

                MMALLOC(msa->sip[i],sizeof(int));
                msa->nsip[i] = 1;
                msa->sip[i][0] = i;
                msa->plen[i] = 0;
        }
        return OK;
ERROR:
        return FAIL;
}



int convert_msa_to_internal(struct msa* msa, int type)
{
        struct alphabet* a = NULL;
        struct msa_seq* seq = NULL;
        int8_t* t = NULL;
        int i,j;

        RUNP(a = create_alphabet(type));

        t = a->to_internal;
        msa->L = a->L;
        for(i = 0; i <  msa->numseq;i++){
                seq = msa->sequences[i];
                for(j =0 ; j < seq->len;j++){
                        if(t[(int) seq->seq[j]] == -1){
                                WARNING_MSG("there should be no character not matching the alphabet");
                                WARNING_MSG("offending character: >>>%c<<<", seq->seq[j]);
                                /* exit(0); */
                        }else{
                                seq->s[j] = t[(int) seq->seq[j]];
                        }
                }

        }
        MFREE(a);
        return OK;
ERROR:
        if(a){
                MFREE(a);
        }
        return FAIL;
}

static int aln_unknown_warning_message_gaps_but_len_diff(struct msa* msa)
{
        int i;
        WARNING_MSG("--------------------------------------------");
        WARNING_MSG("The input sequences contain gap characters: ");

        for(i = 0; i < 128;i++){
                if(msa->letter_freq[i] && ispunct(i)){
                         WARNING_MSG("\"%c\" : %4d found                            ", (char)i,msa->letter_freq[i] );
                }
        }

        WARNING_MSG("BUT the presumably aligned sequences do not have the length.");
        WARNING_MSG("                                            ");
        WARNING_MSG("Kalign will remove the gap characters and   ");
        WARNING_MSG("align the sequences.                        ");
        WARNING_MSG("--------------------------------------------");
        return OK;
}

static int aln_unknown_warning_message_same_len_no_gaps(void)
{
        /* int i; */
        WARNING_MSG("--------------------------------------------");
        WARNING_MSG("All input sequences have the same length.   ");
        WARNING_MSG("BUT there are no gap characters.            ");
        WARNING_MSG("                                            ");
        WARNING_MSG("Unable to determine whether the sequences   ");
        WARNING_MSG("are already aligned.                        ");
        WARNING_MSG("Kalign will align the sequences.            ");
        WARNING_MSG("--------------------------------------------");
        return OK;
}

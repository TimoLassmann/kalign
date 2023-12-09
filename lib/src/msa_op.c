#include "tldevel.h"
#include "msa_struct.h"
#include "alphabet.h"

#include <ctype.h>

#include "msa_alloc.h"
#define MSA_OP_IMPORT
#include "msa_op.h"

static int aln_unknown_warning_message_gaps_but_len_diff(struct msa *msa);
static int aln_unknown_warning_message_same_len_no_gaps(void);

int msa_cpy(struct msa** dest, struct msa* src)
{
        int i;
        struct msa* d = NULL;
        d = *dest;
        if(d == NULL){
                RUN(alloc_msa(&d,src->alloc_numseq));
                /* d = alloc_msa(); */
        }
        if(d->biotype != ALN_BIOTYPE_UNDEF){
        /* if(d->L != ALPHA_UNDEFINED){ */
                if(d->biotype != src->biotype){
                        ERROR_MSG("Input alignments have different alphabets");
                }
        }
        if(d->aligned != 0 && d->aligned != ALN_STATUS_UNKNOWN){
                if(d->aligned != src->aligned){
                        d->aligned = ALN_STATUS_UNKNOWN;
                }
        }

        for(i = 0; i < 128;i++){
                d->letter_freq[i] += src->letter_freq[i];
        }

        d->numseq = 0;
        for(i = 0; i < src->numseq;i++){

                msa_seq_cpy(d->sequences[i], src->sequences[i]);
        }
        d->numseq = src->numseq;
        d->quiet = src->quiet;
        RUN(detect_alphabet(d));
        RUN(detect_aligned(d));
        RUN(set_sip_nsip(d));

        *dest = d;
        return OK;
ERROR:
        return FAIL;
}

int msa_seq_cpy(struct msa_seq *d, struct msa_seq *src)
{
        ASSERT(d != NULL,"No sequence");
        ASSERT(src != NULL,"No sequence");
        while(src->alloc_len > d->alloc_len){
                resize_msa_seq(d);
        }
        snprintf(d->name, MSA_NAME_LEN, "%s", src->name);

        for(int j = 0; j < src->len;j++){
                d->seq[j] = src->seq[j];
                d->s[j] = src->s[j];
                d->gaps[j] = src->gaps[j];
        }
        d->gaps[src->alloc_len] = src->gaps[src->alloc_len];
        d->seq[src->len] = 0;
        d->len = src->len;
        d->rank = src->rank;

        return OK;
ERROR:
        return FAIL;
}


int merge_msa(struct msa** dest, struct msa* src)
{
        int i;
        struct msa* d = NULL;
        d = *dest;
        if(d == NULL){
                RUN(alloc_msa(&d,src->alloc_numseq));
                /* d = alloc_msa(); */
        }
        if(d->biotype != ALN_BIOTYPE_UNDEF){
        /* if(d->L != ALPHA_UNDEFINED){ */
                if(d->biotype != src->biotype){
                        ERROR_MSG("Input alignments have different alphabets");
                }
        }
        if(d->aligned != 0 && d->aligned != ALN_STATUS_UNKNOWN){
                if(d->aligned != src->aligned){
                        d->aligned = ALN_STATUS_UNKNOWN;
                }
        }

        for(i = 0; i < 128;i++){
                d->letter_freq[i] += src->letter_freq[i];
        }

        for(i = 0; i < src->numseq;i++){
                free_msa_seq(d->sequences[d->numseq]);
                d->sequences[d->numseq] = src->sequences[i];
                src->sequences[i] = NULL;
                d->numseq++;
                if(d->alloc_numseq == d->numseq){
                        RUN(resize_msa(d));
                }
        }
        RUN(detect_alphabet(d));
        RUN(detect_aligned(d));
        RUN(set_sip_nsip(d));

        *dest = d;
        return OK;
ERROR:
        return FAIL;
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
                        msa->biotype = ALN_BIOTYPE_DNA;
                        /* msa->L = ALPHA_defDNA; */
                        /* RUN(convert_msa_to_internal(msa, ALPHA_defDNA)); */
                }else if(prot_prob > dna_prob){
                        if(!msa->quiet){
                                LOG_MSG("Detected protein sequences.");
                        }
                        msa->biotype = ALN_BIOTYPE_PROTEIN;
                        /* msa->L = ALPHA_redPROTEIN; */
                        /* RUN(convert_msa_to_internal(msa, ALPHA_redPROTEIN)); */
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

        }
        /* LOG_MSG("%d %d", max_len, min_len); */
        /* exit(0); */
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
        /* LOG_MSG("Aligned: %d gaps: %d",msa->aligned,gaps); */
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

int reformat_settings_msa(struct msa *msa, int rename, int unalign)
{
        for (int i = 0 ;i < msa->numseq;i++){
                        msa->nsip[i] = i;
        }
        if(rename){
                for (int i = 0 ;i < msa->numseq;i++){
                        snprintf(msa->sequences[i]->name, 128, "SEQ%d", i+1);
                }
        }
        if(unalign){
                RUN(dealign_msa(msa));
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

int kalign_msa_to_arr(struct msa* msa, char ***aligned, int *out_aln_len)
{
        ASSERT(msa != NULL,"No MSA!");
        ASSERT(msa->aligned == ALN_STATUS_FINAL,"Sequences are not finalized");

        char** out = NULL;
        int numseq = msa->numseq;
        int aln_len = msa->alnlen;
        MMALLOC(out, sizeof(char*) * numseq);
        for(int i = 0 ; i < numseq;i++){
                out[i] = NULL;
                MMALLOC(out[i], sizeof(char) * (aln_len +1));
        }

        /* galloc(&out, numseq,aln_len+1); */
        for(int i = 0 ; i < numseq;i++){
                for(int j = 0; j < aln_len;j++){
                        out[i][j] = msa->sequences[i]->seq[j];
                }
                out[i][aln_len] = 0;
                /* fprintf(stdout,"IN msa to ALIGNED: %s    %d\n", out[i], msa->alnlen); */
        }

        /* /\* msa->alnlen = aln_len; *\/ */
        /* /\* msa->aligned = ALN_STATUS_FINAL; *\/ */
        /* /\* int aln_len = 0; *\/ */



        /* for (int j = 0; j <= msa->sequences[0]->len;j++){ */
        /*         aln_len += msa->sequences[0]->gaps[j]; */
        /* } */
        /* aln_len += msa->sequences[0]->len; */
        /* aln_len += 1; */

        /* MMALLOC(out, sizeof(char*) * numseq); */
        /* for(int i = 0; i < numseq; i++){ */
        /*         out[i] = 0; */
        /*         MMALLOC(out[i], sizeof(char) * (uint64_t)aln_len); */
        /*         int pos = 0; */
        /*         for(int j = 0;j < msa->sequences[i]->len;j++){ */
        /*                 for(int c = 0;c < msa->sequences[i]->gaps[j];c++){ */
        /*                         out[i][pos] = '-'; */
        /*                         pos++; */
        /*                 } */
        /*                 out[i][pos] = msa->sequences[i]->seq[j]; */
        /*                 pos++; */
        /*         } */
        /*         for(int c = 0;c < msa->sequences[i]->gaps[ msa->sequences[i]->len];c++){ */
        /*                 out[i][pos] = '-'; */
        /*                 pos++; */
        /*         } */
        /*         out[i][pos] = 0; */
        /* } */

        *aligned = out;
        *out_aln_len = aln_len;
        return OK;
ERROR:
        return FAIL;
}


int kalign_arr_to_msa(char** input_sequences, int* len, int numseq,struct msa** multiple_aln)
{
        struct msa* msa = NULL;

        MMALLOC(msa, sizeof(struct msa));
        msa->sequences = NULL;
        msa->alloc_numseq = numseq;
        msa->numseq = numseq;
        msa->num_profiles = 0;
        msa->L = ALPHA_UNDEFINED;
        msa->aligned = 0;
        msa->plen = NULL;
        msa->sip = NULL;
        msa->nsip = NULL;
        msa->quiet = 1;
        MMALLOC(msa->sequences, sizeof(struct msa_seq*) * msa->alloc_numseq);

        for(int i = 0; i < 128; i++){
                msa->letter_freq[i] = 0;
        }

        for(int i = 0; i < msa->alloc_numseq;i++){
                msa->sequences[i] = NULL;
                struct msa_seq* seq = NULL;

                MMALLOC(seq, sizeof(struct msa_seq));
                seq->name = NULL;
                seq->seq = NULL;
                seq->s = NULL;
                seq->gaps = NULL;
                seq->len = len[i];
                seq->alloc_len = len[i]+1;

                MMALLOC(seq->name, sizeof(char)* MSA_NAME_LEN);

                MMALLOC(seq->seq, sizeof(char) * seq->alloc_len);
                MMALLOC(seq->s, sizeof(uint8_t) * seq->alloc_len);
                MMALLOC(seq->gaps, sizeof(int) * (seq->alloc_len+1));
                for(int j = 0;j < seq->alloc_len+1;j++){

                        seq->gaps[j] = 0;

                }
                for(int j = 0; j < len[i];j++){
                        msa->letter_freq[(int)input_sequences[i][j]]++;
                        seq->seq[j] = input_sequences[i][j];
                }
                seq->seq[len[i]] = 0;
                msa->sequences[i] = seq;
                /* LOG_MSG("%s",msa->sequences[i]->seq); */
        }
        RUN(detect_alphabet(msa));
        RUN(detect_aligned(msa));
        RUN(set_sip_nsip(msa));
        *multiple_aln = msa;
        return OK;

ERROR:
        kalign_free_msa(msa);
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

        WARNING_MSG("BUT the presumably aligned sequences do not ");
        WARNING_MSG("have the same length.                       ");
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

int finalise_alignment(struct msa* msa)
{
        ASSERT(msa->aligned == ALN_STATUS_ALIGNED, "Sequences are not aligned");
        struct msa_seq* seq = NULL;
        char* linear_seq = NULL;
        int aln_len = 0;

        
        for(int i = 0; i <= msa->sequences[0]->len;i++){
                aln_len += msa->sequences[0]->gaps[i];
        }
        aln_len += msa->sequences[0]->len;

        for(int i = 0; i < msa->numseq;i++){
                MMALLOC(linear_seq, sizeof(char)* (aln_len+1));
                seq = msa->sequences[i];
                RUN(make_linear_sequence(seq,linear_seq));
                MFREE(seq->seq);
                seq->seq = linear_seq;
                /* seq->len = aln_len; */
                linear_seq = NULL;
        }
        msa->alnlen = aln_len;
        msa->aligned = ALN_STATUS_FINAL;
        return OK;
ERROR:
        return FAIL;
}

int make_linear_sequence(struct msa_seq* seq, char* linear_seq)
{
        int c,j,f;
        f = 0;
        for(j = 0;j < seq->len;j++){
                //LOG_MSG("%d %d",j,seq->gaps[j]);
                for(c = 0;c < seq->gaps[j];c++){
                        linear_seq[f] = '-';
                        f++;

                }
                //LOG_MSG("%d %d %d",j,f,seq->gaps[j]);
                linear_seq[f] = seq->seq[j];
                f++;
        }
        for(c = 0;c < seq->gaps[ seq->len];c++){
                //LOG_MSG("%d %d",j,seq->gaps[seq->len]);
                linear_seq[f] = '-';
                f++;
        }
        linear_seq[f] = 0;
        ///fprintf(stdout,"LINEAR:%s\n",linear_seq);
        return OK;
}

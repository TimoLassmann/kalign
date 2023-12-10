#include "tldevel.h"

#include "msa_struct.h"
#include "alphabet.h"

#define MSA_ALLOC_IMPORT
#include "msa_alloc.h"

int alloc_msa(struct msa** msa, int numseq)
{
        struct msa* m = NULL;
        int i;
        MMALLOC(m, sizeof(struct msa));
        m->sequences = NULL;
        m->alloc_numseq = numseq;
        m->numseq = 0;
        m->num_profiles = 0;
        m->L = ALPHA_UNDEFINED;
        m->biotype = ALN_BIOTYPE_UNDEF;
        m->aligned = 0;
        m->alnlen = 0;
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
        kalign_free_msa(m);
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


void kalign_free_msa(struct msa* msa)
{
        int i;
        if(msa){
                for(i = 0; i < msa->alloc_numseq;i++){
                        if(msa->sequences[i]){
                                free_msa_seq(msa->sequences[i]);
                        }
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
        seq->rank = 0;
        /* seq->name_len = 128; */
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

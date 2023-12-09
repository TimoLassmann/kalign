#include "tldevel.h"
#include "tlrng.h"

#include "msa_struct.h"
#include "msa_alloc.h"
#include "msa_op.h"

#include "alphabet.h"

#define  DSSIM_IMPORT
#include "dssim.h"

typedef enum {
        TMM = 0,
        TMI = 1,
        TMD = 2,
        TII = 3,
        TIM = 4,
        TDD = 5,
        TDM = 6
} trans_type;

typedef enum {
        HMM_MATCH,
        HMM_INSERT,
        HMM_DELETE
} state_type;

struct hmm_param {
        double match_err_p;
        double insert_err_p;
        double indel_p;
        int n_observed_seq;
        int n_sim_seq;
        int len;
        int seed;
        int dna;
};

struct hmm {
        double** match_emit;
        double** insert_emit;
        double** transition;
        struct rng_state* rng;
        int L;
        int len;
};

int hmm_emit_simple(struct hmm *hmm, char **out,int *len);
int sample_pick(double *p, int len, struct rng_state* rng, int *pick);
char emit_letter(int c);
char emit_letter_dna(int c);
/* int hmm_init(struct hmm **hmm, int len, int seed); */
int hmm_init(struct hmm **hmm, struct hmm_param *p);

int get_prior_emit(double **prior);
int get_prior_emit_dna(double** prior);
int get_prior_trans(double** prior);

int pick(double *p, int len, struct rng_state *rng, int *pick);
int hmm_print(struct hmm *m);
int hmm_alloc(struct hmm **hmm, int len, int seed, int dna);
void hmm_free(struct hmm *n);

int hmm_param_init(struct hmm_param **param);
void hmm_param_free(struct hmm_param *p);

int dssim_get_fasta(struct msa **msa, int n_seq, int n_obs, int dna,int len, int seed)
{
        struct msa* m = NULL;
        struct hmm* hmm = NULL;
        struct hmm_param* p = NULL;
        char* tmp_seq = NULL;
        int r_len;
        /* RUN(alloc_msa(&m)); */
        hmm_param_init(&p);
        p->n_sim_seq = n_seq;
        p->n_observed_seq = n_obs;
        p->len = len;
        p->seed = seed;
        p->dna = dna;
        if(n_seq > 100){
                p->indel_p = 0.02;
        }else{
                p->indel_p = 0.04;
        }
        
        RUN(hmm_init(&hmm, p));
        /* hmm_print(hmm); */
        MMALLOC(m, sizeof(struct msa));
        m->sequences = NULL;
        m->alloc_numseq = n_seq;
        m->numseq = n_seq;
        m->num_profiles = 0;
        m->L = ALPHA_UNDEFINED;
        m->aligned = 0;
        m->plen = NULL;
        m->sip = NULL;
        m->nsip = NULL;
        m->quiet = 1;
        MMALLOC(m->sequences, sizeof(struct msa_seq*) * m->alloc_numseq);

        for(int i = 0; i < 128; i++){
                m->letter_freq[i] = 0;
        }


        m->numseq = 0;
        for(int i = 0; i < p->n_sim_seq;i++){
                hmm_emit_simple(hmm, &tmp_seq, &r_len);
                /* LOG_MSG("%s", tmp_seq); */
                m->sequences[i] = NULL;
                struct msa_seq* seq = NULL;

                MMALLOC(seq, sizeof(struct msa_seq));
                seq->name = NULL;
                seq->seq = NULL;
                seq->s = NULL;
                seq->gaps = NULL;
                seq->len = r_len;
                seq->alloc_len = r_len+1;
                seq->rank = i;

                MMALLOC(seq->name, sizeof(char)* MSA_NAME_LEN);
                snprintf(seq->name, MSA_NAME_LEN,"%d",i+1);
                MMALLOC(seq->seq, sizeof(char) * seq->alloc_len);
                MMALLOC(seq->s, sizeof(uint8_t) * seq->alloc_len);
                MMALLOC(seq->gaps, sizeof(int) * (seq->alloc_len+1));
                for(int j = 0;j < seq->alloc_len+1;j++){
                        seq->gaps[j] = 0;
                }
                for(int j = 0; j < r_len;j++){
                        m->letter_freq[(int)tmp_seq[j]]++;
                        seq->seq[j] = tmp_seq[j];
                }
                MFREE(tmp_seq);

                seq->seq[r_len] = 0;
                m->sequences[i] = seq;


                m->numseq++;
        }
        /* tl_random_double(struct rng_state *rng) */
        hmm_param_free(p);
        hmm_free(hmm);

        RUN(detect_alphabet(m));
        /* LOG_MSG("%d biotype",m->biotype); */
        RUN(detect_aligned(m));
        /* LOG_MSG("%d aligned ",m->aligned); */
        RUN(set_sip_nsip(m));
        /* LOG_MSG("%d aligned ",m->aligned); */


        *msa = m;
        return OK;
ERROR:
        return FAIL;
}

/* int main(int argc, char *argv[]) */
/* { */
/*         struct msa* m = NULL; */

/*         dssim_get_fasta(&m, 20, 10, 250, 1); */
/*         for(int i= 0; i < m->numseq;i++){ */
/*                 fprintf(stdout,">Seq_%d\n%s\n",i, m->sequences[i]->seq); */
/*         } */

/*         kalign_free_msa(m); */

/*         return EXIT_SUCCESS; */
/* ERROR: */
/*         kalign_free_msa(m); */
/*         return EXIT_FAILURE; */
/* } */



int hmm_emit_simple(struct hmm *hmm, char **out,int *len)
{
        char* seq = NULL;
        int n_alloc = 256;
        int n = 0;
        int pos = 0;
        int pick;
        state_type state = HMM_INSERT;
        double r;
        double sum;

        state = tl_random_int(hmm->rng, 3);
        MMALLOC(seq, sizeof(char) * n_alloc);
        seq[0] = 0;
        while(pos+1 < hmm->len){
                /* transiation;  */
                switch (state) {
                case HMM_MATCH:
                        r = tl_random_double(hmm->rng);
                        sum = hmm->transition[pos][TMM];
                        if(r < sum){
                                state = HMM_MATCH;
                        }else{
                                sum += hmm->transition[pos][TMI];
                                if(r < sum){
                                        state = HMM_INSERT;
                                }else{
                                        sum += hmm->transition[pos][TMD];
                                        if(r < sum){
                                                state = HMM_DELETE;
                                        }
                                }
                        }
                        pos++;
                        break;
                case HMM_INSERT:
                        r = tl_random_double(hmm->rng);
                        sum = hmm->transition[pos][TII];
                        if(r < sum){
                                state = HMM_INSERT;
                        }else{
                                state = HMM_MATCH;
                                pos++;
                        }
                        break;
                case HMM_DELETE:
                        r = tl_random_double(hmm->rng);
                        sum = hmm->transition[pos][TDD];
                        if(r < sum){
                                state = HMM_DELETE;
                                pos++;
                        }else{
                                state = HMM_MATCH;
                                pos++;
                        }
                        break;
                }

                /* Emit */
                switch (state) {
                case HMM_MATCH:
                        sample_pick(hmm->match_emit[pos], 20, hmm->rng,&pick);
                        if(hmm->L == 20){
                                seq[n] = emit_letter(pick);
                        }else{
                                seq[n] = emit_letter_dna(pick);
                        }
                        n++;
                        if(n == n_alloc){
                                n_alloc = n_alloc + n_alloc / 2;
                                MREALLOC(seq, sizeof(char) * n_alloc);
                        }

                        break;
                case HMM_INSERT:
                        sample_pick(hmm->insert_emit[pos], 20, hmm->rng,&pick);
                        if(hmm->L == 20){
                                seq[n] = emit_letter(pick);
                        }else{
                                seq[n] = emit_letter_dna(pick);
                        }
                        n++;
                        if(n == n_alloc){
                                n_alloc = n_alloc + n_alloc / 2;
                                MREALLOC(seq, sizeof(char) * n_alloc);
                        }
                        break;
                case HMM_DELETE:
                        break;
                }
                seq[n] = 0;
                /* n++; */
                if(n == n_alloc){
                        n_alloc = n_alloc + n_alloc / 2;
                        MREALLOC(seq, sizeof(char) * n_alloc);
                }
                /* LOG_MSG("%d %s n:%d state:%d",pos,seq,n ,state); */
        }
        seq[n] = 0;

        *len = n;

        *out = seq;

        return OK;
ERROR:
        return FAIL;
}

int sample_pick(double *p, int len, struct rng_state* rng, int *pick)
{
        double sum = 0.0;
        double r = 0.0;
        int sel = -1;
        r =tl_random_double(rng);
        sum = 0.0;
        for(int i = 0; i < len;i++){
                sum += p[i];
                if(r < sum){
                        sel = i;
                        break;
                }

        }
        /* LOG_MSG("Selected: %d (%d) r:%f %f", sel, len, r, sum); */
        /* if(sel == -1){ */

        /* } */
        *pick = sel;
        return OK;
}

char emit_letter(int c)
{
        char alphabet[20] = "ACDEFGHIKLMNPQRSTVWY";
        return alphabet[c];
}

char emit_letter_dna(int c)
{
        char alphabet[4] = "ACGT";
        return alphabet[c];
}

int hmm_init(struct hmm **hmm, struct hmm_param* p)
{
        struct hmm* n = NULL;
        double* prior_e = NULL;
        double* prior_t = NULL;
        double sum = 0.0;
        double r = 0.0;
        int pick;

        if(p->dna == 0){
                RUN(get_prior_emit(&prior_e));
                RUN(get_prior_trans(&prior_t));
        }else{
                RUN(get_prior_emit_dna(&prior_e));
                RUN(get_prior_trans(&prior_t));
        }

        prior_t[TMM] = 1.0 - p->indel_p;
        prior_t[TMI] = p->indel_p / 2.0;
        prior_t[TMD] = p->indel_p / 2.0;

        prior_t[TMM] = prior_t[TMM] / (prior_t[TMM] + prior_t[TMI] + prior_t[TMD]);
        prior_t[TMI] = prior_t[TMI] / (prior_t[TMM] + prior_t[TMI] + prior_t[TMD]);
        prior_t[TMD] = prior_t[TMD] / (prior_t[TMM] + prior_t[TMI] + prior_t[TMD]);

        RUN(hmm_alloc(&n, p->len,p->seed,p->dna));

        /* Start state...  */
        /* n->transition[0][TMM] = 0.5; */
        /* n->transition[0][TMI] = 0.5; */
        /* n->transition[0][TMD] = 0.0; */

        for(int i = 0; i < p->len;i++){
                for(int j = 0; j < n->L;j++){
                        n->match_emit[i][j] = 0.0;
                        n->insert_emit[i][j] = 0.0;
                }

                /* match  */
                /* pick =tl_random_int(n->rng, 20); */

                sample_pick(prior_e, n->L, n->rng, &pick);
                for(int j = 0; j < p->n_observed_seq;j++){
                        r = tl_random_double(n->rng);
                        if(r < p->match_err_p){
                                int c = tl_random_int(n->rng, n->L);
                                n->match_emit[i][c] += 1.0;
                        }else{
                                n->match_emit[i][pick] += 1.0;
                        }
                }
                /* insert  */
                /* pick =tl_random_int(n->rng, 20); */
                sample_pick(prior_e, n->L, n->rng, &pick);

                for(int j = 0; j < p->n_observed_seq;j++){
                        r = tl_random_double(n->rng);
                        if(r < p->insert_err_p){
                                int c = tl_random_int(n->rng, n->L);
                                n->insert_emit[i][c] += 1.0;
                        }else{
                                n->insert_emit[i][pick] += 1.0;
                        }
                }
                for(int j = 0; j < n->L;j++){
                        n->match_emit[i][j] += prior_e[j];
                        n->insert_emit[i][j] += prior_e[j];
                }

                sum = 0.0;
                for(int j = 0; j < n->L;j++){
                        sum += n->match_emit[i][j];
                }
                for(int j = 0; j < n->L;j++){
                        n->match_emit[i][j]/=sum;
                }

                sum = 0.0;
                for(int j = 0; j < n->L;j++){
                        sum += n->insert_emit[i][j];
                }
                for(int j = 0; j < n->L;j++){
                        n->insert_emit[i][j]/=sum;

                }

                /* transitions  */
                for(int j = 0; j < 7;j++){
                        n->transition[i][j] = prior_t[j];
                }
        }

        gfree(prior_e);
        gfree(prior_t);
        *hmm = n;
        return OK;
ERROR:
        return FAIL;
}

int get_prior_emit(double** prior)
{
        double* p = NULL;
        double sum = 0.0;
        RUN(galloc(&p , 20));
        p[0]  = 0.075520;                     /* A */
        p[1]  = 0.016973;                     /* C */
        p[2]  = 0.053029;                     /* D */
        p[3]  = 0.063204;                     /* E */
        p[4]  = 0.040762;                     /* F */
        p[5]  = 0.068448;                     /* G */
        p[6]  = 0.022406;                     /* H */
        p[7]  = 0.057284;                     /* I */
        p[8]  = 0.059398;                     /* K */
        p[9]  = 0.093399;                     /* L */
        p[10] = 0.023569;                     /* M */
        p[11] = 0.045293;                     /* N */
        p[12] = 0.049262;                     /* P */
        p[13] = 0.040231;                     /* Q */
        p[14] = 0.051573;                     /* R */
        p[15] = 0.072214;                     /* S */
        p[16] = 0.057454;                     /* T */
        p[17] = 0.065252;                     /* V */
        p[18] = 0.012513;                     /* W */
        p[19] = 0.031985;                     /* Y */

        sum = 0.0;
        for(int i = 0; i < 20;i++){
                sum += p[i];
        }
        for(int i = 0; i < 20;i++){
                p[i]/= sum;
        }
        *prior = p;
        return OK;
ERROR:
        return FAIL;
}

int get_prior_emit_dna(double** prior)
{
        double* p = NULL;
        double sum = 0.0;
        RUN(galloc(&p , 4));
        p[0]  = 0.2;                     /* A */
        p[1]  = 0.3;                     /* C */
        p[2]  = 0.3;                     /* G */
        p[3]  = 0.2;                     /* T */

        sum = 0.0;
        for(int i = 0; i < 4;i++){
                sum += p[i];
        }
        for(int i = 0; i < 4;i++){
                p[i]/= sum;
        }
        *prior = p;
        return OK;
ERROR:
        return FAIL;
}


int get_prior_trans(double** prior)
{
        double* p = NULL;
        RUN(galloc(&p , 7));

        p[TMM] = 0.96;
        p[TMI] = 0.02;
        p[TMD] = 0.02;
        p[TII] = 0.50;
        p[TIM] = 0.50;

        p[TDD] = 0.50;
        p[TDM] = 0.50;

        *prior = p;
        return OK;
ERROR:
        return FAIL;
}

int hmm_print(struct hmm *m)
{
        ASSERT(m == NULL,"No hmm");
        LOG_MSG("Emission: MATCH ");
        for(int i = 0; i < m->len;i++){
                fprintf(stdout,"   %4d ", i);
                for(int j = 0; j < m->L;j++){
                        fprintf(stdout,"%0.2f ", m->match_emit[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        LOG_MSG("Emission: Insert  ");
        for(int i = 0; i < m->len;i++){
                fprintf(stdout,"   %4d ", i);
                for(int j = 0; j < m->L;j++){
                        fprintf(stdout,"%0.2f ", m->insert_emit[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        LOG_MSG("Transition:");
        for(int i = 0; i < m->len;i++){
                fprintf(stdout,"   %4d ", i);
                for(int j = 0; j < 7;j++){
                        fprintf(stdout,"%0.2f ", m->transition[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        return OK;
ERROR:
        return FAIL;
}


int hmm_alloc(struct hmm **hmm, int len, int seed, int dna)
{
        struct hmm* n = NULL;

        ASSERT(len >0 ,"len has to be > 0");
        MMALLOC(n, sizeof(struct hmm));
        n->match_emit = NULL;
        n->insert_emit = NULL;
        n->transition = NULL;
        n->rng = init_rng(seed);
        if(!dna){
                n->L = 20;
        }else{
                n->L = 4;
        }
        n->len = len;

        galloc(&n->match_emit,n->len, n->L);
        galloc(&n->insert_emit,n->len, n->L);
        galloc(&n->transition, n->len, 7);

        *hmm = n;
        return OK;
ERROR:
        return FAIL;
}


void hmm_free(struct hmm *n)
{
        if(n){
                gfree(n->match_emit);
                gfree(n->insert_emit);
                gfree(n->transition);
                free_rng(n->rng);
                MFREE(n);
        }

}

int hmm_param_init(struct hmm_param **param)
{
        struct hmm_param* p = NULL;
        MMALLOC(p, sizeof(struct hmm_param));
        p->insert_err_p = 0.25;
        p->match_err_p = 0.05;
        p->indel_p = 0.04;
        p->n_sim_seq = 20;
        p->n_observed_seq = 10;
        p->len = 250;
        p->seed = 42;
        p->dna = 0;
        *param = p;
        return OK;
ERROR:
        return FAIL;
}

void hmm_param_free(struct hmm_param *p)
{
        if(p){
                MFREE(p);
        }
}

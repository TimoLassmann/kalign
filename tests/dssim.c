#include "tldevel.h"
#include "tlrng.h"

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


struct hmm {
        double** match_emit;
        double** insert_emit;
        double** transition;
        struct rng_state* rng;
        int L;
        int len;
};

int hmm_emit_simple(struct hmm *hmm, char **out);
int sample_pick(double *p, int len, struct rng_state* rng, int *pick);
char emit_letter(int c);
int hmm_init(struct hmm **hmm, int len, int seed);

int get_prior_emit(double **prior);
int get_prior_trans(double** prior);

int pick(double *p, int len, struct rng_state *rng, int *pick);
int hmm_print(struct hmm *m);
int hmm_alloc(struct hmm **hmm, int len, int seed);
void hmm_free(struct hmm *n);

int main(int argc, char *argv[])
{
        /* LOG_MSG("Hmm init"); */
        struct hmm* hmm = NULL;
        char* seq = NULL;
        RUN(hmm_init(&hmm, 350, 0));
        for(int i = 0; i < 20;i++){
                hmm_emit_simple(hmm, &seq);
                fprintf(stdout,">Seq_%d\n%s\n",i, seq);
                MFREE(seq);
        }
        /* tl_random_double(struct rng_state *rng) */
        hmm_free(hmm);
        return EXIT_SUCCESS;
ERROR:
        hmm_free(hmm);
        return EXIT_FAILURE;
}

int hmm_emit_simple(struct hmm *hmm, char **out)
{
        char* seq = NULL;
        int n_alloc = 256;
        int n = 0;
        int pos = 0;
        int pick;
        state_type state = HMM_MATCH;
        double r;
        double sum;
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
                        /* pick(double *p, int len, struct rng_state *rng, int *pick) */
                        sample_pick(hmm->match_emit[pos], 20, hmm->rng,&pick);
                        seq[n] = emit_letter(pick);
                        n++;
                        if(n == n_alloc){
                                n_alloc = n_alloc + n_alloc / 2;
                                MREALLOC(seq, sizeof(char) * n_alloc);
                        }

                        break;
                case HMM_INSERT:
                        sample_pick(hmm->insert_emit[pos], 20, hmm->rng,&pick);
                        seq[n] = emit_letter(pick);
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
        /* LOG_MSG("Selected: %d (%d)", sel, len); */
        *pick = sel;
        return OK;
}

char emit_letter(int c)
{
        char alphabet[20] = "ACDEFGHIKLMNPQRSTVWY";
        return alphabet[c];
}

int hmm_init(struct hmm **hmm, int len, int seed)
{
        struct hmm* n = NULL;
        double* prior_e = NULL;
        double* prior_t = NULL;
        double sum = 0.0;
        double r = 0.0;
        int pick;

        RUN(get_prior_emit(&prior_e));
        RUN(get_prior_trans(&prior_t));


        RUN(hmm_alloc(&n, len,seed));

        /* Start state...  */
        n->transition[0][TMM] = 0.9;
        n->transition[0][TMI] = 0.05;
        n->transition[0][TMD] = 0.05;

        for(int j = 0;j < 20;j++){
                n->match_emit[0][j] = 0.0;
                n->insert_emit[0][j] = 0.0;
        }

        for(int i = 1; i < len;i++){
                for(int j = 0; j < 20;j++){
                        n->match_emit[0][j] = 0.0;
                        n->insert_emit[0][j] = 0.0;
                }

                /* match  */
                /* pick =tl_random_int(n->rng, 20); */

                sample_pick(prior_e, 20, n->rng, &pick);
                for(int j = 0; j < 10;j++){
                        r = tl_random_double(n->rng);
                        if(r < 0.05){
                                int c = tl_random_int(n->rng, 20);
                                n->match_emit[i][c] += 1.0;
                        }else{
                                n->match_emit[i][pick] += 1.0;
                        }
                }
                /* insert  */
                /* pick =tl_random_int(n->rng, 20); */
                sample_pick(prior_e, 20, n->rng, &pick);

                for(int j = 0; j < 10;j++){
                        r = tl_random_double(n->rng);
                        if(r < 0.25){
                                int c = tl_random_int(n->rng, 20);
                                n->insert_emit[i][c] += 1.0;
                        }else{
                                n->insert_emit[i][pick] += 1.0;
                        }
                }
                for(int j = 0; j < 20;j++){
                        n->match_emit[i][j] += prior_e[j];
                        n->insert_emit[i][j] += prior_e[j];
                }
                sum = 0.0;
                for(int j = 0; j < 20;j++){
                        sum += n->match_emit[i][j];
                }
                for(int j = 0; j < 20;j++){
                        n->match_emit[i][j]/=sum;
                }
                sum = 0.0;
                for(int j = 0; j < 20;j++){
                        sum += n->insert_emit[i][j];
                }
                for(int j = 0; j < 20;j++){
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
        galloc(&p , 20);
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

int get_prior_trans(double** prior)
{
        double* p = NULL;
        galloc(&p , 7);

        p[TMM] = 0.9;
        p[TMI] = 0.05;
        p[TMD] = 0.05;
        p[TII] = 0.60;
        p[TIM] = 0.40;

        p[TDD] = 0.70;
        p[TDM] = 0.30;

        *prior = p;
        return OK;
ERROR:
        return FAIL;
}

int hmm_print(struct hmm *m)
{

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


int hmm_alloc(struct hmm **hmm, int len, int seed)
{
        struct hmm* n = NULL;

        ASSERT(len >0 ,"len has to be > 0");
        MMALLOC(n, sizeof(struct hmm));
        n->match_emit = NULL;
        n->insert_emit = NULL;
        n->transition = NULL;
        n->rng = init_rng(seed);
        n->L = 20;
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

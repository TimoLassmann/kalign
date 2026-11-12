#include "tldevel.h"
#include "tlrng.h"

#define BPM_IMPORT
#include "bpm.h"

#include <stdint.h>
#include <string.h>
#include  <stdalign.h>
#include "alphabet.h"
#include "esl_stopwatch.h"

/* Functions needed for the unit test*/

uint8_t dyn_256_print(const uint8_t* t,const uint8_t* p,int n,int m);
int  mutate_seq(uint8_t* s, int len,int k,int L, struct rng_state* rng);

#include "kalign_simd.h"
#include <stdlib.h>

#define SIGMA 13
#define DIV_CEIL(a,b) (a == 0 ? 1 : a/b+(a%b == 0 ? 0 : 1))


struct bpm_struct {
        uint64_t** Peq;
        uint64_t* P;
        uint64_t* M;
        uint64_t HIGH_BIT;
        uint64_t ONE;
        int32_t W;
        int32_t* score;
        int k;
};


static int bpm_block2(const uint8_t *t, const uint8_t *p, int n, int m, int maxd);

static int bpm_precompute(const uint8_t *p, int m, struct bpm_struct **out);
static int bpm_init_block(struct bpm_struct *s, int b);
static inline int advanceBlock(struct bpm_struct *s, int b, uint8_t c, int hIn);

int bpm_block2(const uint8_t *t, const uint8_t *p, int n, int m, int maxd)
{
        struct bpm_struct* s = NULL;
        int w;
        int32_t w_bytes;
        int32_t b_max;
        uint64_t ONE;

        w_bytes = sizeof(uint64_t);
        w =  w_bytes * 8;

        b_max = DIV_CEIL(m,w);

        bpm_precompute(p, m, &s);

        ONE = s->ONE;

        w_bytes = sizeof(uint64_t);

        w =  w_bytes * 8;
        int y = DIV_CEIL(maxd, w) - 1; /* 256 is max edit distance  */
        for (int b = 0; b <= y; b++) {
                bpm_init_block(s,b);
                s->score[b] = (uint32_t)(b + 1) * w;
        }

        for (int i = 0; i < n+s->W; i++) {
                uint8_t c = (uint8_t)t[i];
                int carry = 0;

                for (int b = 0; b <= y; b++) {
                        carry = advanceBlock(s,b, c, carry);
                        s->score[b] += carry;
                }

                if ((s->score[y] - carry <= maxd) && (y < (b_max - 1)) && ((s->Peq[c][y + 1] & ONE) || (carry < 0))) {
                        y += 1;
                        bpm_init_block(s,y);

                        s->score[y] = s->score[y - 1] + w - carry + advanceBlock(s,y, c, carry);
                } else {
                        while (s->score[y] >= (maxd + w)) {
                                if (y == 0) break;
                                y -= 1;
                        }
                }
                if(s->score[y] < s->k){
                        /* LOG_MSG("%d", s->score[y]); */
                        s->k = s->score[y];
                }


                /* if (y == (b_max - 1) && s->score[y] <= maxd) { */
                /*         assert(i - s->W >= 0); */
                /*         /\* result->push_front(SearchResult(score[y], (uint32_t)(i - W))); *\/ */
                /* } */
        }
        int k = s->k;
        for (int c = 0; c < SIGMA; c++) {
                MFREE(s->Peq[c]);
        }
        MFREE(s->Peq);
        MFREE(s->score);
        MFREE(s->P);
        MFREE(s->M);
        MFREE(s);

        return k;
}

inline int advanceBlock(struct bpm_struct *s, int b, uint8_t c, int hIn)
{
        uint64_t Pv = s->P[b];
        uint64_t Mv = s->M[b];
        uint64_t Eq = s->Peq[c][b];
        uint64_t ONE = s->ONE;
        uint64_t HIGH_BIT = s->HIGH_BIT;
        uint64_t Xv, Xh;
        uint64_t Ph, Mh;

        int h_out = 0;

        Xv = Eq | Mv;
        if (hIn < 0) {
                Eq |= ONE;
        }
        Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

        Ph = Mv | (~(Xh | Pv));
        Mh = Pv & Xh;

        if (Ph & HIGH_BIT) {
                h_out += 1;
        }
        if (Mh & HIGH_BIT) {
                h_out -= 1;
        }

        Ph <<= 1;
        Mh <<= 1;

        if (hIn < 0) {
                Mh |= ONE;
        } else if (hIn > 0) {
                Ph |= ONE;
        }

        Pv = Mh | (~(Xv | Ph));
        Mv = Ph & Xv;

        s->P[b] = Pv;
        s->M[b] = Mv;

        return h_out;
}

int bpm_init_block(struct bpm_struct *s, int b)
{
        int w;
        int32_t w_bytes;
        w_bytes = sizeof(uint64_t);

        w =  w_bytes * 8;
        s->P[b] = (uint64_t) -1;//Ones
        s->M[b] = 0;
        s->score[b] = (uint32_t)(b + 1) * w;
        return OK;
}


int bpm_precompute(const uint8_t *p, int m, struct bpm_struct **out)
{
        struct bpm_struct* s = NULL;
        uint64_t** Peq = NULL;
        uint64_t* P = NULL;
        uint64_t* M = NULL;
        int32_t *score = NULL;
        int32_t w_bytes;
        int32_t b_max;
        uint64_t W;
        int w;

        uint64_t HIGH_BIT;
        uint64_t ONE = 1;
        /* m = 1024; *- > 16 */
        w_bytes = sizeof(uint64_t);
        w =  w_bytes * 8;

        b_max = DIV_CEIL(m,w);
        /* LOG_MSG("%d -> %d",m,b_max); */
        W = w * b_max - m;
        HIGH_BIT = ONE << (w - 1);
        MMALLOC(score,b_max * sizeof(int32_t));
        MMALLOC(P ,b_max * w_bytes);
        MMALLOC(M ,b_max * w_bytes);

        MMALLOC(Peq,SIGMA * sizeof(uint64_t *));

        for (int c = 0; c < SIGMA; c++) {
                Peq[c] = NULL;
                MMALLOC(Peq[c], b_max*w_bytes);
                /* Peq[c] = (uint64_t *) malloc(b_max*w_bytes); */

                memset(Peq[c], 0, b_max*w_bytes);
                for (int32_t block = 0; block < b_max; block++) {
                        uint64_t bitPos = (uint64_t) 1;
                        for (int32_t i = block * w; i < (block + 1) * w; ++i) {
                                // fill the remainder after the last block with 1 (>m matches anything)
                                if (i >= m || p[i] == c) {
                                        Peq[c][block] |= bitPos;
                                }
                                bitPos <<= 1;
                        }
                }
        }
        MMALLOC(s ,sizeof(struct bpm_struct));
        s->M = M;
        s->P = P;
        s->Peq = Peq;
        s->W = W;
        s->ONE = ONE;
        s->HIGH_BIT = HIGH_BIT;
        s->score = score;
        s->k = m;
        *out = s;
        return OK;
ERROR:
        return FAIL;
}


/* The actual test.  */
int bpm_test(void);
int bpm_equal_test(void);

int main(void)
{

        /* no-op, must be harmless and repeatable */
        set_broadcast_mask();
        set_broadcast_mask();
        LOG_MSG("SIMD path: %s", kalign_simd_name());
        RUN(bpm_test());
        RUN(bpm_equal_test());
        kalign_simd_force_scalar(1);
        LOG_MSG("SIMD path: %s", kalign_simd_name());
        RUN(bpm_test());
        RUN(bpm_equal_test());
        kalign_simd_force_scalar(0);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int bpm_test(void)
{

        /* idea: start with identical sequences add errors run bpm */
        struct alphabet* alphabet = NULL;
        struct rng_state* rng;
        //long int r;
        int len = 0;
        int i,j,c;
        uint8_t* a = NULL;
        uint8_t* b = NULL;

        int test_iter = 100;

        int dyn_score;
        int bpm_score;
        int bb_score;
        int bb_score2;


        int diff[4][4] = {
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };



        RUNP(rng = init_rng(0));

        RUNP(alphabet = create_alphabet(ALPHA_redPROTEIN));

        char seq[] = "MLLRNSFISQDENDDQLSPTQLGRAAIASSLPPEASLAIFEDLNSASRAIALDTELHMLYLVYFYKNSRAQIIQKIFKIYSIFILKKFKNLEPKFKKKISENITVHITNSIRKKQHFWHVTPINVSVWQECDWHHLFSIFSKLPSDHKRIAKLVGVSEKFILDQLQGRRNDKLLQIHIRFFSALALFDLISEMSIYEVSHKYRIPRGCLQTLQSQSATYAAMIVAFCLRLGWTYLKALLDGFATRLLFGVRSELSELVAIEGIDGQRARILHERGVTCLSHLSACDSSKLAHFLTLAVPYSSSNSNDGLGEWLFGEPRMRVDVAARTLKERARKVLIRRVQELGISVELPKFEENEENIQESCDSGLPDSCEGMEDELEEKENIVKMEEMTKSVTEMSLTDNTISFKSEDDLFKKEIKVEEDEVFIKKEIDEDEEEIVEETVIECLETSLLKLKASTDEVFLRRLSQTFSPIGRSRSILNNSLLEDSFDRPVPRSSIPILNFITPKRESPTPYFEDSFDRPIPGSLPISSSRRKSVLTNIANLDSSRRESINSNASDNNSFDVFVTPPTKSAKEEKRRIAVKHPRVGNIIYSPLTSSPVIKHPKLEINHFYLKDVCHDHNAWNLWTKSSTSTSSCSIRVSDDYTGIAIRTDAGNTFIPLLETFGGENSVPQFFKNISNFYVLKMPLKIFWKQPSPASKYFESFSKCIIPLNTRLEFLKTLAVTVEMYISSMEDAFLIFEKFGIKIFRLKVVRIAAYLNNVIDVEQEENSNFLPILMDRYSILDPEIRKTCSSSLHKAAVEVYSLKPIFEKMCCSGASLQLEMESCQTVLNIFYSGIVFDQALCNSFIYKIRKQIENLEENIWRLAYGKFNIHSSNEVANVLFYRLGLIYPETSGCKPKLRHLPTNKLILEQMNTQHPIVGKILEYRQIQHTLTQCLMPLAKFIGRIHCWFEMCTSTGRILTSVPNLQNVPKRISSDGMSARQLFIANSENLLIGADYKQLELRVLAHLSNDSNLVNLITSDRDLFEELSIQWNFPRDAVKQLCYGLIYGMGAKSLSELTRMSIEDAEKMLKAFFAMFPGVRSYINETKEKVCKEEPISTIIGRRTIIKASGIGEERARIERVAVNYTIQGSASEIFKTAIVDIESKIKEFGAQIVLTIHDEVLVECPEIHVAAASESIENCMQNALSHLLRVPMRVSMKTGRSWADLK";
        len = strlen(seq);
        if(len > 256){
                len = 256;
        }

        MMALLOC(a , sizeof(uint8_t) * len) ;
        MMALLOC(b , sizeof(uint8_t) * len) ;

        for(i = 0;i < len;i++){
                a[i] = alphabet->to_internal[(int)seq[i]];
                b[i] =a[i];
                /* fprintf(stdout,"%d %d\n", a[i],b[i]); */
        }
        fprintf(stdout,"LEN: %d\n", len);


        LOG_MSG("Testing correctness of serial bpm.");
//        len = 63;

        len = 63; //could fix;
        for(i = 0; i < len;i++){
                for (j =0 ; j < test_iter; j++){
                        RUN(mutate_seq(b,len,i,alphabet->L,rng));

                        dyn_score = dyn_256(a,b,len,len);
                        bpm_score = bpm(a,b,len,len);
                        bb_score = bpm_block2(a,b,len,len, len);
                        bb_score2 = bpm_block(a,b,len,len);

                        /* if(j == 0){ */
                                /* LOG_MSG("k:%d %d %d %d %d",i,dyn_score,bpm_score, bb_score , bb_score2); */
                        /* } */
                        if(dyn_score != bpm_score){
                                diff[0][1]++;
                        }
                        if(dyn_score != bb_score){
                                diff[0][2]++;
                        }
                        if(dyn_score != bb_score2){
                                diff[0][3]++;
                        }

                        if(bpm_score != bb_score){
                                diff[1][2]++;
                        }
                        if(bpm_score != bb_score2){
                                diff[1][3]++;
                        }

                        if(bb_score != bb_score2){
                                diff[2][3]++;
                        }

                        /* restore sequence b */
                        for(c = 0;c < len;c++){
                                b[c] = a[c];
                        }
                }
        }
        for(int i = 0; i < 4;i++){
                for(int j = 0; j < 4;j++){
                        fprintf(stdout,"%3d ",diff[i][j]);
                }
                fprintf(stdout, "\n");
        }
        fprintf(stdout, "\n");


        LOG_MSG("Testing correctness of AVX bpm.");
        len = strlen(seq);
        if(len > 255){
                len = 255;
        }
        for(int i = 0; i < 4;i++){
                for(int j = 0; j < 4;j++){
                        diff[i][j] = 0;
                }

        }


        for(i = 0; i < len;i++){
                for (j =0 ; j < test_iter; j++){
                        RUN(mutate_seq(b,len,i,alphabet->L,rng));
                        dyn_score = dyn_256(a,b,len,len);
                        bpm_score = bpm_256(a,b,len,len);


                        bb_score = bpm_block2(a,b,len,len, len);
                        bb_score2 = bpm_block(a,b,len,len);

                        /* if(j == 0){ */
                                /* LOG_MSG("k:%d %d %d %d %d",i,dyn_score,bpm_score, bb_score , bb_score2); */
                        /* } */

                        if(dyn_score != bpm_score){
                                diff[0][1]++;
                        }
                        if(dyn_score != bb_score){
                                diff[0][2]++;
                        }
                        if(dyn_score != bb_score2){
                                diff[0][3]++;
                        }

                        if(bpm_score != bb_score){
                                diff[1][2]++;
                        }
                        if(bpm_score != bb_score2){
                                diff[1][3]++;
                        }

                        if(bb_score != bb_score2){
                                diff[2][3]++;
                        }

                        /* restore sequence b */
                        for(c = 0;c < len;c++){
                                b[c] = a[c];
                        }
                }

        }
        for(int i = 0; i < 4;i++){
                for(int j = 0; j < 4;j++){
                        fprintf(stdout,"%3d ",diff[i][j]);
                }
                fprintf(stdout, "\n");
        }
        fprintf(stdout, "\n");


        len = strlen(seq);
        if(len > 255){
                len = 255;
        }




        double dyn_timing = 0.0;
        double bpm_timing = 0.0;
        double block_timing = 0.0;
        double block2_timing = 0.0;
        int timing_iter = 10000;
        fprintf(stdout,"Dyn\tAVX\tDYN/AVX\n");

        LOG_MSG("Timing DYN 256");
        DECLARE_TIMER(t);


        for(i = 0; i < 1;i++){
                RUN(mutate_seq(b,len,i,alphabet->L,rng));
                START_TIMER(t);
                for(j = 0; j < timing_iter;j++){
                        dyn_score = dyn_256(a,b,len,len);
                }
                STOP_TIMER(t);
                /* GET_TIMING(t); */
                dyn_timing += GET_USERTIME(t);

                START_TIMER(t);
                for(j = 0; j < timing_iter;j++){
                        bpm_score = bpm_256(a,b,len,len);
                }

                STOP_TIMER(t);
                /* GET_TIMING(t); */
                bpm_timing +=GET_USERTIME(t);

                START_TIMER(t);
                for(j = 0; j < timing_iter;j++){
                        bpm_score = bpm_block(a,b,len,len);
                }

                STOP_TIMER(t);
                /* GET_TIMING(t); */
                block2_timing +=GET_USERTIME(t);



                START_TIMER(t);
                for(j = 0; j < timing_iter;j++){
                        bpm_score = bpm_block2(a,b,len,len,len);
                }

                STOP_TIMER(t);
                /* GET_TIMING(t); */
                block_timing +=GET_USERTIME(t);

                for(j = 0;j < len;j++){
                        b[j] = a[j];
                }
                //dyn_timing = GET_TIMING(t);
        }
        LOG_MSG("%f  - dyn ", dyn_timing);
        LOG_MSG("%f  - bpm ", bpm_timing);
        LOG_MSG("%f  - bpm_block2 ", block2_timing);
        LOG_MSG("%f  - bpm_block ", block_timing);
        DESTROY_TIMER(t);
        MFREE(alphabet);
        MFREE(a);
        MFREE(b);
        MFREE(rng);
        return OK;
ERROR:
        return FAIL;
}

int mutate_seq(uint8_t* s, int len,int k,int L, struct rng_state* rng)
{
        int i,j;
        int r;

        for(i = 0; i < k;i++){
                r = tl_random_int(rng,len);

                j = r;
                r = tl_random_int(rng,L);
                s[j] = r;
        }
        return OK;

}



uint8_t dyn_256_print(const uint8_t* t,const uint8_t* p,int n,int m)
{
        uint8_t* prev = NULL;
        uint8_t* cur = NULL;

        uint8_t* tmp = NULL;
        int i,j,c;

        MMALLOC(prev, sizeof(uint8_t)* 257);
        MMALLOC(cur, sizeof(uint8_t)* 257);
        cur[0] = 0;
        fprintf(stdout,"%d ", cur[0]);
        for(j = 1; j <= m;j++){
                cur[j] = cur[j-1] +1;
                fprintf(stdout,"%d ", cur[j]);
        }
        fprintf(stdout,"\n");
        tmp  = cur;
        cur = prev;
        prev = tmp;

        for(i = 1; i <= n;i++){
                cur[0] = prev[0];
                fprintf(stdout,"%d ", cur[0]);
                for(j = 1; j < m;j++){
                        c = 1;
                        if(t[i-1] == p[j-1]){
                                c = 0;

                        }
                        cur[j] = prev[j-1] +c ;

                        cur[j] = MACRO_MIN(cur[j], prev[j]+1);
                        cur[j] = MACRO_MIN(cur[j], cur[j-1]+1);
                        fprintf(stdout,"%d ", cur[j]);
                }
                j =m;
                c = 1;
                if(t[i-1] == p[j-1]){
                        c = 0;

                }
                cur[j] = prev[j-1] +c ;

                cur[j] = MACRO_MIN(cur[j], prev[j]);
                cur[j] = MACRO_MIN(cur[j], cur[j-1]+1);



                fprintf(stdout,"%d ", cur[j]);

                fprintf(stdout,"\n");
                tmp  = cur;
                cur = prev;
                prev = tmp;

        }
        c = prev[m];
        MFREE(prev);
        MFREE(cur);
        return c;
ERROR:
        return 255;

}


/* Plain and dispatched kernels must agree on every input: random pairs of
   random and related sequences, many lengths including the block edges. */
int bpm_equal_test(void)
{
        struct rng_state* rng = NULL;
        static const int lens[] = {0,1,2,31,32,33,63,64,65,100,127,128,129,191,192,193,
                255,256,257,300,383,384,385,511,512,513,700,1023,1024,1025,1500,2000};
        const int nl = (int)(sizeof(lens) / sizeof(lens[0]));
        uint8_t* a = NULL;
        uint8_t* b = NULL;
        long checked = 0;
        int li,lj,it,i;

        RUNP(rng = init_rng(42));
        MMALLOC(a, 2100);
        MMALLOC(b, 2100);
        for(li = 0; li < nl; li++){
                for(lj = 0; lj < nl; lj++){
                        /* all pairs for short lengths, a sample for long ones */
                        const int iters = (lens[li] > 700 || lens[lj] > 700) ? 2 : 6;
                        for(it = 0; it < iters; it++){
                                const int alpha = (it % 3 == 0) ? 2 : ((it % 3 == 1) ? 4 : 13);
                                const int n = lens[li];
                                const int m = lens[lj];
                                for(i = 0; i < n; i++){
                                        a[i] = (uint8_t)tl_random_int(rng, alpha);
                                }
                                if(it >= 3){
                                        /* pattern related to the text: copy with edits */
                                        for(i = 0; i < m; i++){
                                                if(i < n && tl_random_int(rng, 10) != 0){
                                                        b[i] = a[i];
                                                }else{
                                                        b[i] = (uint8_t)tl_random_int(rng, alpha);
                                                }
                                        }
                                }else{
                                        for(i = 0; i < m; i++){
                                                b[i] = (uint8_t)tl_random_int(rng, alpha);
                                        }
                                }
                                {
                                        const int r_plain = bpm_block_plain(a,b,n,m);
                                        const int r_disp = bpm_block(a,b,n,m);
                                        if(r_plain != r_disp){
                                                ERROR_MSG("bpm_block differs: n %d m %d plain %d dispatched %d",n,m,r_plain,r_disp);
                                        }
                                        if(n > 0 && m > 0){
                                                const int p256 = bpm_256_plain(a,b,n,m);
                                                const int d256 = bpm_256(a,b,n,m);
                                                if(p256 != d256){
                                                        ERROR_MSG("bpm_256 differs: n %d m %d plain %d dispatched %d",n,m,p256,d256);
                                                }
                                                if(m <= 255){
                                                        /* independent check against the DP */
                                                        const int dyn = dyn_256(a,b,n,m);
                                                        if(dyn != d256){
                                                                ERROR_MSG("bpm_256 vs dyn_256: n %d m %d dyn %d bpm %d",n,m,dyn,d256);
                                                        }
                                                }
                                        }else if(n == 0 || m == 0){
                                                /* degenerate: only require identical plain and dispatched (done above) */
                                        }
                                        checked++;
                                }
                        }
                }
        }
        LOG_MSG("bpm_equal_test: %ld pairs identical (%s)", checked, kalign_simd_name());
        MFREE(a);
        MFREE(b);
        free_rng(rng);
        return OK;
ERROR:
        return FAIL;
}

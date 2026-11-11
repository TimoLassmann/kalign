/* Hirschberg driver with pluggable kernels.
 *
 * Line-for-line copy of aln_runner_serial() + aln_continue() from
 * lib/src/aln_controller.c (Kalign 3.7.0), with the nine DP calls replaced by
 * function pointers.  The serial and the threaded aln_runner() perform the
 * same computations, so the path is identical to a real run's (verified by
 * bench_dp --selfcheck against paths recorded in the dumps).
 */
#include <float.h>
#include <stdint.h>

#include "kdp_kernel.h"

#ifdef KDP_HIRSCH_PROF
/* -DKDP_HIRSCH_PROF: time spent in the passes (by sub-call width) and in
   the meetup inside full alignments, printed at exit. */
#include <stdio.h>
#include <time.h>
static double hp_ns[4];
static double hp_wns[6]; static long long hp_wcells[6], hp_wcalls[6];   /* passes by width bucket */         /* 0: passes W >= 16, 1: passes W < 16, 2: meetup, 3: driver total */
static long long hp_calls[3], hp_cells[2];
static inline double hp_now(void)
{
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return ts.tv_sec * 1e9 + ts.tv_nsec;
}
__attribute__((destructor)) static void hp_print(void)
{
        if(hp_ns[3] == 0.0) return;
        fprintf(stderr, "hirsch_prof total %.3f s | passes W>=16 %.1f%% (%lld calls, %lld cells) | passes W<16 %.1f%% (%lld calls, %lld cells) | meetup %.1f%% (%lld) | driver rest %.1f%%\n",
                hp_ns[3] / 1e9, 100 * hp_ns[0] / hp_ns[3], hp_calls[0], hp_cells[0],
                100 * hp_ns[1] / hp_ns[3], hp_calls[1], hp_cells[1],
                100 * hp_ns[2] / hp_ns[3], hp_calls[2],
                100 * (hp_ns[3] - hp_ns[0] - hp_ns[1] - hp_ns[2]) / hp_ns[3]);
        const char* lab[6] = {"W<16", "16-31", "32-63", "64-127", "128-255", ">=256"};
        for(int b = 0; b < 6; b++){
                if(hp_wcells[b]) fprintf(stderr, "   W %-8s %5.1f%% of time  %5.1f%% of cells  %.3f ns/cell  %lld calls\n", lab[b],
                                         100 * hp_wns[b] / hp_ns[3], 100.0 * hp_wcells[b] / (hp_cells[0] + hp_cells[1]),
                                         hp_wns[b] / hp_wcells[b], hp_wcalls[b]);
        }
}
#endif

struct hctx {
        kdp_pass_fn fwd;
        kdp_pass_fn bwd;
        kdp_meet_fn meet;
};

static int runner(const struct hctx* h, struct aln_mem* m);

static int cont(const struct hctx* h, struct aln_mem* m, float input_states[], int old_cor[], int meet, int transition)
{
        int* path = m->path;
        switch(transition){
        case 1:
                path[old_cor[4]] = meet;
                path[old_cor[4] + 1] = meet + 1;
                m->f[0].a = input_states[0]; m->f[0].ga = input_states[1]; m->f[0].gb = input_states[2];
                m->b[0].a = 0.0F; m->b[0].ga = -FLT_MAX; m->b[0].gb = -FLT_MAX;
                m->starta = old_cor[0]; m->enda = old_cor[4] - 1;
                m->startb = old_cor[2]; m->endb = meet - 1;
                runner(h, m);
                m->starta = old_cor[4] + 1; m->enda = old_cor[1];
                m->startb = meet + 1; m->endb = old_cor[3];
                m->f[0].a = 0.0F; m->f[0].ga = -FLT_MAX; m->f[0].gb = -FLT_MAX;
                m->b[0].a = input_states[3]; m->b[0].ga = input_states[4]; m->b[0].gb = input_states[5];
                runner(h, m);
                break;
        case 2:
                path[old_cor[4]] = meet;
                m->f[0].a = input_states[0]; m->f[0].ga = input_states[1]; m->f[0].gb = input_states[2];
                m->b[0].a = 0.0F; m->b[0].ga = -FLT_MAX; m->b[0].gb = -FLT_MAX;
                m->starta = old_cor[0]; m->enda = old_cor[4] - 1;
                m->startb = old_cor[2]; m->endb = meet - 1;
                runner(h, m);
                m->starta = old_cor[4]; m->enda = old_cor[1];
                m->startb = meet + 1; m->endb = old_cor[3];
                m->f[0].a = -FLT_MAX; m->f[0].ga = 0.0F; m->f[0].gb = -FLT_MAX;
                m->b[0].a = input_states[3]; m->b[0].ga = input_states[4]; m->b[0].gb = input_states[5];
                runner(h, m);
                break;
        case 3:
                path[old_cor[4]] = meet;
                m->f[0].a = input_states[0]; m->f[0].ga = input_states[1]; m->f[0].gb = input_states[2];
                m->b[0].a = 0.0F; m->b[0].ga = -FLT_MAX; m->b[0].gb = -FLT_MAX;
                m->starta = old_cor[0]; m->enda = old_cor[4] - 1;
                m->startb = old_cor[2]; m->endb = meet - 1;
                runner(h, m);
                m->starta = old_cor[4] + 1; m->enda = old_cor[1];
                m->startb = meet; m->endb = old_cor[3];
                m->f[0].a = -FLT_MAX; m->f[0].ga = -FLT_MAX; m->f[0].gb = 0.0;
                m->b[0].a = input_states[3]; m->b[0].ga = input_states[4]; m->b[0].gb = input_states[5];
                runner(h, m);
                break;
        case 5:
                path[old_cor[4] + 1] = meet + 1;
                m->f[0].a = input_states[0]; m->f[0].ga = input_states[1]; m->f[0].gb = input_states[2];
                m->b[0].a = -FLT_MAX; m->b[0].ga = 0.0F; m->b[0].gb = -FLT_MAX;
                m->starta = old_cor[0]; m->enda = old_cor[4];
                m->startb = old_cor[2]; m->endb = meet - 1;
                runner(h, m);
                m->starta = old_cor[4] + 1; m->enda = old_cor[1];
                m->startb = meet + 1; m->endb = old_cor[3];
                m->f[0].a = 0.0F; m->f[0].ga = -FLT_MAX; m->f[0].gb = -FLT_MAX;
                m->b[0].a = input_states[3]; m->b[0].ga = input_states[4]; m->b[0].gb = input_states[5];
                runner(h, m);
                break;
        case 6:
                m->f[0].a = input_states[0]; m->f[0].ga = input_states[1]; m->f[0].gb = input_states[2];
                m->b[0].a = -FLT_MAX; m->b[0].ga = -FLT_MAX; m->b[0].gb = 0.0F;
                m->starta = old_cor[0]; m->enda = old_cor[4] - 1;
                m->startb = old_cor[2]; m->endb = meet;
                runner(h, m);
                m->starta = old_cor[4] + 1; m->enda = old_cor[1];
                m->startb = meet; m->endb = old_cor[3];
                m->f[0].a = -FLT_MAX; m->f[0].ga = -FLT_MAX; m->f[0].gb = 0.0F;
                m->b[0].a = input_states[3]; m->b[0].ga = input_states[4]; m->b[0].gb = input_states[5];
                runner(h, m);
                break;
        case 7:
                path[old_cor[4] + 1] = meet + 1;
                m->f[0].a = input_states[0]; m->f[0].ga = input_states[1]; m->f[0].gb = input_states[2];
                m->b[0].a = -FLT_MAX; m->b[0].ga = -FLT_MAX; m->b[0].gb = 0.0F;
                m->starta = old_cor[0]; m->enda = old_cor[4] - 1;
                m->startb = old_cor[2]; m->endb = meet;
                runner(h, m);
                m->starta = old_cor[4] + 1; m->enda = old_cor[1];
                m->startb = meet + 1; m->endb = old_cor[3];
                m->f[0].a = 0.0F; m->f[0].ga = -FLT_MAX; m->f[0].gb = -FLT_MAX;
                m->b[0].a = input_states[3]; m->b[0].ga = input_states[4]; m->b[0].gb = input_states[5];
                runner(h, m);
                break;
        default:
                break;
        }
        return 0;
}

static int runner(const struct hctx* h, struct aln_mem* m)
{
        float input_states[6];
        int old_cor[5];
        float score = 0.0F;
        int mid;
        int meet = 0;
        int transition = 0;

        if(m->starta >= m->enda){
                return 0;
        }
        if(m->startb >= m->endb){
                return 0;
        }
        input_states[0] = m->f[0].a;
        input_states[1] = m->f[0].ga;
        input_states[2] = m->f[0].gb;
        input_states[3] = m->b[0].a;
        input_states[4] = m->b[0].ga;
        input_states[5] = m->b[0].gb;

        mid = ((m->enda - m->starta) / 2) + m->starta;
        old_cor[0] = m->starta;
        old_cor[1] = m->enda;
        old_cor[2] = m->startb;
        old_cor[3] = m->endb;
        old_cor[4] = mid;

        m->enda = mid;
        m->starta_2 = mid;
        m->enda_2 = old_cor[1];

#ifdef KDP_HIRSCH_PROF
        {
                const int wide = (m->endb - m->startb - 1) >= 16;
                double t0 = hp_now();
                h->fwd(m);
                h->bwd(m);
                double t1 = hp_now();
                h->meet(m, old_cor, &meet, &transition, &score);
                double t2 = hp_now();
                hp_ns[wide ? 0 : 1] += t1 - t0;
                hp_calls[wide ? 0 : 1] += 2;
                hp_cells[wide ? 0 : 1] += (long long)(old_cor[1] - old_cor[0]) * (old_cor[3] - old_cor[2]);
                hp_ns[2] += t2 - t1;
                hp_calls[2]++;
                {
                        const int W = old_cor[3] - old_cor[2] - 1;
                        const int b = W < 16 ? 0 : W < 32 ? 1 : W < 64 ? 2 : W < 128 ? 3 : W < 256 ? 4 : 5;
                        hp_wns[b] += t1 - t0;
                        hp_wcells[b] += (long long)(old_cor[1] - old_cor[0]) * (old_cor[3] - old_cor[2]);
                        hp_wcalls[b] += 2;
                }
        }
#else
        h->fwd(m);
        h->bwd(m);
        h->meet(m, old_cor, &meet, &transition, &score);
#endif

        if(m->mode == ALN_MODE_SCORE_ONLY){
                m->score = score;
        }else{
                cont(h, m, input_states, old_cor, meet, transition);
        }
        return 0;
}

int kdp_hirschberg(const struct kdp_resolved* r, int kind, struct aln_mem* m)
{
        struct hctx h;
        h.fwd = r->fwd[kind];
        h.bwd = r->bwd[kind];
        h.meet = r->meet[kind];
#ifdef KDP_HIRSCH_PROF
        double t0 = hp_now();
        int ret = runner(&h, m);
        hp_ns[3] += hp_now() - t0;
        return ret;
#else
        return runner(&h, m);
#endif
}

int kdp_align(const struct kdp_resolved* r, const struct kdp_problem* p, struct aln_mem* m)
{
        kdp_mem_reset(p, m);
        if(r->k->align){
                return r->k->align(m, p->h.kind);
        }
        return kdp_hirschberg(r, p->h.kind, m);
}

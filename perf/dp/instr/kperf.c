/* kperf.c -- runtime of the KALIGN_PERF_TIMERS instrumentation and the
 * opt-in DP problem dump (KALIGN_DP_DUMP=<dir>).  Compiled only into the
 * instrumented library variant (kalign_perf_static, see perf/dp/CMakeLists.txt).
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <errno.h>
#include <sys/stat.h>

#include "aln_struct.h"
#include "aln_param.h"
#include "anchor_consistency.h"

#include "kperf.h"
#include "kdp_format.h"

#define KPERF_MAX_DEPTH 512
#define KPERF_MAX_RUNNER 64

struct kp_acc {
        uint64_t self_ns;
        uint64_t incl_ns;
        uint64_t calls;
        uint64_t cells;
};

static struct kp_acc acc[KPH_N][KP_N];

struct kp_frame {
        int id;
        int saved_phase;
        int outermost;
        uint64_t t0;
        uint64_t child;
        int64_t cells;
};

static _Thread_local struct kp_frame kp_stack[KPERF_MAX_DEPTH];
static _Thread_local int kp_depth = 0;
static _Thread_local int kp_phase = KPH_OTHER;
static _Thread_local int kp_id_active[KP_N];

static _Thread_local struct aln_mem* kp_runner_stack[KPERF_MAX_RUNNER];
static _Thread_local int kp_runner_depth = 0;

static const char* id_names[KP_N] = {
        "driver", "io_read", "io_write", "run_setup", "ensemble_other",
        "ens_poar", "ens_score", "ens_consensus",
        "tree_other", "tree_pick_anchor", "tree_dist_seq", "tree_dist_pairaln",
        "guide_ctx", "tree_kmeans",
        "cons_build", "cons_bonus",
        "progressive_other", "refine_other", "realign_dist", "finalise",
        "dp_runner",
        "dp_ss_fwd", "dp_ss_bwd", "dp_ss_meet",
        "dp_sp_fwd", "dp_sp_bwd", "dp_sp_meet",
        "dp_pp_fwd", "dp_pp_bwd", "dp_pp_meet",
        "dpc_ss_fwd", "dpc_ss_bwd", "dpc_ss_meet",
        "dpc_sp_fwd", "dpc_sp_bwd", "dpc_sp_meet",
        "dpc_pp_fwd", "dpc_pp_bwd", "dpc_pp_meet",
        "alnmem_init", "prof_make", "prof_gaps", "prof_update",
        "path_conv", "make_seq", "vsm", "sp_score"
};

static const char* phase_names[KPH_N] = {
        "other", "tree", "consbuild", "progressive", "refine", "realign", "ensemble"
};

const char* kperf_id_name(int id)
{
        if(id < 0 || id >= KP_N) return "?";
        return id_names[id];
}

const char* kperf_phase_name(int ph)
{
        if(ph < 0 || ph >= KPH_N) return "?";
        return phase_names[ph];
}

uint64_t kperf_now_ns(void)
{
#ifdef __APPLE__
        return clock_gettime_nsec_np(CLOCK_UPTIME_RAW);
#else
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
#endif
}

int kperf_push(int id, int phase, int64_t cells)
{
        struct kp_frame* f;
        if(kp_depth >= KPERF_MAX_DEPTH){
                return 0;
        }
        f = &kp_stack[kp_depth];
        f->id = id;
        f->saved_phase = kp_phase;
        f->child = 0;
        f->cells = cells;
        f->outermost = (kp_id_active[id] == 0);
        kp_id_active[id]++;
        if(phase >= 0){
                kp_phase = phase;
        }
        kp_depth++;
        f->t0 = kperf_now_ns();
        return kp_depth;
}

void kperf_pop_cleanup(int* tok)
{
        uint64_t t1 = kperf_now_ns();
        struct kp_frame* f;
        uint64_t elapsed;
        struct kp_acc* a;
        if(*tok <= 0){
                return;
        }
        /* Regions are strictly nested (scope based), so tok == kp_depth. */
        kp_depth = *tok - 1;
        f = &kp_stack[kp_depth];
        elapsed = t1 - f->t0;
        a = &acc[kp_phase][f->id];
        __atomic_fetch_add(&a->self_ns, elapsed - (f->child < elapsed ? f->child : elapsed), __ATOMIC_RELAXED);
        if(f->outermost){
                __atomic_fetch_add(&a->incl_ns, elapsed, __ATOMIC_RELAXED);
        }
        __atomic_fetch_add(&a->calls, 1, __ATOMIC_RELAXED);
        if(f->cells > 0){
                __atomic_fetch_add(&a->cells, (uint64_t)f->cells, __ATOMIC_RELAXED);
        }
        kp_id_active[f->id]--;
        kp_phase = f->saved_phase;
        if(kp_depth > 0){
                kp_stack[kp_depth - 1].child += elapsed;
        }
}

void kperf_reset(void)
{
        memset(acc, 0, sizeof(acc));
}

int kperf_write_json(FILE* f)
{
        int ph, id, first = 1;
        fprintf(f, "{\"records\": [\n");
        for(ph = 0; ph < KPH_N; ph++){
                for(id = 0; id < KP_N; id++){
                        struct kp_acc* a = &acc[ph][id];
                        if(a->calls == 0){
                                continue;
                        }
                        fprintf(f, "%s  {\"phase\": \"%s\", \"region\": \"%s\", \"self_ns\": %llu, "
                                "\"incl_ns\": %llu, \"calls\": %llu, \"cells\": %llu}",
                                first ? "" : ",\n",
                                phase_names[ph], id_names[id],
                                (unsigned long long)a->self_ns,
                                (unsigned long long)a->incl_ns,
                                (unsigned long long)a->calls,
                                (unsigned long long)a->cells);
                        first = 0;
                }
        }
        fprintf(f, "\n]}\n");
        return 0;
}

/* ------------------------------------------------------------------------ */
/* DP problem dump                                                           */
/* ------------------------------------------------------------------------ */

static int dump_init_done = 0;
static const char* dump_dir = NULL;
static const char* dump_tag = "";
static long dump_every = 1;
static long dump_max = 50;
static long long dump_min_cells = 0;
static uint64_t dump_seen[KPH_N][3];
static uint64_t dump_written[KPH_N][3];

static void dump_init(void)
{
        const char* s;
        if(__atomic_load_n(&dump_init_done, __ATOMIC_ACQUIRE)){
                return;
        }
        dump_dir = getenv("KALIGN_DP_DUMP");
        if(dump_dir && dump_dir[0] == 0){
                dump_dir = NULL;
        }
        if((s = getenv("KALIGN_DP_DUMP_TAG")) != NULL) dump_tag = s;
        if((s = getenv("KALIGN_DP_DUMP_EVERY")) != NULL){ dump_every = atol(s); if(dump_every < 1) dump_every = 1; }
        if((s = getenv("KALIGN_DP_DUMP_MAX")) != NULL){ dump_max = atol(s); }
        if((s = getenv("KALIGN_DP_DUMP_MIN_CELLS")) != NULL){ dump_min_cells = atoll(s); }
        if(dump_dir){
                if(mkdir(dump_dir, 0755) != 0 && errno != EEXIST){
                        fprintf(stderr, "kperf: cannot create KALIGN_DP_DUMP dir %s: %s (dumping disabled)\n",
                                dump_dir, strerror(errno));
                        dump_dir = NULL;
                }
        }
        __atomic_store_n(&dump_init_done, 1, __ATOMIC_RELEASE);
}

static int kind_of(const struct aln_mem* m)
{
        if(m->seq1) return KDP_SS;
        if(m->prof2) return KDP_PP;
        return KDP_SP;
}

static const char* kind_names[3] = { "ss", "sp", "pp" };

static void dump_write(const struct kperf_runner_tok* t)
{
        const struct aln_mem* m = t->m;
        struct kdp_header h;
        char path[4096];
        FILE* f;
        int i, j;

        memset(&h, 0, sizeof(h));
        memcpy(h.magic, KDP_MAGIC, 8);
        h.version = KDP_VERSION;
        h.kind = t->kind;
        h.phase = t->phase;
        h.mode = t->mode;
        h.len_a = m->len_a;
        h.len_b = m->len_b;
        h.starta = t->sa;
        h.enda = t->ea;
        h.startb = t->sb;
        h.endb = t->eb;
        h.sip = m->sip;
        h.consistency = m->consistency ? 1 : 0;
        h.flip = m->flip_threshold > 0.0F ? 1 : 0;
        h.has_path = (t->mode == ALN_MODE_FULL) ? 1 : 0;
        h.gpo = m->ap->gpo;
        h.gpe = m->ap->gpe;
        h.tgpe = m->ap->tgpe;
        h.subm_offset = m->ap->subm_offset;
        for(i = 0; i < 3; i++){
                h.f0[i] = t->f0[i];
                h.b0[i] = t->b0[i];
        }
        h.score = (t->mode == ALN_MODE_SCORE_ONLY) ? m->score : 0.0F;
        for(i = 0; i < 23; i++){
                for(j = 0; j < 23; j++){
                        h.subm[i * 23 + j] = m->ap->subm[i][j];
                }
        }
        h.call_index = t->idx;
        if(m->consistency){
                h.bonus_rows = m->consistency->n_rows;
                h.bonus_K = m->consistency->K;
        }

        snprintf(path, sizeof(path), "%s/%s%s_%s_%06llu_%dx%d.kdp",
                 dump_dir, dump_tag, phase_names[t->phase], kind_names[t->kind],
                 (unsigned long long)t->idx, m->len_a, m->len_b);
        f = fopen(path, "wb");
        if(!f){
                fprintf(stderr, "kperf: cannot write %s: %s\n", path, strerror(errno));
                return;
        }
        fwrite(&h, sizeof(h), 1, f);
        if(t->kind == KDP_SS){
                fwrite(m->seq1, 1, (size_t)m->len_a, f);
        }
        if(t->kind == KDP_SS || t->kind == KDP_SP){
                fwrite(m->seq2, 1, (size_t)m->len_b, f);
        }
        if(t->kind == KDP_SP || t->kind == KDP_PP){
                fwrite(m->prof1, sizeof(float), (size_t)(m->len_a + 2) * 64, f);
        }
        if(t->kind == KDP_PP){
                fwrite(m->prof2, sizeof(float), (size_t)(m->len_b + 2) * 64, f);
        }
        if(h.has_path){
                for(i = 0; i < m->len_a + 2; i++){
                        int32_t v = (int32_t)m->path[i];
                        fwrite(&v, sizeof(v), 1, f);
                }
        }
        if(m->consistency){
                size_t nb = (size_t)h.bonus_rows * (size_t)h.bonus_K;
                for(i = 0; i < (int)nb; i++){
                        int32_t v = (int32_t)m->consistency->cols[i];
                        fwrite(&v, sizeof(v), 1, f);
                }
                fwrite(m->consistency->vals, sizeof(float), nb, f);
        }
        fclose(f);
}

struct kperf_runner_tok kperf_runner_enter(struct aln_mem* m)
{
        struct kperf_runner_tok t;
        memset(&t, 0, sizeof(t));
        t.m = m;
        if(kp_runner_depth > 0 && kp_runner_stack[kp_runner_depth - 1] == m){
                /* Recursive Hirschberg call (or aln_runner -> aln_runner_serial). */
                t.top = 0;
                if(kp_runner_depth < KPERF_MAX_RUNNER){
                        kp_runner_stack[kp_runner_depth] = m;
                }
                kp_runner_depth++;
                return t;
        }
        t.top = 1;
        if(kp_runner_depth < KPERF_MAX_RUNNER){
                kp_runner_stack[kp_runner_depth] = m;
        }
        kp_runner_depth++;

        t.kind = kind_of(m);
        t.phase = kp_phase;
        t.sa = m->starta;
        t.ea = m->enda;
        t.sb = m->startb;
        t.eb = m->endb;
        t.mode = m->mode;
        t.f0[0] = m->f[0].a; t.f0[1] = m->f[0].ga; t.f0[2] = m->f[0].gb;
        t.b0[0] = m->b[0].a; t.b0[1] = m->b[0].ga; t.b0[2] = m->b[0].gb;

        dump_init();
        if(dump_dir){
                uint64_t seen = __atomic_fetch_add(&dump_seen[t.phase][t.kind], 1, __ATOMIC_RELAXED);
                long long cells = (long long)(t.ea - t.sa) * (long long)(t.eb - t.sb);
                t.idx = seen;
                if((seen % (uint64_t)dump_every) == 0 && cells >= dump_min_cells){
                        uint64_t w = __atomic_fetch_add(&dump_written[t.phase][t.kind], 1, __ATOMIC_RELAXED);
                        if(dump_max < 0 || (long)w < dump_max){
                                t.dump = 1;
                        }
                }
        }
        t.tok = kperf_push(KP_DP_RUNNER, -1, (int64_t)(t.ea - t.sa) * (int64_t)(t.eb - t.sb));
        return t;
}

void kperf_runner_exit(struct kperf_runner_tok* t)
{
        kp_runner_depth--;
        if(!t->top){
                return;
        }
        kperf_pop_cleanup(&t->tok);
        if(t->dump){
                dump_write(t);
        }
}

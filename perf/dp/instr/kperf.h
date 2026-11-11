/* kperf.h -- compile-time gated timer instrumentation for Kalign.
 *
 * Only included when KALIGN_PERF_TIMERS is defined (see the guarded include
 * at the end of lib/src/tldevel.h, which otherwise defines every KPERF_*
 * macro as empty, so a normal build contains no instrumentation at all).
 *
 * Model: a per-thread stack of timed regions.  Each region records
 *   self  : time spent in the region excluding nested timed regions
 *   incl  : inclusive time (outermost activation of that id only)
 *   calls : number of activations
 *   cells : work units supplied by the region (DP cells for DP kernels)
 * keyed by (phase, region id).  The phase is the pipeline stage that is
 * active (guide tree, consistency build, progressive, refine, ...), so the
 * same DP kernel is accounted separately per stage.  Self times of all
 * regions under the driver's top-level region sum to its inclusive time.
 *
 * Regions are opened with a cleanup-attribute variable so every return path
 * closes them; no edits to return statements are needed.
 *
 * Threading: the stack and phase are thread local, the accumulators are
 * updated atomically.  Phases are not propagated to thread-pool workers, so
 * profile with a single thread (-n 1) for a clean per-phase breakdown.
 */
#ifndef KPERF_H
#define KPERF_H

#include <stdint.h>
#include <stdio.h>

struct aln_mem;

enum kperf_phase {
        KPH_OTHER = 0,
        KPH_TREE,          /* guide tree: anchors, distances, k-means/UPGMA */
        KPH_CONSBUILD,     /* anchor consistency table build (pairwise DPs) */
        KPH_PROGRESSIVE,   /* progressive alignment along the guide tree */
        KPH_REFINE,        /* two-pass refinement (refine_alignment) */
        KPH_REALIGN,       /* realign: distances from the alignment */
        KPH_ENSEMBLE,      /* ensemble bookkeeping (POAR, scoring, consensus) */
        KPH_N
};

enum kperf_id {
        KP_DRIVER = 0,          /* whole driver (perf/dp/kalign_prof) */
        KP_IO_READ,
        KP_IO_WRITE,
        KP_RUN,                 /* kalign_single_run self: input check, sort, alphabet, params */
        KP_ENSEMBLE,            /* kalign_ensemble_from_configs self: msa copies etc. */
        KP_ENS_POAR,            /* pos_matrix_from_msa + extract_poars */
        KP_ENS_SCORE,           /* score_alignments */
        KP_ENS_CONSENSUS,       /* build_consensus */
        KP_TREE,                /* build_tree_kmeans* self */
        KP_TREE_PICK_ANCHOR,    /* pick_anchor */
        KP_TREE_DIST_SEQ,       /* d_estimation: BPM / k-mer distances */
        KP_TREE_DIST_PAIRALN,   /* guide_dist_*_matrix self (DP counted separately) */
        KP_GUIDE_CTX,           /* guide_dist_ctx_create */
        KP_TREE_KMEANS,         /* bisecting_kmeans self (clustering, UPGMA) */
        KP_CONS_BUILD,          /* anchor_consistency_build self */
        KP_CONS_BONUS,          /* per-merge sparse bonus construction */
        KP_PROGRESSIVE,         /* create_msa_tree* self */
        KP_REFINE,              /* refine_alignment self */
        KP_REALIGN_DIST,        /* pick_anchor_from_alignment, compute_aln_anchor_dist */
        KP_FINALISE,            /* finalise_alignment */
        KP_DP_RUNNER,           /* top-level aln_runner self: Hirschberg recursion overhead */
        /* DP kernels: base + 9*consistency + 3*kind + op */
        KP_DP_SS_FWD, KP_DP_SS_BWD, KP_DP_SS_MEET,
        KP_DP_SP_FWD, KP_DP_SP_BWD, KP_DP_SP_MEET,
        KP_DP_PP_FWD, KP_DP_PP_BWD, KP_DP_PP_MEET,
        KP_DPC_SS_FWD, KP_DPC_SS_BWD, KP_DPC_SS_MEET,
        KP_DPC_SP_FWD, KP_DPC_SP_BWD, KP_DPC_SP_MEET,
        KP_DPC_PP_FWD, KP_DPC_PP_BWD, KP_DPC_PP_MEET,
        KP_ALNMEM,              /* init_alnmem */
        KP_PROF_MAKE,           /* make_profile_n */
        KP_PROF_GAPS,           /* set_gap_penalties_n */
        KP_PROF_UPDATE,         /* update_n */
        KP_PATH_CONV,           /* add_gap_info_to_path_n, mirror_path_n, convert_raw_path */
        KP_MAKE_SEQ,            /* make_seq (gap bookkeeping) */
        KP_VSM,                 /* vsm_merge_offset, vsm_(un)shift_profiles */
        KP_SP_SCORE,            /* compute_sp_score */
        KP_N
};

#define KPERF_DP_FWD 0
#define KPERF_DP_BWD 1
#define KPERF_DP_MEET 2
#define KPERF_DP_SS 0
#define KPERF_DP_SP 1
#define KPERF_DP_PP 2

struct kperf_runner_tok {
        struct aln_mem* m;
        int top;
        int tok;
        int dump;
        int kind;
        int phase;
        int sa, ea, sb, eb;
        int mode;
        float f0[3];
        float b0[3];
        uint64_t idx;
};

int  kperf_push(int id, int phase, int64_t cells);
void kperf_pop_cleanup(int* tok);
struct kperf_runner_tok kperf_runner_enter(struct aln_mem* m);
void kperf_runner_exit(struct kperf_runner_tok* t);

/* Reporting API (used by the driver). */
void kperf_reset(void);
int  kperf_write_json(FILE* f);
const char* kperf_id_name(int id);
const char* kperf_phase_name(int ph);
uint64_t kperf_now_ns(void);

#define KPERF_CAT2(a, b) a##b
#define KPERF_CAT(a, b) KPERF_CAT2(a, b)

/* Open a timed region that lasts until the enclosing scope exits. */
#define KPERF_REGION(id) \
        int KPERF_CAT(kperf_tok_, __LINE__) __attribute__((cleanup(kperf_pop_cleanup))) = kperf_push((id), -1, 0)
/* Same, and switch the phase for the duration of the region. */
#define KPERF_PHASE_REGION(id, ph) \
        int KPERF_CAT(kperf_tok_, __LINE__) __attribute__((cleanup(kperf_pop_cleanup))) = kperf_push((id), (ph), 0)
/* DP kernel region: kind KPERF_DP_SS/SP/PP, op KPERF_DP_FWD/BWD/MEET. */
#define KPERF_DP(kind, op, m, cells) \
        int KPERF_CAT(kperf_tok_, __LINE__) __attribute__((cleanup(kperf_pop_cleanup))) = \
                kperf_push(KP_DP_SS_FWD + ((m)->consistency ? 9 : 0) + 3 * (kind) + (op), -1, (int64_t)(cells))
/* Top-level Hirschberg call: timing of the recursion and the opt-in dump. */
#define KPERF_RUNNER(m) \
        struct kperf_runner_tok kperf_runner_tok_ __attribute__((cleanup(kperf_runner_exit))) = kperf_runner_enter(m)

#endif

/* kdp_kernel.h -- pluggable DP kernel interface for bench_dp.
 *
 * A kernel provides any subset of the nine Hirschberg building blocks
 * (forward pass, backward pass, meetup; for seq-seq, seq-profile and
 * profile-profile), with exactly the semantics and signatures of Kalign's
 * aln_<kind>_foward/_backward/_meetup in lib/src/aln_*.c.  Missing entries
 * fall back to the reference kernel ("ref"), so a prototype can replace
 * only, say, the profile-profile forward and backward passes.
 *
 * Alternatively a kernel may provide .align, a complete aligner that takes
 * the whole problem and writes the raw Hirschberg-style path into m->path
 * (path[i] = j (1-based) if row i is matched to column j, -1 if row i is
 * gapped, i = 1..len_a), bypassing the Hirschberg driver.
 *
 * To add a prototype: drop a .c file into perf/dp/kernels/ that defines a
 * struct kdp_kernel and registers it:
 *
 *     #include "kdp_kernel.h"
 *     static int my_pp_fwd(struct aln_mem* m) { ... }
 *     static const struct kdp_kernel k = {
 *             .name = "pp_neon",
 *             .description = "NEON profile-profile forward",
 *             .fwd = { [KDP_PP] = my_pp_fwd },
 *     };
 *     KDP_REGISTER_KERNEL(k)
 *
 * then re-run cmake / make bench_dp and use `bench_dp --kernel pp_neon`.
 *
 * Data available to a kernel (struct aln_mem, lib/src/aln_struct.h):
 *   m->prof1/m->prof2  profiles, 64 floats per column, rows 0..len+1
 *                      ([0..22] residue counts, [27..29] gap open/ext/terminal
 *                      for this column, [32..54] pre-weighted substitution
 *                      scores); m->seq1/m->seq2 encoded residues (0..22);
 *   m->ap->subm[23][23], m->ap->gpo/gpe/tgpe/subm_offset, m->sip;
 *   m->f / m->b        state rows (struct states {a, ga, gb}), size >= len+2;
 *   m->starta/enda/startb/endb, m->starta_2/enda_2, m->len_a/len_b.
 *   m->consistency     the sparse consistency bonus recorded in the dump
 *                      (struct sparse_bonus, anchor_consistency.h; add
 *                      sparse_bonus_lookup(m->consistency, i, j) to the
 *                      match score as lib/src/aln_*.c do), or NULL when the
 *                      real call had none or bench_dp runs with --no-bonus.
 * kdp_kernel_state is set to whatever .prepare returned (or NULL) before
 * any kernel function of that problem is called.
 */
#ifndef KDP_KERNEL_H
#define KDP_KERNEL_H

#include <stdint.h>
#include "kdp_format.h"
#include "aln_struct.h"
#include "aln_param.h"
#include "anchor_consistency.h"

typedef int (*kdp_pass_fn)(struct aln_mem* m);
typedef int (*kdp_meet_fn)(struct aln_mem* m, int old_cor[], int* meet, int* t, float* score);

struct kdp_kernel {
        const char* name;
        const char* description;
        kdp_pass_fn fwd[3];     /* indexed by enum kdp_kind: KDP_SS, KDP_SP, KDP_PP */
        kdp_pass_fn bwd[3];
        kdp_meet_fn meet[3];
        /* Optional complete aligner (FULL mode): fills m->path, see above. */
        int (*align)(struct aln_mem* m, int kind);
        /* Optional per-problem setup/teardown (e.g. transposed or packed
           profiles).  Timed separately (prepare_ns), not inside ns/cell. */
        int (*prepare)(struct aln_mem* m, int kind, void** state);
        void (*release)(void* state);
};

extern void* kdp_kernel_state;

void kdp_register_kernel(const struct kdp_kernel* k);
const struct kdp_kernel* kdp_find_kernel(const char* name);
int kdp_n_kernels(void);
const struct kdp_kernel* kdp_kernel_at(int i);

/* Resolved kernel: every slot filled (fallback to "ref"). */
struct kdp_resolved {
        const struct kdp_kernel* k;
        kdp_pass_fn fwd[3];
        kdp_pass_fn bwd[3];
        kdp_meet_fn meet[3];
};
int kdp_resolve(const struct kdp_kernel* k, struct kdp_resolved* r);

#define KDP_REGISTER_KERNEL(kvar) \
        __attribute__((constructor)) static void kdp_register_##kvar(void) { kdp_register_kernel(&(kvar)); }

/* ---- problems ---- */

struct kdp_problem {
        char path[1024];
        struct kdp_header h;
        uint8_t* seq1;
        uint8_t* seq2;
        float* prof1;
        float* prof2;
        int32_t* ref_path;      /* path from the real run (may be NULL) */
        struct aln_param ap;    /* subm points into subm_rows */
        float* subm_rows[23];
        struct sparse_bonus bonus;      /* valid if h.consistency */
        int use_bonus;                  /* attach bonus to m->consistency */
};

int  kdp_load(const char* path, struct kdp_problem** p);
void kdp_free(struct kdp_problem* p);

/* aln_mem set up for problem p; f/b/path buffers sized for p. */
int  kdp_mem_create(const struct kdp_problem* p, struct aln_mem** m);
void kdp_mem_free(struct aln_mem* m);
/* Reset m for a full-range call: starts/ends, f[0]/b[0] from the dump. */
void kdp_mem_reset(const struct kdp_problem* p, struct aln_mem* m);

/* Hirschberg driver using the resolved kernel (copy of aln_runner_serial +
   aln_continue from lib/src/aln_controller.c with function pointers). */
int kdp_hirschberg(const struct kdp_resolved* r, int kind, struct aln_mem* m);
/* Full alignment of p into m->path using r (driver or r->k->align). */
int kdp_align(const struct kdp_resolved* r, const struct kdp_problem* p, struct aln_mem* m);

/* Independent scorer: score of a raw path under Kalign's top-level affine
   model (terminal gaps [29] without open/close, internal gaps open [27],
   extend [28], gap->match close = [27] of the previous column), in double. */
double kdp_path_score(const struct kdp_problem* p, const int* path);

#endif

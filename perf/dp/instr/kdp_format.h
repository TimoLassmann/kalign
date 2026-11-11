/* kdp_format.h -- on-disk format of a dumped Kalign DP problem (.kdp).
 *
 * One file = one top-level call of aln_runner()/aln_runner_serial() (the
 * full Hirschberg alignment of two sequences/profiles), captured during a
 * real Kalign run by the instrumented library (KALIGN_DP_DUMP=<dir>).
 *
 * Layout (native endianness, little-endian on arm64/x86-64):
 *   struct kdp_header
 *   uint8_t seq1[len_a]                 if kind == KDP_SS
 *   uint8_t seq2[len_b]                 if kind == KDP_SS or KDP_SP
 *   float   prof1[(len_a + 2) * 64]     if kind == KDP_SP or KDP_PP
 *   float   prof2[(len_b + 2) * 64]     if kind == KDP_PP
 *   int32_t path[len_a + 2]             if has_path (raw Hirschberg path
 *                                       produced by the real run)
 *   int32_t bonus_cols[bonus_rows * bonus_K]  if consistency (the sparse
 *   float   bonus_vals[bonus_rows * bonus_K]  consistency bonus, struct
 *                                       sparse_bonus in anchor_consistency.h)
 *
 * Kinds follow Kalign's internal orientation after do_align() swapped the
 * inputs: "a" is always the DP row dimension (outer loop), "b" the column
 * dimension (inner loop).  For KDP_SP the profile is prof1 (rows) and the
 * sequence is seq2 (columns), exactly as in aln_seqprofile.c.
 */
#ifndef KDP_FORMAT_H
#define KDP_FORMAT_H

#include <stdint.h>

#define KDP_MAGIC "KALDP02"
#define KDP_VERSION 2

enum kdp_kind { KDP_SS = 0, KDP_SP = 1, KDP_PP = 2 };

struct kdp_header {
        char magic[8];          /* KDP_MAGIC, NUL terminated */
        int32_t version;
        int32_t kind;           /* enum kdp_kind */
        int32_t phase;          /* enum kperf_phase of the caller (tree, progressive, ...) */
        int32_t mode;           /* ALN_MODE_FULL (1) or ALN_MODE_SCORE_ONLY (2) */
        int32_t len_a;
        int32_t len_b;
        int32_t starta;         /* range of the top-level call (normally 0,len_a,0,len_b) */
        int32_t enda;
        int32_t startb;
        int32_t endb;
        int32_t sip;            /* m->sip (seq-profile gap scaling) */
        int32_t consistency;    /* 1 if the real run used a sparse consistency bonus
                                   (dumped after the path, see bonus_rows/bonus_K) */
        int32_t flip;           /* 1 if the real run used flip_threshold > 0 (inline refine trials) */
        int32_t has_path;       /* 1 if path[] follows */
        float gpo;              /* m->ap->gpo/gpe/tgpe/subm_offset at call time */
        float gpe;
        float tgpe;
        float subm_offset;
        float f0[3];            /* m->f[0] a, ga, gb at entry */
        float b0[3];            /* m->b[0] a, ga, gb at entry */
        float score;            /* m->score after the call (score-only mode), else 0 */
        float subm[23 * 23];    /* m->ap->subm (seq-seq kernel reads it) */
        uint64_t call_index;    /* per (phase,kind) top-level call counter */
        int32_t bonus_rows;     /* sparse_bonus n_rows (0 if no bonus) */
        int32_t bonus_K;        /* sparse_bonus K */
        int32_t reserved[6];
};

#endif

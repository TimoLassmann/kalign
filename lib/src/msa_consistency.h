#ifndef MSA_CONSISTENCY_H
#define MSA_CONSISTENCY_H

#ifdef MSA_CONSISTENCY_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

struct msa;
struct poar_table;
struct sparse_bonus;

/* Context for POAR-based consistency scoring.
 * Stored as a non-owning reference in msa->poar_consistency. */
struct poar_consistency_ctx {
        struct poar_table* poar;  /* non-owning reference */
        int n_runs;
        float weight;
};

/* Compute consistency bonus for a pair of tree nodes from POAR data.
 *
 * For each member sequence pair (sa in node_a, sb in node_b), looks up
 * the POAR entries to find how often each residue pair was co-aligned
 * across the N ensemble runs.  Maps ungapped positions to DP positions
 * using the current gap structure and accumulates bonuses into a
 * sparse_bonus matrix.
 *
 * node_a/node_b: tree node indices (leaf = sequence index)
 * len_a/len_b:   DP dimensions (profile widths)
 */
EXTERN int poar_consistency_get_bonus(struct poar_consistency_ctx* ctx,
                                       struct msa* msa,
                                       int node_a, int len_a,
                                       int node_b, int len_b,
                                       struct sparse_bonus** bonus_out);

#undef MSA_CONSISTENCY_IMPORT
#undef EXTERN

#endif

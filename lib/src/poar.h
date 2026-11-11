#ifndef POAR_H
#define POAR_H

#ifdef POAR_IMPORT
#define EXTERN
#else
#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif
#endif

#include <stdint.h>

/* Position matrix: maps (seq, col) -> residue position or -1 (gap) */
struct pos_matrix {
        int** col_to_res;   /* [seq][col] -> residue pos or -1 */
        int numseq;
        int alnlen;
};

/* POAR entry: packed (pos_i, pos_j) key + support bitmask */
struct poar_entry {
        uint64_t key;       /* POAR_KEY(pos_i, pos_j) = (pos_i << 32) | pos_j */
        uint32_t support;   /* bitmask: bit k set if alignment k has this POAR */
};

/* Residue positions of the lower (i) and upper (j) sequence of a pair.  The
   key orders entries by (pos_i, pos_j).  (Kalign up to 4.0.0 packed
   pos_i << 20 | pos_j into 32 bits, so residues i and i + 4096 of a sequence
   longer than 4096 shared a key and their support was merged.) */
static inline uint64_t POAR_KEY(int pos_i, int pos_j)
{
        return ((uint64_t)(uint32_t)pos_i << 32) | (uint64_t)(uint32_t)pos_j;
}
static inline int POAR_KEY_I(uint64_t key) { return (int)(key >> 32); }
static inline int POAR_KEY_J(uint64_t key) { return (int)(key & 0xFFFFFFFFu); }

/* Stored sorted array of entries for one sequence pair (legacy backend:
   tables read from a POAR file). */
struct poar_pair {
        struct poar_entry* entries;
        int n_entries;
        int alloc_entries;
};

#define POAR_MAX_RUNS 32

/* POAR table.
 *
 * Implicit backend (ensemble): for every run k the table keeps only the
 * column of each residue, col[k][off[s] + r] = alignment column of residue r
 * of sequence s in run k.  The entries of a pair (i,j), i < j, are derived on
 * demand by poar_pair_entries(); they are exactly the sorted, bit-OR'd
 * POAR_KEY(pos_i, pos_j) keys a stored per-pair table would hold.  Memory is
 * O(n_runs * total residues).
 *
 * Stored backend (kalign_consensus_from_poar): pairs[] holds the entries read from a file.
 * pairs != NULL selects it. */
struct poar_table {
        int32_t* col[POAR_MAX_RUNS];  /* per run: residue -> column (NULL = run absent) */
        int run_ids[POAR_MAX_RUNS];   /* indices k with col[k] != NULL, ascending */
        int alnlen[POAR_MAX_RUNS];    /* per run: alignment length (columns) */
        int n_present;
        size_t* off;                  /* numseq + 1 residue offsets (NULL until first run) */
        int64_t* row_hist;            /* implicit backend, optional: [numseq][POAR_MAX_RUNS + 1], the number
                                         of entries of the pairs (i, j > i) of row i that s runs support
                                         (filled by score_runs_poar, dropped when a run is added) */
        struct poar_pair** pairs;     /* stored backend, indexed by pair_index */
        int64_t n_pairs;              /* numseq*(numseq-1)/2 (stored backend) */
        int numseq;
        int n_alignments;
};

/* Reusable buffer for poar_pair_entries; one per thread. */
struct poar_scratch {
        struct poar_entry* buf;
        size_t alloc;
};

/* Column -> residue map of one sequence j for every run (a thread-local
 * helper for the pair enumeration): after poar_dense_load(j), the pairs
 * (i,j) for any i < j are enumerated with plain array lookups. */
struct poar_dense {
        int32_t* d[POAR_MAX_RUNS];    /* per present run: column -> residue of the loaded sequence, -1 = gap */
        int j;                        /* loaded sequence or -1 */
};

/* Empty implicit table; runs are added with extract_poars. */
EXTERN int poar_table_alloc(struct poar_table** table, int numseq);
EXTERN void poar_table_free(struct poar_table* table);

EXTERN void poar_scratch_init(struct poar_scratch* sc);
EXTERN void poar_scratch_free(struct poar_scratch* sc);

/* Entries of pair (i,j), i < j: sorted by key, support = OR of (1 << k) over
 * the runs holding that key.  *entries points into the table (stored backend)
 * or into sc->buf (valid until the next call with the same scratch). */
EXTERN int poar_pair_entries(const struct poar_table* table, int i, int j,
                             struct poar_scratch* sc,
                             const struct poar_entry** entries, size_t* n);

EXTERN int poar_dense_init(const struct poar_table* table, struct poar_dense* dn);
EXTERN void poar_dense_free(const struct poar_table* table, struct poar_dense* dn);
/* Make sequence j the loaded one (releases the previous one). */
EXTERN void poar_dense_load(const struct poar_table* table, struct poar_dense* dn, int j);
/* Same result as poar_pair_entries(table, i, j, ...); faster when sequence j
 * is the one loaded in dn. */
EXTERN int poar_pair_entries_dense(const struct poar_table* table, struct poar_dense* dn,
                                   int i, int j, struct poar_scratch* sc,
                                   const struct poar_entry** entries, size_t* n);

/* The entries of poar_pair_entries_dense that lo..hi runs support (the
 * support has lo..hi bits set).  When sequence j is not the loaded one, or
 * the table is a stored one, all entries are returned, so the caller still
 * has to check the support. */
EXTERN int poar_pair_entries_dense_levels(const struct poar_table* table, struct poar_dense* dn,
                                          int i, int j, struct poar_scratch* sc, int lo, int hi,
                                          const struct poar_entry** entries, size_t* n);

/* Support bitmask of the key POAR_KEY(ri of i, rj of j) in the pair table: bit k
 * set if run k aligns the two residues.  i and j may be given in either order
 * (the lower-indexed sequence supplies the high key half). */
EXTERN uint32_t poar_support(const struct poar_table* table,
                             int i, int ri, int j, int rj);

/* 1 if the support of a residue pair of (i,j) is a plain column-equality
 * test on the run maps (implicit backend), 0 for a stored table, whose
 * entries are looked up through poar_pair_entries. */
EXTERN int poar_pair_is_exact(const struct poar_table* table, int i, int j);

EXTERN int pos_matrix_from_msa(struct pos_matrix** pm, char** seqs, int numseq, int alnlen);
EXTERN void pos_matrix_free(struct pos_matrix* pm);

/* Record the residue -> column map of alignment aln_idx in the table. */
EXTERN int extract_poars(struct poar_table* table, struct pos_matrix* pm,
                         int aln_idx);

EXTERN int poar_table_write(struct poar_table* table, const char* path);
EXTERN int poar_table_read(struct poar_table** table, const char* path);

#undef POAR_IMPORT
#undef EXTERN

#endif

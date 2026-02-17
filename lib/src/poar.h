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
        uint32_t key;       /* (pos_i << 20) | pos_j */
        uint32_t support;   /* bitmask: bit k set if alignment k has this POAR */
};

/* One sorted array per sequence pair */
struct poar_pair {
        struct poar_entry* entries;
        int n_entries;
        int alloc_entries;
};

/* Table of all POAR pairs */
struct poar_table {
        struct poar_pair** pairs;  /* indexed by flattened (i,j) pair */
        int n_pairs;               /* numseq*(numseq-1)/2 */
        int numseq;
        int n_alignments;
};

EXTERN int poar_table_alloc(struct poar_table** table, int numseq);
EXTERN void poar_table_free(struct poar_table* table);

EXTERN int pos_matrix_from_msa(struct pos_matrix** pm, char** seqs, int numseq, int alnlen);
EXTERN void pos_matrix_free(struct pos_matrix* pm);

EXTERN int extract_poars(struct poar_table* table, struct pos_matrix* pm, int aln_idx);

EXTERN int poar_table_write(struct poar_table* table, const char* path);
EXTERN int poar_table_read(struct poar_table** table, const char* path);

#undef POAR_IMPORT
#undef EXTERN

#endif

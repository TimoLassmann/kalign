#include "tldevel.h"

#include <ctype.h>
#include <string.h>

#define POAR_IMPORT
#include "poar.h"

static inline int pair_index(int i, int j, int numseq)
{
        /* Flat index for pair (i,j) where i < j */
        return i * numseq - (i * (i + 1)) / 2 + (j - i - 1);
}

static inline uint32_t pack_key(int pos_i, int pos_j)
{
        return ((uint32_t)pos_i << 20) | (uint32_t)pos_j;
}

static int poar_pair_alloc(struct poar_pair** pp)
{
        struct poar_pair* p = NULL;
        MMALLOC(p, sizeof(struct poar_pair));
        p->entries = NULL;
        p->n_entries = 0;
        p->alloc_entries = 64;
        MMALLOC(p->entries, sizeof(struct poar_entry) * p->alloc_entries);
        *pp = p;
        return OK;
ERROR:
        return FAIL;
}

static void poar_pair_free(struct poar_pair* p)
{
        if(p){
                if(p->entries){
                        MFREE(p->entries);
                }
                MFREE(p);
        }
}

/* Insert key into sorted array, or update support if already present. */
static int poar_pair_insert(struct poar_pair* p, uint32_t key, int aln_idx)
{
        /* Find insertion point */
        int lo = 0;
        int hi = p->n_entries;
        while(lo < hi){
                int mid = lo + (hi - lo) / 2;
                if(p->entries[mid].key < key){
                        lo = mid + 1;
                }else if(p->entries[mid].key == key){
                        /* Already exists, set bit */
                        p->entries[mid].support |= (1u << aln_idx);
                        return OK;
                }else{
                        hi = mid;
                }
        }

        /* Need to insert at position lo */
        if(p->n_entries >= p->alloc_entries){
                p->alloc_entries *= 2;
                MREALLOC(p->entries, sizeof(struct poar_entry) * p->alloc_entries);
        }

        /* Shift entries right */
        if(lo < p->n_entries){
                memmove(&p->entries[lo + 1], &p->entries[lo],
                        sizeof(struct poar_entry) * (p->n_entries - lo));
        }

        p->entries[lo].key = key;
        p->entries[lo].support = (1u << aln_idx);
        p->n_entries++;
        return OK;
ERROR:
        return FAIL;
}

int poar_table_alloc(struct poar_table** table, int numseq)
{
        struct poar_table* t = NULL;
        int n_pairs;
        int i;

        MMALLOC(t, sizeof(struct poar_table));
        t->pairs = NULL;
        t->numseq = numseq;
        t->n_alignments = 0;
        n_pairs = numseq * (numseq - 1) / 2;
        t->n_pairs = n_pairs;

        MMALLOC(t->pairs, sizeof(struct poar_pair*) * n_pairs);
        for(i = 0; i < n_pairs; i++){
                t->pairs[i] = NULL;
                RUN(poar_pair_alloc(&t->pairs[i]));
        }

        *table = t;
        return OK;
ERROR:
        poar_table_free(t);
        return FAIL;
}

void poar_table_free(struct poar_table* table)
{
        if(table){
                if(table->pairs){
                        for(int i = 0; i < table->n_pairs; i++){
                                poar_pair_free(table->pairs[i]);
                        }
                        MFREE(table->pairs);
                }
                MFREE(table);
        }
}

int pos_matrix_from_msa(struct pos_matrix** pm, char** seqs, int numseq, int alnlen)
{
        struct pos_matrix* m = NULL;
        int i, j;

        MMALLOC(m, sizeof(struct pos_matrix));
        m->col_to_res = NULL;
        m->numseq = numseq;
        m->alnlen = alnlen;

        MMALLOC(m->col_to_res, sizeof(int*) * numseq);
        for(i = 0; i < numseq; i++){
                m->col_to_res[i] = NULL;
                MMALLOC(m->col_to_res[i], sizeof(int) * alnlen);

                int res_pos = -1;
                for(j = 0; j < alnlen; j++){
                        if(isalpha((int)seqs[i][j])){
                                res_pos++;
                                m->col_to_res[i][j] = res_pos;
                        }else{
                                m->col_to_res[i][j] = -1;
                        }
                }
        }

        *pm = m;
        return OK;
ERROR:
        pos_matrix_free(m);
        return FAIL;
}

void pos_matrix_free(struct pos_matrix* pm)
{
        if(pm){
                if(pm->col_to_res){
                        for(int i = 0; i < pm->numseq; i++){
                                if(pm->col_to_res[i]){
                                        MFREE(pm->col_to_res[i]);
                                }
                        }
                        MFREE(pm->col_to_res);
                }
                MFREE(pm);
        }
}

int extract_poars(struct poar_table* table, struct pos_matrix* pm, int aln_idx)
{
        int i, j, col;
        int numseq = pm->numseq;
        int alnlen = pm->alnlen;

        ASSERT(aln_idx < 32, "Maximum 32 alignments supported in ensemble");

        for(i = 0; i < numseq - 1; i++){
                for(j = i + 1; j < numseq; j++){
                        int pidx = pair_index(i, j, numseq);
                        struct poar_pair* pp = table->pairs[pidx];

                        for(col = 0; col < alnlen; col++){
                                int ri = pm->col_to_res[i][col];
                                int rj = pm->col_to_res[j][col];
                                if(ri >= 0 && rj >= 0){
                                        uint32_t key = pack_key(ri, rj);
                                        RUN(poar_pair_insert(pp, key, aln_idx));
                                }
                        }
                }
        }

        if(aln_idx >= table->n_alignments){
                table->n_alignments = aln_idx + 1;
        }
        return OK;
ERROR:
        return FAIL;
}

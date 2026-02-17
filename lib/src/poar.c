#include "tldevel.h"

#include <ctype.h>
#include <string.h>
#include <stdio.h>

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

/* Binary POAR file format:
   4 bytes: "POAR" magic
   4 bytes: version (1)
   4 bytes: numseq
   4 bytes: n_alignments
   For each pair (n_pairs = numseq*(numseq-1)/2):
     4 bytes: n_entries
     n_entries * 8 bytes: (key, support) pairs
*/

#define POAR_MAGIC 0x524F4150  /* "POAR" in little-endian */
#define POAR_VERSION 1

int poar_table_write(struct poar_table* table, const char* path)
{
        FILE* fp = NULL;
        uint32_t magic = POAR_MAGIC;
        uint32_t version = POAR_VERSION;
        uint32_t numseq = (uint32_t)table->numseq;
        uint32_t n_alignments = (uint32_t)table->n_alignments;
        int i;

        ASSERT(table != NULL, "No POAR table");
        ASSERT(path != NULL, "No output path");

        fp = fopen(path, "wb");
        if(!fp){
                ERROR_MSG("Cannot open %s for writing", path);
        }

        fwrite(&magic, 4, 1, fp);
        fwrite(&version, 4, 1, fp);
        fwrite(&numseq, 4, 1, fp);
        fwrite(&n_alignments, 4, 1, fp);

        for(i = 0; i < table->n_pairs; i++){
                struct poar_pair* pp = table->pairs[i];
                uint32_t n_entries = (uint32_t)pp->n_entries;
                fwrite(&n_entries, 4, 1, fp);
                if(n_entries > 0){
                        fwrite(pp->entries, sizeof(struct poar_entry), n_entries, fp);
                }
        }

        fclose(fp);
        return OK;
ERROR:
        if(fp) fclose(fp);
        return FAIL;
}

int poar_table_read(struct poar_table** out_table, const char* path)
{
        FILE* fp = NULL;
        struct poar_table* t = NULL;
        uint32_t magic, version, numseq, n_alignments;
        int i;

        ASSERT(path != NULL, "No input path");

        fp = fopen(path, "rb");
        if(!fp){
                ERROR_MSG("Cannot open %s for reading", path);
        }

        if(fread(&magic, 4, 1, fp) != 1 || magic != POAR_MAGIC){
                ERROR_MSG("Invalid POAR file magic in %s", path);
        }
        if(fread(&version, 4, 1, fp) != 1 || version != POAR_VERSION){
                ERROR_MSG("Unsupported POAR file version %u in %s", version, path);
        }
        if(fread(&numseq, 4, 1, fp) != 1){
                ERROR_MSG("Failed to read numseq from %s", path);
        }
        if(fread(&n_alignments, 4, 1, fp) != 1){
                ERROR_MSG("Failed to read n_alignments from %s", path);
        }

        MMALLOC(t, sizeof(struct poar_table));
        t->pairs = NULL;
        t->numseq = (int)numseq;
        t->n_alignments = (int)n_alignments;
        t->n_pairs = (int)(numseq * (numseq - 1) / 2);

        MMALLOC(t->pairs, sizeof(struct poar_pair*) * t->n_pairs);
        for(i = 0; i < t->n_pairs; i++){
                t->pairs[i] = NULL;
        }

        for(i = 0; i < t->n_pairs; i++){
                uint32_t n_entries;
                struct poar_pair* pp = NULL;

                if(fread(&n_entries, 4, 1, fp) != 1){
                        ERROR_MSG("Failed to read pair %d entries count", i);
                }

                MMALLOC(pp, sizeof(struct poar_pair));
                pp->entries = NULL;
                pp->n_entries = (int)n_entries;
                pp->alloc_entries = n_entries > 0 ? (int)n_entries : 1;

                MMALLOC(pp->entries, sizeof(struct poar_entry) * pp->alloc_entries);

                if(n_entries > 0){
                        if(fread(pp->entries, sizeof(struct poar_entry), n_entries, fp) != n_entries){
                                MFREE(pp->entries);
                                MFREE(pp);
                                ERROR_MSG("Failed to read pair %d entries", i);
                        }
                }

                t->pairs[i] = pp;
        }

        fclose(fp);
        *out_table = t;
        return OK;
ERROR:
        if(fp) fclose(fp);
        if(t) poar_table_free(t);
        return FAIL;
}

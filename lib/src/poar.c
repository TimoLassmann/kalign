#include "tldevel.h"

#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define POAR_IMPORT
#include "poar.h"

static inline int64_t pair_index(int i, int j, int numseq)
{
        /* Flat index for pair (i,j) where i < j */
        return (int64_t)i * (int64_t)numseq - ((int64_t)i * ((int64_t)i + 1)) / 2
                + (int64_t)(j - i - 1);
}


int poar_table_alloc(struct poar_table** table, int numseq)
{
        struct poar_table* t = NULL;

        MMALLOC(t, sizeof(struct poar_table));
        for(int k = 0; k < POAR_MAX_RUNS; k++){
                t->col[k] = NULL;
                t->run_ids[k] = -1;
                t->alnlen[k] = 0;
        }
        t->n_present = 0;
        t->off = NULL;
        t->row_hist = NULL;
        t->pairs = NULL;
        t->n_pairs = 0;
        t->numseq = numseq;
        t->n_alignments = 0;

        *table = t;
        return OK;
ERROR:
        return FAIL;
}

void poar_table_free(struct poar_table* table)
{
        if(table){
                for(int k = 0; k < POAR_MAX_RUNS; k++){
                        if(table->col[k]){
                                MFREE(table->col[k]);
                        }
                }
                if(table->off){
                        MFREE(table->off);
                }
                if(table->row_hist){
                        MFREE(table->row_hist);
                }
                if(table->pairs){
                        for(int64_t i = 0; i < table->n_pairs; i++){
                                struct poar_pair* p = table->pairs[i];
                                if(p){
                                        if(p->entries){
                                                MFREE(p->entries);
                                        }
                                        MFREE(p);
                                }
                        }
                        MFREE(table->pairs);
                }
                MFREE(table);
        }
}

void poar_scratch_init(struct poar_scratch* sc)
{
        sc->buf = NULL;
        sc->alloc = 0;
}

void poar_scratch_free(struct poar_scratch* sc)
{
        if(sc->buf){
                free(sc->buf);
        }
        sc->buf = NULL;
        sc->alloc = 0;
}

static int poar_scratch_reserve(struct poar_scratch* sc, size_t n)
{
        if(n > sc->alloc){
                size_t na = sc->alloc ? sc->alloc : 256;
                struct poar_entry* nb = NULL;
                while(na < n){
                        na *= 2;
                }
                nb = realloc(sc->buf, sizeof(struct poar_entry) * na);
                if(nb == NULL){
                        ERROR_MSG("POAR scratch allocation of %zu entries failed", na);
                }
                sc->buf = nb;
                sc->alloc = na;
        }
        return OK;
ERROR:
        return FAIL;
}

int pos_matrix_from_msa(struct pos_matrix** pm, char** seqs, int numseq, int alnlen)
{ KPERF_REGION(KP_ENS_POAR);
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


/* Record the residue -> column map of alignment aln_idx.  Residues are the
   isalpha characters of the aligned sequences (as in pos_matrix_from_msa). */
int extract_poars(struct poar_table* table, struct pos_matrix* pm, int aln_idx)
{ KPERF_REGION(KP_ENS_POAR);
        int numseq = pm->numseq;
        int alnlen = pm->alnlen;
        int32_t* col = NULL;

        ASSERT(table != NULL, "No POAR table");
        ASSERT(table->pairs == NULL, "Cannot add runs to a stored POAR table");
        ASSERT(aln_idx >= 0 && aln_idx < POAR_MAX_RUNS, "Maximum 32 alignments supported in ensemble");
        ASSERT(table->col[aln_idx] == NULL, "POAR run %d already extracted", aln_idx);
        ASSERT(numseq == table->numseq, "POAR run has %d sequences, table %d", numseq, table->numseq);

        if(table->off == NULL){
                /* first run fixes the residue count of every sequence */
                MMALLOC(table->off, sizeof(size_t) * ((size_t)numseq + 1));
                table->off[0] = 0;
                for(int i = 0; i < numseq; i++){
                        size_t n = 0;
                        for(int c = 0; c < alnlen; c++){
                                if(pm->col_to_res[i][c] >= 0){
                                        n++;
                                }
                        }
                        table->off[i + 1] = table->off[i] + n;
                }
        }

        {
                size_t total = table->off[numseq];
                MMALLOC(col, sizeof(int32_t) * (total > 0 ? total : 1));
                for(int i = 0; i < numseq; i++){
                        size_t n = table->off[i + 1] - table->off[i];
                        size_t seen = 0;
                        int32_t* dst = col + table->off[i];
                        for(int c = 0; c < alnlen; c++){
                                int r = pm->col_to_res[i][c];
                                if(r >= 0){
                                        if((size_t)r >= n){
                                                ERROR_MSG("POAR run %d: sequence %d has more residues than in the first run", aln_idx, i);
                                        }
                                        dst[r] = c;
                                        seen++;
                                }
                        }
                        if(seen != n){
                                ERROR_MSG("POAR run %d: sequence %d has %zu residues, expected %zu",
                                          aln_idx, i, seen, n);
                        }
                }
        }

        table->col[aln_idx] = col;
        table->alnlen[aln_idx] = alnlen;
        col = NULL;
        if(table->row_hist){
                MFREE(table->row_hist);     /* counted without this run */
                table->row_hist = NULL;
        }
        table->n_present = 0;
        for(int k = 0; k < POAR_MAX_RUNS; k++){
                if(table->col[k]){
                        table->run_ids[table->n_present++] = k;
                }
        }
        if(aln_idx >= table->n_alignments){
                table->n_alignments = aln_idx + 1;
        }
        return OK;
ERROR:
        if(col) MFREE(col);
        return FAIL;
}

int poar_pair_is_exact(const struct poar_table* table, int i, int j)
{
        (void)i;
        (void)j;
        return table->pairs == NULL && table->off != NULL;
}

/* Next aligned residue pair of two strictly increasing residue -> column
   maps, from cursors (*pa, *pb).  Returns 1 and leaves the cursors on the
   match, or 0 when exhausted. */
static inline int next_match(const int32_t* a, size_t la,
                             const int32_t* b, size_t lb,
                             size_t* pa, size_t* pb)
{
        size_t x = *pa;
        size_t y = *pb;
        while(x < la && y < lb){
                int32_t u = a[x];
                int32_t v = b[y];
                if(u == v){
                        *pa = x;
                        *pb = y;
                        return 1;
                }
                if(u < v){
                        x++;
                }else{
                        y++;
                }
        }
        *pa = la;
        *pb = lb;
        return 0;
}

static int cmp_entry_key(const void* a, const void* b)
{
        uint64_t x = ((const struct poar_entry*)a)->key;
        uint64_t y = ((const struct poar_entry*)b)->key;
        return (x > y) - (x < y);
}

int poar_pair_entries(const struct poar_table* table, int i, int j,
                      struct poar_scratch* sc,
                      const struct poar_entry** entries, size_t* n)
{
        size_t li, lj, cap, ne = 0;
        int nr = table->n_present;

        if(table->pairs){
                struct poar_pair* pp = table->pairs[pair_index(i, j, table->numseq)];
                *entries = pp ? pp->entries : NULL;
                *n = pp ? (size_t)pp->n_entries : 0;
                return OK;
        }

        li = table->off ? table->off[i + 1] - table->off[i] : 0;
        lj = table->off ? table->off[j + 1] - table->off[j] : 0;
        cap = (li < lj ? li : lj) * (size_t)nr;
        if(cap == 0){
                *entries = sc->buf;
                *n = 0;
                return OK;
        }
        RUN(poar_scratch_reserve(sc, cap));

        /* Each run's matches come out in ascending (ri, rj) order, which
           is ascending key order.  K-way union with OR'd run bits. */
        const int32_t* a[POAR_MAX_RUNS];
        const int32_t* b[POAR_MAX_RUNS];
        size_t pa[POAR_MAX_RUNS];
        size_t pb[POAR_MAX_RUNS];
        uint64_t key[POAR_MAX_RUNS];
        const uint64_t DONE = UINT64_MAX;

        for(int r = 0; r < nr; r++){
                int k = table->run_ids[r];
                a[r] = table->col[k] + table->off[i];
                b[r] = table->col[k] + table->off[j];
                pa[r] = 0;
                pb[r] = 0;
                key[r] = next_match(a[r], li, b[r], lj, &pa[r], &pb[r])
                        ? POAR_KEY((int)pa[r], (int)pb[r]) : DONE;
        }
        for(;;){
                uint64_t m = DONE;
                uint32_t bits = 0;
                for(int r = 0; r < nr; r++){
                        if(key[r] < m){
                                m = key[r];
                        }
                }
                if(m == DONE){
                        break;
                }
                for(int r = 0; r < nr; r++){
                        if(key[r] == m){
                                bits |= 1u << table->run_ids[r];
                                pa[r]++;
                                pb[r]++;
                                key[r] = next_match(a[r], li, b[r], lj, &pa[r], &pb[r])
                                        ? POAR_KEY((int)pa[r], (int)pb[r]) : DONE;
                        }
                }
                sc->buf[ne].key = m;
                sc->buf[ne].support = bits;
                ne++;
        }
        *entries = sc->buf;
        *n = ne;
        return OK;
ERROR:
        return FAIL;
}

int poar_dense_init(const struct poar_table* table, struct poar_dense* dn)
{
        dn->j = -1;
        for(int k = 0; k < POAR_MAX_RUNS; k++){
                dn->d[k] = NULL;
        }
        if(table->pairs){
                return OK;
        }
        for(int r = 0; r < table->n_present; r++){
                int k = table->run_ids[r];
                size_t n = (size_t)table->alnlen[k];
                dn->d[k] = malloc(sizeof(int32_t) * (n ? n : 1));
                if(dn->d[k] == NULL){
                        ERROR_MSG("POAR dense map allocation failed");
                }
                for(size_t c = 0; c < n; c++){
                        dn->d[k][c] = -1;
                }
        }
        return OK;
ERROR:
        poar_dense_free(table, dn);
        return FAIL;
}

void poar_dense_free(const struct poar_table* table, struct poar_dense* dn)
{
        (void)table;
        for(int k = 0; k < POAR_MAX_RUNS; k++){
                if(dn->d[k]){
                        free(dn->d[k]);
                }
                dn->d[k] = NULL;
        }
        dn->j = -1;
}

static void dense_clear(const struct poar_table* table, struct poar_dense* dn)
{
        int j = dn->j;
        if(j >= 0){
                size_t lj = table->off[j + 1] - table->off[j];
                for(int r = 0; r < table->n_present; r++){
                        int k = table->run_ids[r];
                        const int32_t* b = table->col[k] + table->off[j];
                        int32_t* d = dn->d[k];
                        for(size_t y = 0; y < lj; y++){
                                d[b[y]] = -1;
                        }
                }
        }
        dn->j = -1;
}

void poar_dense_load(const struct poar_table* table, struct poar_dense* dn, int j)
{
        size_t lj;
        if(table->pairs || table->off == NULL){
                return;
        }
        dense_clear(table, dn);
        lj = table->off[j + 1] - table->off[j];
        for(int r = 0; r < table->n_present; r++){
                int k = table->run_ids[r];
                const int32_t* b = table->col[k] + table->off[j];
                int32_t* d = dn->d[k];
                for(size_t y = 0; y < lj; y++){
                        d[b[y]] = (int32_t)y;
                }
        }
        dn->j = j;
}

/* Entries of the pair, generated residue by residue of the lower sequence:
   for each residue x the partners y_r in the loaded sequence (one per run)
   are looked up in the column -> residue maps; distinct partners are emitted
   in ascending order with the OR of their run bits.  Same entries, in the
   same order, as the k-way union of the runs' match lists, restricted to
   the entries supported by lo..hi runs.  A residue whose partner is the
   same in every run gives one entry with all NR runs; any other residue
   gives entries with fewer than NR runs, so a pass that only wants the
   full level (lo = NR) skips those residues and a pass below it (hi < NR)
   skips the others. */
static inline __attribute__((always_inline))
size_t dense_enum(const struct poar_table* table, const struct poar_dense* dn,
                  int i, size_t li, struct poar_entry* out, const int NR,
                  const int lo, const int hi)
{
        const int32_t* a[POAR_MAX_RUNS];
        const int32_t* d[POAR_MAX_RUNS];
        uint32_t bit[POAR_MAX_RUNS];
        uint32_t full = 0;
        size_t ne = 0;
        const int want_full = NR >= lo && NR <= hi;
        const int want_part = lo < NR;
        int ident = 1;

        for(int r = 0; r < NR; r++){
                int k = table->run_ids[r];
                a[r] = table->col[k] + table->off[i];
                d[r] = dn->d[k];
                bit[r] = 1u << k;
                full |= bit[r];
                ident &= (k == r);      /* run r is alignment r: the run mask is the support */
        }
        for(size_t x = 0; x < li; x++){
                int32_t y[POAR_MAX_RUNS];
                int32_t diff = 0;
                for(int r = 0; r < NR; r++){
                        y[r] = d[r][a[r][x]];
                }
                for(int r = 1; r < NR; r++){
                        diff |= y[r] ^ y[0];
                }
                if(diff == 0){
                        if(want_full && y[0] >= 0){
                                out[ne].key = POAR_KEY((int)x, y[0]);
                                out[ne].support = full;
                                ne++;
                        }
                }else if(want_part){
                        /* eq[r]: the runs with the same partner as run r.  A partner is emitted once, by
                           its first run, if it is a residue and lo..hi runs share it; the selected
                           partners go out in ascending order. */
                        uint32_t eq[POAR_MAX_RUNS];
                        uint32_t sel = 0;
                        for(int r = 0; r < NR; r++){
                                eq[r] = 1u << r;
                        }
                        for(int r = 1; r < NR; r++){
                                for(int q = 0; q < r; q++){
                                        const uint32_t e = (uint32_t)(y[r] == y[q]);
                                        eq[r] |= e << q;
                                        eq[q] |= e << r;
                                }
                        }
                        for(int r = 0; r < NR; r++){
                                const int runs = __builtin_popcount(eq[r]);
                                const uint32_t ok = (uint32_t)(y[r] >= 0) & (uint32_t)((eq[r] & ((1u << r) - 1u)) == 0)
                                        & (uint32_t)(runs >= lo) & (uint32_t)(runs <= hi);
                                sel |= ok << r;
                        }
                        while(sel){
                                int best = __builtin_ctz(sel);
                                uint32_t rest = sel & (sel - 1u);
                                uint32_t bits;
                                while(rest){
                                        const int r = __builtin_ctz(rest);
                                        rest &= rest - 1u;
                                        if(y[r] < y[best]){
                                                best = r;
                                        }
                                }
                                sel &= ~(1u << best);
                                if(ident){
                                        bits = eq[best];
                                }else{
                                        uint32_t e = eq[best];
                                        bits = 0;
                                        while(e){
                                                bits |= bit[__builtin_ctz(e)];
                                                e &= e - 1u;
                                        }
                                }
                                out[ne].key = POAR_KEY((int)x, y[best]);
                                out[ne].support = bits;
                                ne++;
                        }
                }
        }
        return ne;
}

int poar_pair_entries_dense_levels(const struct poar_table* table, struct poar_dense* dn,
                                   int i, int j, struct poar_scratch* sc, int lo, int hi,
                                   const struct poar_entry** entries, size_t* n)
{
        size_t li, lj, cap;
        if(dn->j != j || table->pairs || !poar_pair_is_exact(table, i, j)){
                return poar_pair_entries(table, i, j, sc, entries, n);
        }
        li = table->off[i + 1] - table->off[i];
        lj = table->off[j + 1] - table->off[j];
        cap = (li < lj ? li : lj) * (size_t)table->n_present;
        if(cap == 0){
                *entries = sc->buf;
                *n = 0;
                return OK;
        }
        RUN(poar_scratch_reserve(sc, cap));
        switch(table->n_present){
        case 1: *n = dense_enum(table, dn, i, li, sc->buf, 1, lo, hi); break;
        case 2: *n = dense_enum(table, dn, i, li, sc->buf, 2, lo, hi); break;
        case 3: *n = dense_enum(table, dn, i, li, sc->buf, 3, lo, hi); break;
        case 4: *n = dense_enum(table, dn, i, li, sc->buf, 4, lo, hi); break;
        case 5: *n = dense_enum(table, dn, i, li, sc->buf, 5, lo, hi); break;
        default: *n = dense_enum(table, dn, i, li, sc->buf, table->n_present, lo, hi); break;
        }
        *entries = sc->buf;
        return OK;
ERROR:
        return FAIL;
}

int poar_pair_entries_dense(const struct poar_table* table, struct poar_dense* dn,
                            int i, int j, struct poar_scratch* sc,
                            const struct poar_entry** entries, size_t* n)
{
        return poar_pair_entries_dense_levels(table, dn, i, j, sc, 1, POAR_MAX_RUNS, entries, n);
}

uint32_t poar_support(const struct poar_table* table,
                      int i, int ri, int j, int rj)
{
        uint32_t bits = 0;
        size_t li, lj;
        int nr = table->n_present;

        if(i > j){
                int t = i; i = j; j = t;
                t = ri; ri = rj; rj = t;
        }

        if(table->pairs){
                const struct poar_pair* pp = table->pairs[pair_index(i, j, table->numseq)];
                uint64_t key = POAR_KEY(ri, rj);
                int lo = 0;
                int hi = pp ? pp->n_entries : 0;
                while(lo < hi){
                        int mid = lo + (hi - lo) / 2;
                        if(pp->entries[mid].key < key){
                                lo = mid + 1;
                        }else if(pp->entries[mid].key == key){
                                return pp->entries[mid].support;
                        }else{
                                hi = mid;
                        }
                }
                return 0;
        }

        if(table->off == NULL || ri < 0 || rj < 0){
                return 0;
        }
        li = table->off[i + 1] - table->off[i];
        lj = table->off[j + 1] - table->off[j];
        if((size_t)ri >= li || (size_t)rj >= lj){
                return 0;
        }

        for(int r = 0; r < nr; r++){
                int k = table->run_ids[r];
                const int32_t* c = table->col[k];
                if(c[table->off[i] + (size_t)ri] == c[table->off[j] + (size_t)rj]){
                        bits |= 1u << k;
                }
        }
        return bits;
}

/* Binary POAR file format:
   4 bytes: "POAR" magic
   4 bytes: version (2)
   4 bytes: numseq
   4 bytes: n_alignments
   For each pair (n_pairs = numseq*(numseq-1)/2):
     4 bytes: n_entries
     n_entries * 12 bytes: 8-byte key POAR_KEY(pos_i, pos_j), 4-byte support
   Version 1 (Kalign up to 4.0.0) stored 4-byte keys pos_i << 20 | pos_j;
   they are read and converted.  In a version 1 file, residue pairs whose
   lower position is 4096 or more were already merged when it was written. */

#define POAR_MAGIC 0x524F4150  /* "POAR" in little-endian */
#define POAR_VERSION 2

int poar_table_write(struct poar_table* table, const char* path)
{
        FILE* fp = NULL;
        uint32_t magic = POAR_MAGIC;
        uint32_t version = POAR_VERSION;
        uint32_t numseq;
        uint32_t n_alignments;
        struct poar_scratch sc;

        poar_scratch_init(&sc);
        ASSERT(table != NULL, "No POAR table");
        ASSERT(path != NULL, "No output path");
        numseq = (uint32_t)table->numseq;
        n_alignments = (uint32_t)table->n_alignments;

        fp = fopen(path, "wb");
        if(!fp){
                ERROR_MSG("Cannot open %s for writing", path);
        }

        fwrite(&magic, 4, 1, fp);
        fwrite(&version, 4, 1, fp);
        fwrite(&numseq, 4, 1, fp);
        fwrite(&n_alignments, 4, 1, fp);

        for(int i = 0; i < table->numseq - 1; i++){
                for(int j = i + 1; j < table->numseq; j++){
                        const struct poar_entry* e = NULL;
                        size_t n = 0;
                        uint32_t n_entries;
                        RUN(poar_pair_entries(table, i, j, &sc, &e, &n));
                        n_entries = (uint32_t)n;
                        fwrite(&n_entries, 4, 1, fp);
                        for(uint32_t x = 0; x < n_entries; x++){
                                fwrite(&e[x].key, 8, 1, fp);
                                fwrite(&e[x].support, 4, 1, fp);
                        }
                }
        }

        poar_scratch_free(&sc);
        fclose(fp);
        return OK;
ERROR:
        poar_scratch_free(&sc);
        if(fp) fclose(fp);
        return FAIL;
}

int poar_table_read(struct poar_table** out_table, const char* path)
{
        FILE* fp = NULL;
        struct poar_table* t = NULL;
        uint32_t magic, version, numseq, n_alignments;
        int64_t i;

        ASSERT(path != NULL, "No input path");

        fp = fopen(path, "rb");
        if(!fp){
                ERROR_MSG("Cannot open %s for reading", path);
        }

        if(fread(&magic, 4, 1, fp) != 1 || magic != POAR_MAGIC){
                ERROR_MSG("Invalid POAR file magic in %s", path);
        }
        if(fread(&version, 4, 1, fp) != 1 || (version != 1 && version != POAR_VERSION)){
                ERROR_MSG("Unsupported POAR file version %u in %s", version, path);
        }
        if(fread(&numseq, 4, 1, fp) != 1){
                ERROR_MSG("Failed to read numseq from %s", path);
        }
        if(fread(&n_alignments, 4, 1, fp) != 1){
                ERROR_MSG("Failed to read n_alignments from %s", path);
        }
        if(numseq > (uint32_t)INT_MAX){
                ERROR_MSG("POAR file %s: numseq=%u exceeds INT_MAX", path, numseq);
        }

        RUN(poar_table_alloc(&t, (int)numseq));
        t->n_alignments = (int)n_alignments;
        t->n_pairs = numseq < 2 ? 0 : (int64_t)((uint64_t)numseq * ((uint64_t)numseq - 1) / 2);

        MMALLOC(t->pairs, sizeof(struct poar_pair*) * (size_t)(t->n_pairs > 0 ? t->n_pairs : 1));
        for(i = 0; i < t->n_pairs; i++){
                t->pairs[i] = NULL;
        }

        for(i = 0; i < t->n_pairs; i++){
                uint32_t n_entries;
                struct poar_pair* pp = NULL;

                if(fread(&n_entries, 4, 1, fp) != 1){
                        ERROR_MSG("Failed to read pair %lld entries count", (long long)i);
                }

                /* Cap so the cast to int doesn't wrap negative and the per-entry
                   allocation size doesn't overflow size_t. */
                if(n_entries > (uint32_t)(INT_MAX / sizeof(struct poar_entry))){
                        ERROR_MSG("POAR file: pair %lld has n_entries=%u, exceeds bounds",
                                  (long long)i, n_entries);
                }

                MMALLOC(pp, sizeof(struct poar_pair));
                pp->entries = NULL;
                pp->n_entries = (int)n_entries;
                pp->alloc_entries = n_entries > 0 ? (int)n_entries : 1;

                MMALLOC(pp->entries, sizeof(struct poar_entry) * pp->alloc_entries);

                for(uint32_t x = 0; x < n_entries; x++){
                        struct poar_entry* pe = &pp->entries[x];
                        int ok;
                        if(version == 1){
                                uint32_t k32;
                                ok = fread(&k32, 4, 1, fp) == 1;
                                pe->key = POAR_KEY((int)(k32 >> 20), (int)(k32 & 0xFFFFFu));
                        }else{
                                ok = fread(&pe->key, 8, 1, fp) == 1;
                        }
                        ok = ok && fread(&pe->support, 4, 1, fp) == 1;
                        if(!ok){
                                MFREE(pp->entries);
                                MFREE(pp);
                                ERROR_MSG("Failed to read pair %lld entries", (long long)i);
                        }
                }
                if(version == 1){
                        qsort(pp->entries, n_entries, sizeof(struct poar_entry), cmp_entry_key);
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

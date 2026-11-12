#include "tldevel.h"

#include <string.h>

#include "msa_struct.h"
#include "msa_alloc.h"
#include "poar.h"
#include "kalign_bits.h"

#include "threadpool/threadpool.h"

#define CONSENSUS_MSA_IMPORT
#include "consensus_msa.h"

/* Union-find data structure with per-set sequence membership tracking
   and element linked lists for efficient set enumeration.

   Column DAG edges.  A column (set) has an edge to the column of the next
   residue of each of its members.  Many members share one successor column,
   so each root keeps a list of the members whose successors are known to be
   in distinct columns: succ_head/succ_tail/succ_next (and pred_* for the
   previous residues).  The lists start with every member that has a
   successor (predecessor) and are thinned while they are scanned
   (uf_scan_*): a member whose neighbour column was already seen in the same
   scan is unlinked, for good, because two columns that are equal stay equal.
   A scan therefore costs the number of distinct neighbour columns, not the
   number of members, and visits exactly the same set of columns. */
struct uf_set {
        int* parent;
        int* rank;
        int* elem_seq;          /* element -> sequence index */
        uint64_t** seq_mask;    /* per-root: bitmask of sequences in set */
        int* set_head;          /* root -> first element in set (-1 if none) */
        int* set_tail;          /* root -> last element in set */
        int* next_in_set;       /* element -> next element in same set (-1 if last) */
        int* succ_head;         /* root -> first member with a successor edge (-1 if none) */
        int* succ_tail;
        int* succ_next;
        int* pred_head;         /* root -> first member with a predecessor edge (-1 if none) */
        int* pred_tail;
        int* pred_next;
        uint32_t* dup;          /* root -> stamp of the scan that last saw it as a neighbour */
        uint32_t dup_stamp;
        int n_elements;
        int numseq;
        int mask_words;         /* number of uint64_t words per bitmask */
};

static void uf_free(struct uf_set* u);

static int uf_alloc(struct uf_set** uf, int n, int* seq_offsets,
                    int* seq_lengths, int numseq)
{
        struct uf_set* u = NULL;
        int mw = (numseq + 63) / 64;

        MMALLOC(u, sizeof(struct uf_set));
        u->parent = NULL;
        u->rank = NULL;
        u->elem_seq = NULL;
        u->seq_mask = NULL;
        u->set_head = NULL;
        u->set_tail = NULL;
        u->next_in_set = NULL;
        u->succ_head = NULL;
        u->succ_tail = NULL;
        u->succ_next = NULL;
        u->pred_head = NULL;
        u->pred_tail = NULL;
        u->pred_next = NULL;
        u->dup = NULL;
        u->dup_stamp = 0;
        u->n_elements = n;
        u->numseq = numseq;
        u->mask_words = mw;

        MMALLOC(u->parent, sizeof(int) * n);
        MMALLOC(u->rank, sizeof(int) * n);
        MMALLOC(u->elem_seq, sizeof(int) * n);
        MMALLOC(u->set_head, sizeof(int) * n);
        MMALLOC(u->next_in_set, sizeof(int) * n);
        MMALLOC(u->seq_mask, sizeof(uint64_t*) * n);
        /* seq_mask first: uf_free walks it */
        for(int i = 0; i < n; i++){
                u->seq_mask[i] = NULL;
        }
        MMALLOC(u->set_tail, sizeof(int) * n);
        MMALLOC(u->succ_head, sizeof(int) * n);
        MMALLOC(u->succ_tail, sizeof(int) * n);
        MMALLOC(u->succ_next, sizeof(int) * n);
        MMALLOC(u->pred_head, sizeof(int) * n);
        MMALLOC(u->pred_tail, sizeof(int) * n);
        MMALLOC(u->pred_next, sizeof(int) * n);
        MMALLOC(u->dup, sizeof(uint32_t) * n);
        memset(u->dup, 0, sizeof(uint32_t) * n);

        for(int i = 0; i < n; i++){
                u->parent[i] = i;
                u->rank[i] = 0;
                u->set_head[i] = i;
                u->set_tail[i] = i;
                u->next_in_set[i] = -1;
                /* Bitmasks are materialised only for merged sets; a
                   singleton's mask is just the bit of its own sequence. */
                u->seq_mask[i] = NULL;
        }

        /* Fill elem_seq */
        for(int s = 0; s < numseq; s++){
                for(int p = 0; p < seq_lengths[s]; p++){
                        int elem = seq_offsets[s] + p;
                        u->elem_seq[elem] = s;
                        /* every residue starts as its own column: one
                           successor edge unless it is the last residue,
                           one predecessor edge unless it is the first */
                        u->succ_next[elem] = -1;
                        u->pred_next[elem] = -1;
                        u->succ_head[elem] = u->succ_tail[elem] = (p + 1 < seq_lengths[s]) ? elem : -1;
                        u->pred_head[elem] = u->pred_tail[elem] = (p > 0) ? elem : -1;
                }
        }

        *uf = u;
        return OK;
ERROR:
        uf_free(u);
        return FAIL;
}

static void uf_free(struct uf_set* u)
{
        if(u){
                if(u->seq_mask){
                        for(int i = 0; i < u->n_elements; i++){
                                if(u->seq_mask[i]) MFREE(u->seq_mask[i]);
                        }
                        MFREE(u->seq_mask);
                }
                if(u->next_in_set) MFREE(u->next_in_set);
                if(u->set_head) MFREE(u->set_head);
                if(u->set_tail) MFREE(u->set_tail);
                if(u->succ_head) MFREE(u->succ_head);
                if(u->succ_tail) MFREE(u->succ_tail);
                if(u->succ_next) MFREE(u->succ_next);
                if(u->pred_head) MFREE(u->pred_head);
                if(u->pred_tail) MFREE(u->pred_tail);
                if(u->pred_next) MFREE(u->pred_next);
                if(u->dup) MFREE(u->dup);
                if(u->elem_seq) MFREE(u->elem_seq);
                if(u->parent) MFREE(u->parent);
                if(u->rank) MFREE(u->rank);
                MFREE(u);
        }
}

static int uf_find(struct uf_set* u, int x)
{
        while(u->parent[x] != x){
                u->parent[x] = u->parent[u->parent[x]]; /* path halving */
                x = u->parent[x];
        }
        return x;
}

/* A fresh stamp for one neighbour scan of one column. */
static inline uint32_t uf_scan_stamp(struct uf_set* u)
{
        if(++u->dup_stamp == 0){
                memset(u->dup, 0, sizeof(uint32_t) * (size_t)u->n_elements);
                u->dup_stamp = 1;
        }
        return u->dup_stamp;
}

/* Online topological order over consensus columns (Pearce & Kelly, 2006).

   Columns are union-find roots. A directed edge col_A -> col_B exists
   when some sequence has consecutive residues in cols A and B. ord[]
   holds a topological order of this DAG and is updated incrementally as
   columns merge, so a cycle check only has to search the columns that
   lie between the two merge candidates in the current order. */
struct ord_node {
        int ord;
        int node;
};

struct topo_order {
        int* ord;                /* uf root -> position in topological order */
        int* mark;               /* uf root -> visit stamp */
        int stamp;
        int* stack;              /* DFS stack */
        struct ord_node* fwd;    /* roots reachable from lo inside the window */
        struct ord_node* bwd;    /* roots reaching hi inside the window */
        int* slots;              /* order positions freed by fwd + bwd */
        struct ord_node* tmp;    /* radix sort buffer */
        int key_bits;            /* bits needed for an order position */
        int n_fwd;
        int n_bwd;
};

static int topo_alloc(struct topo_order** topo, int n)
{
        struct topo_order* t = NULL;

        MMALLOC(t, sizeof(struct topo_order));
        t->ord = NULL;
        t->mark = NULL;
        t->stack = NULL;
        t->fwd = NULL;
        t->bwd = NULL;
        t->slots = NULL;
        t->tmp = NULL;
        t->key_bits = 1;
        while(t->key_bits < 31 && (1 << t->key_bits) < n){
                t->key_bits++;
        }
        t->stamp = 0;
        t->n_fwd = 0;
        t->n_bwd = 0;

        MMALLOC(t->ord, sizeof(int) * n);
        MMALLOC(t->mark, sizeof(int) * n);
        MMALLOC(t->stack, sizeof(int) * n);
        MMALLOC(t->fwd, sizeof(struct ord_node) * n);
        MMALLOC(t->bwd, sizeof(struct ord_node) * n);
        MMALLOC(t->slots, sizeof(int) * n);
        MMALLOC(t->tmp, sizeof(struct ord_node) * n);
        memset(t->mark, 0, sizeof(int) * n);

        *topo = t;
        return OK;
ERROR:
        return FAIL;
}

static void topo_free(struct topo_order* t)
{
        if(t){
                if(t->ord) MFREE(t->ord);
                if(t->mark) MFREE(t->mark);
                if(t->stack) MFREE(t->stack);
                if(t->fwd) MFREE(t->fwd);
                if(t->bwd) MFREE(t->bwd);
                if(t->slots) MFREE(t->slots);
                if(t->tmp) MFREE(t->tmp);
                MFREE(t);
        }
}

struct elem_key {
        double key;
        int elem;
};

static int cmp_elem_key(const void* a, const void* b)
{
        const struct elem_key* x = (const struct elem_key*)a;
        const struct elem_key* y = (const struct elem_key*)b;
        if(x->key < y->key) return -1;
        if(x->key > y->key) return 1;
        return x->elem - y->elem;
}

/* Initial order: every residue is its own column, so the DAG is a set of
   disjoint chains. Interleave the chains by relative sequence position so
   that likely merge partners start close together in the order, which
   keeps the Pearce-Kelly search windows small. */
static int topo_init(struct topo_order* t, int* seq_offsets,
                     int* seq_lengths, int numseq, int n)
{
        struct elem_key* keys = NULL;

        MMALLOC(keys, sizeof(struct elem_key) * n);
        for(int s = 0; s < numseq; s++){
                for(int p = 0; p < seq_lengths[s]; p++){
                        int e = seq_offsets[s] + p;
                        keys[e].key = ((double)p + 0.5) / (double)seq_lengths[s];
                        keys[e].elem = e;
                }
        }
        qsort(keys, n, sizeof(struct elem_key), cmp_elem_key);
        for(int i = 0; i < n; i++){
                t->ord[keys[i].elem] = i;
        }
        MFREE(keys);
        return OK;
ERROR:
        return FAIL;
}

/* Collect roots reachable from lo whose order is below ord[hi]. Only
   those can lie on a path lo -> hi. Returns 1 if hi is reachable. */
static int topo_search_fwd(struct uf_set* u, struct topo_order* t,
                           int lo, int hi)
{
        int ub = t->ord[hi];
        int sp = 0;

        t->stamp++;
        t->n_fwd = 0;
        t->mark[lo] = t->stamp;
        t->stack[sp++] = lo;

        while(sp > 0){
                int cur = t->stack[--sp];
                t->fwd[t->n_fwd].ord = t->ord[cur];
                t->fwd[t->n_fwd].node = cur;
                t->n_fwd++;

                /* distinct successor columns of cur; duplicates are unlinked */
                const uint32_t ds = uf_scan_stamp(u);
                int prev_elem = -1;
                int nx;
                for(int elem = u->succ_head[cur]; elem >= 0; elem = nx){
                        nx = u->succ_next[elem];
                        int next = uf_find(u, elem + 1);
                        if(u->dup[next] == ds){
                                if(prev_elem < 0){
                                        u->succ_head[cur] = nx;
                                }else{
                                        u->succ_next[prev_elem] = nx;
                                }
                                if(nx < 0){
                                        u->succ_tail[cur] = prev_elem;
                                }
                                continue;
                        }
                        u->dup[next] = ds;
                        prev_elem = elem;
                        if(next == hi){
                                return 1;
                        }
                        if(t->ord[next] < ub && t->mark[next] != t->stamp){
                                t->mark[next] = t->stamp;
                                t->stack[sp++] = next;
                        }
                }
        }
        return 0;
}

/* Collect roots that reach hi and whose order is above ord[lo]. */
static void topo_search_bwd(struct uf_set* u, struct topo_order* t,
                            int lo, int hi)
{
        int lb = t->ord[lo];
        int sp = 0;

        t->stamp++;
        t->n_bwd = 0;
        t->mark[hi] = t->stamp;
        t->stack[sp++] = hi;

        while(sp > 0){
                int cur = t->stack[--sp];
                t->bwd[t->n_bwd].ord = t->ord[cur];
                t->bwd[t->n_bwd].node = cur;
                t->n_bwd++;

                /* distinct predecessor columns of cur; duplicates are unlinked */
                const uint32_t ds = uf_scan_stamp(u);
                int prev_elem = -1;
                int nx;
                for(int elem = u->pred_head[cur]; elem >= 0; elem = nx){
                        nx = u->pred_next[elem];
                        int prev = uf_find(u, elem - 1);
                        if(u->dup[prev] == ds){
                                if(prev_elem < 0){
                                        u->pred_head[cur] = nx;
                                }else{
                                        u->pred_next[prev_elem] = nx;
                                }
                                if(nx < 0){
                                        u->pred_tail[cur] = prev_elem;
                                }
                                continue;
                        }
                        u->dup[prev] = ds;
                        prev_elem = elem;
                        if(t->ord[prev] > lb && t->mark[prev] != t->stamp){
                                t->mark[prev] = t->stamp;
                                t->stack[sp++] = prev;
                        }
                }
        }
}

/* Sort by ord, ascending.  The ord values of distinct columns are distinct,
   so every correct sort gives the same array; this one avoids the indirect
   comparison calls of qsort (it runs once or twice per merge). */
static void ord_sort(struct ord_node* a, int n)
{
        while(n > 24){
                const int x = a[0].ord;
                const int y = a[n >> 1].ord;
                const int z = a[n - 1].ord;
                const int piv = x < y ? (y < z ? y : (x < z ? z : x))
                                      : (x < z ? x : (y < z ? z : y));
                int i = 0;
                int j = n - 1;
                while(i <= j){
                        while(a[i].ord < piv) i++;
                        while(a[j].ord > piv) j--;
                        if(i <= j){
                                struct ord_node tmp = a[i];
                                a[i] = a[j];
                                a[j] = tmp;
                                i++;
                                j--;
                        }
                }
                /* recurse into the smaller part, loop on the larger */
                if(j + 1 < n - i){
                        ord_sort(a, j + 1);
                        a += i;
                        n -= i;
                }else{
                        ord_sort(a + i, n - i);
                        n = j + 1;
                }
        }
        for(int i = 1; i < n; i++){
                struct ord_node tmp = a[i];
                int j = i - 1;
                while(j >= 0 && a[j].ord > tmp.ord){
                        a[j + 1] = a[j];
                        j--;
                }
                a[j + 1] = tmp;
        }
}

/* ord_sort for large arrays: least-significant-digit radix sort on the
   order position (key_bits bits), 11 bits per pass. */
#define ORD_RADIX_BITS 11
#define ORD_RADIX_MIN 512
static void ord_sort_big(struct ord_node* a, int n, struct ord_node* tmp, int key_bits)
{
        struct ord_node* src = a;
        struct ord_node* dst = tmp;
        if(n < ORD_RADIX_MIN){
                ord_sort(a, n);
                return;
        }
        for(int shift = 0; shift < key_bits; shift += ORD_RADIX_BITS){
                int cnt[1 << ORD_RADIX_BITS];
                int sum = 0;
                memset(cnt, 0, sizeof(cnt));
                for(int i = 0; i < n; i++){
                        cnt[(src[i].ord >> shift) & ((1 << ORD_RADIX_BITS) - 1)]++;
                }
                for(int d = 0; d < (1 << ORD_RADIX_BITS); d++){
                        int c = cnt[d];
                        cnt[d] = sum;
                        sum += c;
                }
                for(int i = 0; i < n; i++){
                        dst[cnt[(src[i].ord >> shift) & ((1 << ORD_RADIX_BITS) - 1)]++] = src[i];
                }
                struct ord_node* sw = src;
                src = dst;
                dst = sw;
        }
        if(src != a){
                memcpy(a, src, sizeof(struct ord_node) * (size_t)n);
        }
}

/* Reassign the order positions held by bwd and fwd so that every bwd
   root precedes every fwd root, keeping relative order within each set.
   hi ends up in the last bwd slot, directly before lo's successors. */
static void topo_reorder(struct topo_order* t)
{
        int i = 0, j = 0, k = 0;

        ord_sort_big(t->bwd, t->n_bwd, t->tmp, t->key_bits);
        ord_sort_big(t->fwd, t->n_fwd, t->tmp, t->key_bits);

        while(i < t->n_bwd && j < t->n_fwd){
                if(t->bwd[i].ord < t->fwd[j].ord){
                        t->slots[k++] = t->bwd[i++].ord;
                }else{
                        t->slots[k++] = t->fwd[j++].ord;
                }
        }
        while(i < t->n_bwd) t->slots[k++] = t->bwd[i++].ord;
        while(j < t->n_fwd) t->slots[k++] = t->fwd[j++].ord;

        k = 0;
        for(i = 0; i < t->n_bwd; i++){
                t->ord[t->bwd[i].node] = t->slots[k++];
        }
        for(j = 0; j < t->n_fwd; j++){
                t->ord[t->fwd[j].node] = t->slots[k++];
        }
}

/* Merge two sets. Returns 1 if merge succeeded, 0 if blocked by:
   - same-sequence conflict (two residues from same seq in one column)
   - ordering cycle (merging would create a cycle in the column DAG) */
static int uf_union_safe(struct uf_set* u, struct topo_order* t, int a, int b)
{
        int ra = uf_find(u, a);
        int rb = uf_find(u, b);
        if(ra == rb) return 1;

        /* Check for conflict: do the two sets share any sequence? */
        if(u->seq_mask[ra] == NULL && u->seq_mask[rb] == NULL){
                if(u->elem_seq[ra] == u->elem_seq[rb]){
                        return 0;
                }
        }else if(u->seq_mask[ra] == NULL || u->seq_mask[rb] == NULL){
                int single = u->seq_mask[ra] == NULL ? ra : rb;
                const uint64_t* m = u->seq_mask[ra] == NULL ? u->seq_mask[rb] : u->seq_mask[ra];
                int sq = u->elem_seq[single];
                if(m[sq / 64] & (1ULL << (sq % 64))){
                        return 0;
                }
        }else{
                for(int w = 0; w < u->mask_words; w++){
                        if(u->seq_mask[ra][w] & u->seq_mask[rb][w]){
                                return 0; /* conflict: same sequence in both sets */
                        }
                }
        }

        /* In a topological order a path can only run from the earlier
           column to the later one, and merging the two creates a cycle
           exactly when such a path exists. */
        int lo = ra;
        int hi = rb;
        if(t->ord[lo] > t->ord[hi]){
                lo = rb;
                hi = ra;
        }
        if(topo_search_fwd(u, t, lo, hi)){
                return 0; /* would create cycle */
        }
        topo_search_bwd(u, t, lo, hi);
        topo_reorder(t);
        int merged_ord = t->ord[hi];

        /* Safe to merge */
        int new_root;
        int old_root;
        if(u->rank[ra] < u->rank[rb]){
                u->parent[ra] = rb;
                new_root = rb;
                old_root = ra;
        }else if(u->rank[ra] > u->rank[rb]){
                u->parent[rb] = ra;
                new_root = ra;
                old_root = rb;
        }else{
                u->parent[rb] = ra;
                u->rank[ra]++;
                new_root = ra;
                old_root = rb;
        }
        t->ord[new_root] = merged_ord;

        /* Merge sequence bitmasks into new root (materialising singleton
           masks on demand) and release the old root's mask */
        if(u->seq_mask[new_root] == NULL){
                int sq = u->elem_seq[new_root];
                MMALLOC(u->seq_mask[new_root], sizeof(uint64_t) * u->mask_words);
                memset(u->seq_mask[new_root], 0, sizeof(uint64_t) * u->mask_words);
                u->seq_mask[new_root][sq / 64] |= (1ULL << (sq % 64));
        }
        if(u->seq_mask[old_root] == NULL){
                int sq = u->elem_seq[old_root];
                u->seq_mask[new_root][sq / 64] |= (1ULL << (sq % 64));
        }else{
                for(int w = 0; w < u->mask_words; w++){
                        u->seq_mask[new_root][w] |= u->seq_mask[old_root][w];
                }
                MFREE(u->seq_mask[old_root]);
        }

        /* Concatenate the member lists and the edge lists: append
           old_root's to new_root's */
        if(u->set_head[old_root] >= 0){
                if(u->set_head[new_root] < 0){
                        u->set_head[new_root] = u->set_head[old_root];
                }else{
                        u->next_in_set[u->set_tail[new_root]] = u->set_head[old_root];
                }
                u->set_tail[new_root] = u->set_tail[old_root];
        }
        u->set_head[old_root] = -1;
        if(u->succ_head[old_root] >= 0){
                if(u->succ_head[new_root] < 0){
                        u->succ_head[new_root] = u->succ_head[old_root];
                }else{
                        u->succ_next[u->succ_tail[new_root]] = u->succ_head[old_root];
                }
                u->succ_tail[new_root] = u->succ_tail[old_root];
                u->succ_head[old_root] = -1;
        }
        if(u->pred_head[old_root] >= 0){
                if(u->pred_head[new_root] < 0){
                        u->pred_head[new_root] = u->pred_head[old_root];
                }else{
                        u->pred_next[u->pred_tail[new_root]] = u->pred_head[old_root];
                }
                u->pred_tail[new_root] = u->pred_tail[old_root];
                u->pred_head[old_root] = -1;
        }

        return 1;
ERROR:
        return -1;
}

static inline int popcount32(uint32_t x)
{
        return kalign_popcount32(x);
}

#define MAX_SUPPORT 32
#if MAX_SUPPORT != POAR_MAX_RUNS
#error "MAX_SUPPORT must equal POAR_MAX_RUNS (row_hist layout)"
#endif
/* Candidates are held in memory for at most this many at a time; a support
   level larger than that is processed in blocks of consecutive rows. */
#ifndef CAND_BUDGET
#define CAND_BUDGET ((int64_t)1 << 21)
#endif

struct cand_pair {
        int elem_i;
        int elem_j;
};

static void for_rows(threadpool_t* pool, int numseq, int r0, int r1,
                     void (*fn)(int, int, void*), void* arg)
{
        if(pool != NULL && numseq > 16){
                tp_parallel_for(pool, r0, r1, fn, arg);
        }else{
                fn(r0, r1, arg);
        }
}
#define FOR_ROWS(pool, numseq, r0, r1, fn, arg) for_rows(pool, numseq, r0, r1, fn, arg)

struct cand_ctx {
        const struct poar_table* table;
        int numseq;
        int min_level;              /* hist: lowest level counted */
        int lo;                     /* store: levels [lo, hi] have a bucket */
        int hi;
        const int* seq_offsets;
        const int* uf_parent;       /* store: union-find parents, not modified during the pass */
        int64_t* row_hist;          /* [numseq][MAX_SUPPORT + 1] */
        const int64_t* row_start;   /* store: [numseq][MAX_SUPPORT + 1] bucket offsets */
        struct cand_pair** bucket;  /* [MAX_SUPPORT + 1] */
        int error;
};

/* Count, per row and support level, the POAR entries of a pair.  Sequence j
   is loaded once and paired with every row i < j of the chunk. */
static void cand_hist_chunk(int start, int end, void* arg)
{
        struct cand_ctx* c = (struct cand_ctx*)arg;
        struct poar_scratch sc;
        struct poar_dense dn;
        poar_scratch_init(&sc);
        if(poar_dense_init(c->table, &dn) != OK){
                c->error = 1;
                poar_scratch_free(&sc);
                return;
        }
        for(int j = start + 1; j < c->numseq; j++){
                int i1 = j < end ? j : end;
                poar_dense_load(c->table, &dn, j);
                for(int i = start; i < i1; i++){
                        int64_t* h = c->row_hist + (size_t)i * (MAX_SUPPORT + 1);
                        const struct poar_entry* e = NULL;
                        size_t n = 0;
                        if(poar_pair_entries_dense_levels(c->table, &dn, i, j, &sc, c->min_level, MAX_SUPPORT,
                                                          &e, &n) != OK){
                                c->error = 1;
                                poar_dense_free(c->table, &dn);
                                poar_scratch_free(&sc);
                                return;
                        }
                        for(size_t x = 0; x < n; x++){
                                int pc = popcount32(e[x].support);
                                if(pc >= c->min_level){
                                        h[pc]++;
                                }
                        }
                }
        }
        poar_dense_free(c->table, &dn);
        poar_scratch_free(&sc);
}

/* Root of x without path compression: safe for concurrent readers while
   nothing merges. */
static inline int uf_find_ro(const int* parent, int x)
{
        while(parent[x] != x){
                x = parent[x];
        }
        return x;
}

/* Fill the level buckets: within a level, candidates land in (i asc, j asc,
   key asc) order.  A pair whose residues are already in one column when
   the pass runs can never change anything (columns only grow), so its
   slot is marked empty (elem_i = -1) and the serial merge loop skips it;
   about nine in ten candidates of a large input end this way. */
static void cand_store_chunk(int start, int end, void* arg)
{
        struct cand_ctx* c = (struct cand_ctx*)arg;
        struct poar_scratch sc;
        struct poar_dense dn;
        int64_t* pos = NULL;
        poar_scratch_init(&sc);
        if(end <= start){
                return;
        }
        pos = malloc(sizeof(int64_t) * (size_t)(end - start) * (MAX_SUPPORT + 1));
        if(pos == NULL || poar_dense_init(c->table, &dn) != OK){
                c->error = 1;
                free(pos);
                return;
        }
        for(int i = start; i < end; i++){
                for(int v = 0; v <= MAX_SUPPORT; v++){
                        pos[(size_t)(i - start) * (MAX_SUPPORT + 1) + v] =
                                c->row_start[(size_t)i * (MAX_SUPPORT + 1) + v];
                }
        }
        for(int j = start + 1; j < c->numseq; j++){
                int i1 = j < end ? j : end;
                poar_dense_load(c->table, &dn, j);
                for(int i = start; i < i1; i++){
                        int64_t* p = pos + (size_t)(i - start) * (MAX_SUPPORT + 1);
                        const struct poar_entry* e = NULL;
                        size_t n = 0;
                        if(poar_pair_entries_dense_levels(c->table, &dn, i, j, &sc, c->lo, c->hi, &e, &n) != OK){
                                c->error = 1;
                                poar_dense_free(c->table, &dn);
                                poar_scratch_free(&sc);
                                free(pos);
                                return;
                        }
                        for(size_t x = 0; x < n; x++){
                                int pc = popcount32(e[x].support);
                                if(pc >= c->lo && c->bucket[pc] != NULL){
                                        struct cand_pair* cp = &c->bucket[pc][p[pc]++];
                                        const int ei = c->seq_offsets[i] + POAR_KEY_I(e[x].key);
                                        const int ej = c->seq_offsets[j] + POAR_KEY_J(e[x].key);
                                        if(uf_find_ro(c->uf_parent, ei) == uf_find_ro(c->uf_parent, ej)){
                                                cp->elem_i = -1;
                                                cp->elem_j = -1;
                                        }else{
                                                cp->elem_i = ei;
                                                cp->elem_j = ej;
                                        }
                                }
                        }
                }
        }
        poar_dense_free(c->table, &dn);
        poar_scratch_free(&sc);
        free(pos);
}

/* One merge attempt.  Below min_support a pair may only attach a residue
   that is still in a column of its own; it never joins two columns. */
static int apply_candidate(struct uf_set* uf, struct topo_order* topo,
                            int elem_i, int elem_j, int support, int min_support)
{
        if(support < min_support){
                int ra = uf_find(uf, elem_i);
                int rb = uf_find(uf, elem_j);
                if(uf->next_in_set[uf->set_head[ra]] != -1 &&
                   uf->next_in_set[uf->set_head[rb]] != -1){
                        return 0;
                }
        }
        return uf_union_safe(uf, topo, elem_i, elem_j) < 0 ? -1 : 0;
}

int build_consensus(struct poar_table* table,
                    int* seq_lengths, int numseq,
                    int min_support,
                    int attach_support,
                    struct msa* out_msa
                    , threadpool_t* pool
                    )
{ KPERF_REGION(KP_ENS_CONSENSUS);
        struct uf_set* uf = NULL;
        int* seq_offsets = NULL;
        int* col_id = NULL;        /* element -> column */
        int* root_to_col = NULL;   /* uf root -> column index */
        struct topo_order* topo = NULL;
        struct ord_node* roots = NULL;
        int total_residues = 0;
        int n_cols = 0;
        int i, s, pos;
        char** out_seqs = NULL;
        int64_t* row_hist = NULL;
        int64_t* row_start = NULL;
        struct cand_pair* arena = NULL;

        ASSERT(table != NULL, "No POAR table");
        ASSERT(out_msa != NULL, "No output MSA");

        /* Candidates below min_support are only collected for the
           singleton-attachment pass (attach_support in [1, min_support)). */
        int collect_support = min_support;
        if(attach_support > 0 && attach_support < min_support){
                collect_support = attach_support;
        }

        /* Compute offsets */
        MMALLOC(seq_offsets, sizeof(int) * numseq);
        for(s = 0; s < numseq; s++){
                seq_offsets[s] = total_residues;
                total_residues += seq_lengths[s];
        }

        /* Initialize union-find with sequence membership tracking */
        RUN(uf_alloc(&uf, total_residues, seq_offsets, seq_lengths, numseq));

        /* Topological order of the column DAG, kept valid across merges */
        RUN(topo_alloc(&topo, total_residues));
        RUN(topo_init(topo, seq_offsets, seq_lengths, numseq, total_residues));

        /* Merge POAR entries from the highest support level down, so
           higher-confidence pairs merge first; within a level in (i asc,
           j asc, key asc) order.  Entries are derived from the table on the
           fly: one counting pass, then one pass per group of levels whose
           candidates fit in CAND_BUDGET (a larger level is split into blocks
           of consecutive rows). */
        {
                int64_t total[MAX_SUPPORT + 1];
                int64_t bucket_n[MAX_SUPPORT + 1];
                struct cand_pair* bucket[MAX_SUPPORT + 1];
                int64_t arena_n = 0;
                int64_t arena_cap = 0;   /* first allocation: min(budget, all candidates) */
                struct cand_ctx cc;
                int min_level = collect_support < 1 ? 1 : collect_support;
                int lv;

                for(lv = 0; lv <= MAX_SUPPORT; lv++){
                        total[lv] = 0;
                        bucket[lv] = NULL;
                }

                MMALLOC(row_hist, sizeof(int64_t) * (size_t)numseq * (MAX_SUPPORT + 1));
                cc.table = table;
                cc.numseq = numseq;
                cc.min_level = min_level;
                cc.lo = 0;
                cc.hi = MAX_SUPPORT;
                cc.seq_offsets = seq_offsets;
                cc.row_hist = row_hist;
                cc.row_start = NULL;
                cc.uf_parent = uf->parent;
                cc.bucket = bucket;
                cc.error = 0;
                if(table->row_hist != NULL && table->pairs == NULL){
                        /* counted with the run scores (score_runs_poar); levels below min_level are
                           in there too and are not read */
                        memcpy(row_hist, table->row_hist, sizeof(int64_t) * (size_t)numseq * (MAX_SUPPORT + 1));
                }else{
                        memset(row_hist, 0, sizeof(int64_t) * (size_t)numseq * (MAX_SUPPORT + 1));
                        FOR_ROWS(pool, numseq, 0, numseq - 1, cand_hist_chunk, &cc);
                        if(cc.error){
                                ERROR_MSG("POAR candidate counting failed");
                        }
                }
                for(i = 0; i < numseq; i++){
                        for(lv = min_level; lv <= MAX_SUPPORT; lv++){
                                total[lv] += row_hist[(size_t)i * (MAX_SUPPORT + 1) + lv];
                        }
                }

                for(lv = min_level; lv <= MAX_SUPPORT; lv++){
                        arena_cap += total[lv];
                }
                if(arena_cap > CAND_BUDGET) arena_cap = CAND_BUDGET;
                MMALLOC(row_start, sizeof(int64_t) * (size_t)numseq * (MAX_SUPPORT + 1));
                lv = MAX_SUPPORT;
                while(lv >= min_level){
                        int hi, lo, r0, r1;
                        int64_t acc;

                        if(total[lv] == 0){
                                lv--;
                                continue;
                        }

                        /* A group of consecutive levels is collected in one
                           pass over all rows if it fits; a level that does
                           not fit is processed in blocks of consecutive rows,
                           which keeps the (i, j, key) order. */
                        hi = lv;
                        acc = 0;
                        if(total[lv] <= CAND_BUDGET){
                                while(lv >= min_level && acc + total[lv] <= CAND_BUDGET){
                                        acc += total[lv];
                                        lv--;
                                }
                        }else{
                                lv--;
                        }
                        lo = lv + 1;

                        r0 = 0;
                        while(r0 < numseq - 1){
                                /* row block [r0, r1) */
                                r1 = r0;
                                if(lo == hi && total[hi] > CAND_BUDGET){
                                        acc = 0;
                                        do{
                                                acc += row_hist[(size_t)r1 * (MAX_SUPPORT + 1) + hi];
                                                r1++;
                                        }while(r1 < numseq - 1 &&
                                               acc + row_hist[(size_t)r1 * (MAX_SUPPORT + 1) + hi] <= CAND_BUDGET);
                                }else{
                                        r1 = numseq - 1;
                                }

                                {
                                        int64_t need = 0;
                                        int64_t used = 0;
                                        for(int v = lo; v <= hi; v++){
                                                int64_t run = 0;
                                                for(i = r0; i < r1; i++){
                                                        size_t idx = (size_t)i * (MAX_SUPPORT + 1) + v;
                                                        row_start[idx] = run;
                                                        run += row_hist[idx];
                                                }
                                                bucket_n[v] = run;
                                                need += run;
                                        }
                                        /* One candidate arena, reused by every pass */
                                        if(need > arena_n){
                                                int64_t cap = need;
                                                if(arena) MFREE(arena);
                                                if(cap < arena_cap) cap = arena_cap;
                                                MMALLOC(arena, sizeof(struct cand_pair) * (size_t)cap);
                                                arena_n = cap;
                                        }
                                        for(int v = lo; v <= hi; v++){
                                                bucket[v] = NULL;
                                                if(bucket_n[v] > 0){
                                                        bucket[v] = arena + used;
                                                        used += bucket_n[v];
                                                }
                                        }
                                }
                                cc.lo = lo;
                                cc.hi = hi;
                                cc.row_start = row_start;
                                FOR_ROWS(pool, numseq, r0, r1, cand_store_chunk, &cc);
                                if(cc.error){
                                        ERROR_MSG("POAR candidate collection failed");
                                }
                                for(int v = hi; v >= lo; v--){
                                        if(bucket[v] == NULL) continue;
                                        for(int64_t x = 0; x < bucket_n[v]; x++){
                                                if(bucket[v][x].elem_i < 0){
                                                        continue;   /* already one column at collection time */
                                                }
                                                RUN(apply_candidate(uf, topo, bucket[v][x].elem_i,
                                                                    bucket[v][x].elem_j, v, min_support));
                                        }
                                        bucket[v] = NULL;
                                }
                                r0 = r1;
                        }
                }

                MFREE(row_hist);
                if(row_start) MFREE(row_start);
                if(arena) MFREE(arena);
        }

        /* Number columns by their position in the topological order */
        MMALLOC(root_to_col, sizeof(int) * total_residues);
        MMALLOC(col_id, sizeof(int) * total_residues);
        MMALLOC(roots, sizeof(struct ord_node) * total_residues);

        for(i = 0; i < total_residues; i++){
                root_to_col[i] = -1;
        }
        n_cols = 0;
        for(i = 0; i < total_residues; i++){
                int root = uf_find(uf, i);
                if(root_to_col[root] == -1){
                        root_to_col[root] = 0;
                        roots[n_cols].ord = topo->ord[root];
                        roots[n_cols].node = root;
                        n_cols++;
                }
        }
        ord_sort(roots, n_cols);
        for(i = 0; i < n_cols; i++){
                root_to_col[roots[i].node] = i;
        }
        for(i = 0; i < total_residues; i++){
                col_id[i] = root_to_col[uf_find(uf, i)];
        }

        MFREE(roots);
        MFREE(root_to_col);
        topo_free(topo);
        topo = NULL;

        /* Every sequence must occupy strictly increasing columns */
        for(s = 0; s < numseq; s++){
                for(pos = 1; pos < seq_lengths[s]; pos++){
                        int elem = seq_offsets[s] + pos;
                        if(col_id[elem] <= col_id[elem - 1]){
                                ERROR_MSG("Consensus column order violates residue order in sequence %d", s);
                        }
                }
        }

        /* Build output alignment */
        MMALLOC(out_seqs, sizeof(char*) * numseq);
        for(s = 0; s < numseq; s++){
                out_seqs[s] = NULL;
        }
        for(s = 0; s < numseq; s++){
                MMALLOC(out_seqs[s], sizeof(char) * (n_cols + 1));
                memset(out_seqs[s], '-', n_cols);
                out_seqs[s][n_cols] = '\0';
        }

        /* Place residues */
        for(s = 0; s < numseq; s++){
                for(pos = 0; pos < seq_lengths[s]; pos++){
                        int elem = seq_offsets[s] + pos;
                        /* Use original residue character from the MSA */
                        out_seqs[s][col_id[elem]] = out_msa->sequences[s]->seq[pos];
                }
        }

        /* Write back into MSA: replace seq pointers */
        for(s = 0; s < numseq; s++){
                MFREE(out_msa->sequences[s]->seq);
                out_msa->sequences[s]->seq = out_seqs[s];
                out_msa->sequences[s]->len = n_cols;
                out_seqs[s] = NULL;
        }
        out_msa->alnlen = n_cols;
        out_msa->aligned = ALN_STATUS_FINAL;

        MFREE(out_seqs);
        MFREE(col_id);
        MFREE(seq_offsets);
        uf_free(uf);
        return OK;
ERROR:
        if(out_seqs){
                for(s = 0; s < numseq; s++){
                        if(out_seqs[s]) MFREE(out_seqs[s]);
                }
                MFREE(out_seqs);
        }
        if(col_id) MFREE(col_id);
        if(root_to_col) MFREE(root_to_col);
        if(roots) MFREE(roots);
        if(row_hist) MFREE(row_hist);
        if(row_start) MFREE(row_start);
        if(arena) MFREE(arena);
        topo_free(topo);
        if(seq_offsets) MFREE(seq_offsets);
        uf_free(uf);
        return FAIL;
}

/* Compute per-residue and per-column confidence from POAR table.
   For each residue at alignment position pos in sequence i:
   - Count residue-residue pairs with other sequences j that also have
     a residue (not gap) at the same column.
   - For each such pair, look up POAR support (popcount).
   - confidence = sum(supports) / (n_residue_pairs * n_alignments)
   Gaps get confidence 0.0.
   Column confidence = mean of residue confidences in that column. */
/* 1 if every residue pair of the table is a plain column-equality test
   (implicit table). */
static int table_all_exact(const struct poar_table* table, int numseq)
{
        (void)numseq;
        return table->pairs == NULL && table->off != NULL;
}

struct confidence_ctx {
        struct poar_table* table;
        struct pos_matrix* pm;
        struct msa* msa;
        int numseq;
        int alnlen;
        int n_alignments;
        int all_exact;
        int error;
};

/* Number of runs aligning the residue at table index ia with the one at ib,
   for every b after a: isum[a] and isum[b] both gain it. */
static inline __attribute__((always_inline))
void conf_pairs_fast(const struct poar_table* table, const size_t* idx, int m,
                     int* isum, const int NR)
{
        const int32_t* cr[POAR_MAX_RUNS];
        for(int r = 0; r < NR; r++){
                cr[r] = table->col[table->run_ids[r]];
        }
        for(int a = 0; a < m - 1; a++){
                int32_t va[POAR_MAX_RUNS];
                int acc = 0;
                for(int r = 0; r < NR; r++){
                        va[r] = cr[r][idx[a]];
                }
                for(int b = a + 1; b < m; b++){
                        int sup = 0;
                        for(int r = 0; r < NR; r++){
                                sup += (va[r] == cr[r][idx[b]]);
                        }
                        acc += sup;
                        isum[b] += sup;
                }
                isum[a] += acc;
        }
}

/* Per-residue confidence of the alignment columns [start, end): the
   residues of a column are paired among themselves once. */
static void confidence_chunk(int start, int end, void* arg)
{
        struct confidence_ctx* c = (struct confidence_ctx*)arg;
        const struct poar_table* table = c->table;
        int* mj = NULL;
        int* mr = NULL;
        int* isum = NULL;
        size_t* idx = NULL;

        mj = malloc(sizeof(int) * (size_t)c->numseq);
        mr = malloc(sizeof(int) * (size_t)c->numseq);
        isum = malloc(sizeof(int) * (size_t)c->numseq);
        idx = malloc(sizeof(size_t) * (size_t)c->numseq);
        if(!mj || !mr || !isum || !idx){
                c->error = 1;
                goto done;
        }
        for(int col = start; col < end; col++){
                int m = 0;
                for(int j = 0; j < c->numseq; j++){
                        int rj = c->pm->col_to_res[j][col];
                        if(rj >= 0){
                                mj[m] = j;
                                mr[m] = rj;
                                m++;
                        }
                }
                for(int a = 0; a < m; a++){
                        isum[a] = 0;
                }
                if(c->all_exact){
                        for(int a = 0; a < m; a++){
                                idx[a] = table->off[mj[a]] + (size_t)mr[a];
                        }
                        switch(table->n_present){
                        case 1: conf_pairs_fast(table, idx, m, isum, 1); break;
                        case 2: conf_pairs_fast(table, idx, m, isum, 2); break;
                        case 3: conf_pairs_fast(table, idx, m, isum, 3); break;
                        case 4: conf_pairs_fast(table, idx, m, isum, 4); break;
                        case 5: conf_pairs_fast(table, idx, m, isum, 5); break;
                        default: conf_pairs_fast(table, idx, m, isum, table->n_present); break;
                        }
                }else{
                        for(int a = 0; a < m - 1; a++){
                                for(int b = a + 1; b < m; b++){
                                        int sup = popcount32(poar_support(table, mj[a], mr[a], mj[b], mr[b]));
                                        isum[a] += sup;
                                        isum[b] += sup;
                                }
                        }
                }
                for(int a = 0; a < m; a++){
                        if(m > 1 && c->n_alignments > 0){
                                c->msa->sequences[mj[a]]->confidence[col] =
                                        (float)((double)isum[a] / ((double)(m - 1) * (double)c->n_alignments));
                        }else{
                                c->msa->sequences[mj[a]]->confidence[col] = 0.0f;
                        }
                }
        }
done:
        free(mj);
        free(mr);
        free(isum);
        free(idx);
}

int compute_residue_confidence(struct poar_table* table,
                               struct msa* aligned_msa
                               , threadpool_t* pool
                               )
{
        struct pos_matrix* pm = NULL;
        int numseq = aligned_msa->numseq;
        int alnlen = aligned_msa->alnlen;
        int n_alignments = table->n_alignments;
        int i, col;
        char** seqs = NULL;

        ASSERT(table != NULL, "No POAR table");
        ASSERT(aligned_msa != NULL, "No aligned MSA");
        ASSERT(alnlen > 0, "Alignment length is 0");

        /* Build position matrix from aligned MSA */
        MMALLOC(seqs, sizeof(char*) * numseq);
        for(i = 0; i < numseq; i++){
                seqs[i] = aligned_msa->sequences[i]->seq;
        }
        RUN(pos_matrix_from_msa(&pm, seqs, numseq, alnlen));

        /* Allocate per-residue confidence arrays */
        for(i = 0; i < numseq; i++){
                if(aligned_msa->sequences[i]->confidence){
                        MFREE(aligned_msa->sequences[i]->confidence);
                }
                aligned_msa->sequences[i]->confidence = NULL;
                MMALLOC(aligned_msa->sequences[i]->confidence,
                        sizeof(float) * alnlen);
                for(col = 0; col < alnlen; col++){
                        aligned_msa->sequences[i]->confidence[col] = 0.0f;
                }
        }

        /* Allocate per-column confidence */
        if(aligned_msa->col_confidence){
                MFREE(aligned_msa->col_confidence);
        }
        aligned_msa->col_confidence = NULL;
        MMALLOC(aligned_msa->col_confidence, sizeof(float) * alnlen);

        /* Compute per-residue confidence */
        {
                struct confidence_ctx ctx = {
                        table, pm, aligned_msa, numseq, alnlen, n_alignments,
                        table_all_exact(table, numseq), 0
                };
                if(pool != NULL && numseq > 16){
                        tp_parallel_for(pool, 0, alnlen, confidence_chunk, &ctx);
                }else
                {
                        confidence_chunk(0, alnlen, &ctx);
                }
                if(ctx.error){
                        ERROR_MSG("POAR confidence allocation failed");
                }
        }

        /* Compute per-column confidence: mean over non-gap residues */
        for(col = 0; col < alnlen; col++){
                double sum = 0.0;
                int count = 0;
                for(i = 0; i < numseq; i++){
                        if(pm->col_to_res[i][col] >= 0){
                                sum += aligned_msa->sequences[i]->confidence[col];
                                count++;
                        }
                }
                if(count > 0){
                        aligned_msa->col_confidence[col] = (float)(sum / count);
                }else{
                        aligned_msa->col_confidence[col] = 0.0f;
                }
        }

        pos_matrix_free(pm);
        MFREE(seqs);
        return OK;
ERROR:
        if(pm) pos_matrix_free(pm);
        if(seqs) MFREE(seqs);
        return FAIL;
}

/* Score one alignment against the POAR table.
   Returns expected number of correct pairs: for each aligned pair,
   (support - 1) / (m - 1) gives the fraction of OTHER alignments
   agreeing. Summing these gives expected correct pairs.
   This rewards both high recall (many pairs) and high precision
   (pairs with broad agreement). */
/* Add (support - 1) / denom for every column where both sequences have a
   residue, in column order, to *acc.  Supports come from the run column
   maps directly (implicit table); a stored table's entries are scanned with
   a forward cursor. */
static int score_pair(const struct poar_table* table, const struct pos_matrix* pm,
                      int i, int j, double denom,
                      struct poar_scratch* sc, double* acc)
{
        double sum = *acc;
        int alnlen = pm->alnlen;
        int col;

        if(poar_pair_is_exact(table, i, j)){
                for(col = 0; col < alnlen; col++){
                        int ri = pm->col_to_res[i][col];
                        int rj = pm->col_to_res[j][col];
                        if(ri >= 0 && rj >= 0){
                                int support = popcount32(poar_support(table, i, ri, j, rj));
                                sum += (double)(support - 1) / denom;
                        }
                }
        }else{
                const struct poar_entry* e = NULL;
                size_t n = 0;
                size_t cur = 0;
                RUN(poar_pair_entries(table, i, j, sc, &e, &n));
                for(col = 0; col < alnlen; col++){
                        int ri = pm->col_to_res[i][col];
                        int rj = pm->col_to_res[j][col];
                        if(ri >= 0 && rj >= 0){
                                uint64_t key = POAR_KEY(ri, rj);
                                int support = 0;
                                while(cur < n && e[cur].key < key) cur++;
                                if(cur < n && e[cur].key == key){
                                        support = popcount32(e[cur].support);
                                }
                                sum += (double)(support - 1) / denom;
                        }
                }
        }
        *acc = sum;
        return OK;
ERROR:
        return FAIL;
}

/* Same terms as the exact branch of score_pair, in the same order, walking the
   residues of i (rc_i: residue -> column in the scored alignment) instead of
   all columns; tbl[s] = (s - 1) / denom. */
static inline __attribute__((always_inline))
void score_pair_fast(const struct poar_table* table, const struct pos_matrix* pm,
                     const int* rc_i, size_t li, int i, int j,
                     const double* tbl, double* acc, const int NR)
{
        const int32_t* ci[POAR_MAX_RUNS];
        const int32_t* cj[POAR_MAX_RUNS];
        const int* colres_j = pm->col_to_res[j];
        double sum = *acc;
        for(int r = 0; r < NR; r++){
                int k = table->run_ids[r];
                ci[r] = table->col[k] + table->off[i];
                cj[r] = table->col[k] + table->off[j];
        }
        for(size_t x = 0; x < li; x++){
                int rj = colres_j[rc_i[x]];
                if(rj >= 0){
                        int sup = 0;
                        for(int r = 0; r < NR; r++){
                                sup += (ci[r][x] == cj[r][rj]);
                        }
                        sum += tbl[sup];
                }
        }
        *acc = sum;
}

struct score_fast {
        const int* rc;              /* residue -> column of the scored alignment, by table->off */
        double tbl[POAR_MAX_RUNS + 1];
};

static int score_pair_dispatch(const struct poar_table* table, const struct pos_matrix* pm,
                               const struct score_fast* sf,
                               int i, int j, double denom,
                               struct poar_scratch* sc, double* acc)
{
        if(sf != NULL && poar_pair_is_exact(table, i, j)){
                const int* rc_i = sf->rc + table->off[i];
                size_t li = table->off[i + 1] - table->off[i];
                switch(table->n_present){
                case 1: score_pair_fast(table, pm, rc_i, li, i, j, sf->tbl, acc, 1); break;
                case 2: score_pair_fast(table, pm, rc_i, li, i, j, sf->tbl, acc, 2); break;
                case 3: score_pair_fast(table, pm, rc_i, li, i, j, sf->tbl, acc, 3); break;
                case 4: score_pair_fast(table, pm, rc_i, li, i, j, sf->tbl, acc, 4); break;
                case 5: score_pair_fast(table, pm, rc_i, li, i, j, sf->tbl, acc, 5); break;
                default: score_pair_fast(table, pm, rc_i, li, i, j, sf->tbl, acc, table->n_present); break;
                }
                return OK;
        }
        return score_pair(table, pm, i, j, denom, sc, acc);
}

/* ---- all runs of the table scored in one pass ----
   score_alignment_poar on run k adds, for every sequence pair (i, j) and every residue of i that run k
   aligns to a residue of j, tbl[number of runs with that same pairing], in ascending residue order.  Here
   sequence i is the loaded one and the residues y of j are walked: x[r] is the residue of i that run r puts
   in the column of y.  For a fixed run the aligned pairs of (i, j) rise in both sequences, so walking y
   visits them in the same order as walking the residues of i, and every run's sum receives the same terms
   in the same order: the scores are the same doubles, bit for bit. */
static inline __attribute__((always_inline))
void score_runs_pair(const struct poar_table* table, const struct poar_dense* dn,
                     int j, const double* tbl, double* acc, int64_t* hist, const int NR)
{
        const int32_t* cj[POAR_MAX_RUNS];
        const int32_t* d[POAR_MAX_RUNS];
        double a[POAR_MAX_RUNS];
        int64_t h[POAR_MAX_RUNS + 1];
        const size_t lj = table->off[j + 1] - table->off[j];
        if(NR < 1 || NR > POAR_MAX_RUNS){
                return;                 /* no runs: nothing to score (score_runs_poar rejects such a table) */
        }
        const double v_all = tbl[NR];
        for(int r = 0; r < NR; r++){
                int k = table->run_ids[r];
                cj[r] = table->col[k] + table->off[j];
                d[r] = dn->d[k];
                a[r] = acc[r];
        }
        for(int r = 0; r <= NR; r++){
                h[r] = 0;
        }
        for(size_t y = 0; y < lj; y++){
                int32_t x[POAR_MAX_RUNS];
                int32_t diff = 0;
                for(int r = 0; r < NR; r++){
                        x[r] = d[r][cj[r][y]];
                }
                for(int r = 1; r < NR; r++){
                        diff |= x[r] ^ x[0];
                }
                if(diff == 0){
                        if(x[0] >= 0){
                                for(int r = 0; r < NR; r++){
                                        a[r] += v_all;
                                }
                                h[NR]++;
                        }
                }else{
                        for(int r = 0; r < NR; r++){
                                if(x[r] >= 0){
                                        int sup = 0;
                                        int first = 1;
                                        for(int q = 0; q < NR; q++){
                                                const int same = (x[q] == x[r]);
                                                sup += same;
                                                first &= !(same && q < r);
                                        }
                                        a[r] += tbl[sup];
                                        h[sup] += first;        /* one entry per distinct partner */
                                }
                        }
                }
        }
        for(int r = 0; r < NR; r++){
                acc[r] = a[r];
        }
        for(int r = 1; r <= NR; r++){
                hist[r] += h[r];
        }
}

/* adds the terms of the pairs (i, j), j = i + 1 .. numseq - 1, to acc[r] (run table->run_ids[r]) and counts
   the entries of these pairs by the number of runs that support them: hist[s], s = 1 .. n_present */
static void score_runs_row(const struct poar_table* table, struct poar_dense* dn,
                           int i, int numseq, const double* tbl, double* acc, int64_t* hist)
{
        poar_dense_load(table, dn, i);
        for(int j = i + 1; j < numseq; j++){
                switch(table->n_present){
                case 1: score_runs_pair(table, dn, j, tbl, acc, hist, 1); break;
                case 2: score_runs_pair(table, dn, j, tbl, acc, hist, 2); break;
                case 3: score_runs_pair(table, dn, j, tbl, acc, hist, 3); break;
                case 4: score_runs_pair(table, dn, j, tbl, acc, hist, 4); break;
                case 5: score_runs_pair(table, dn, j, tbl, acc, hist, 5); break;
                default: score_runs_pair(table, dn, j, tbl, acc, hist, table->n_present); break;
                }
        }
}

struct score_runs_ctx {
        const struct poar_table* table;
        const double* tbl;
        double* row_scores;         /* [numseq][n_present] */
        int64_t* row_hist;          /* [numseq][POAR_MAX_RUNS + 1] */
        int numseq;
        int error;
};

static void score_runs_chunk(int start, int end, void* arg)
{
        struct score_runs_ctx* c = (struct score_runs_ctx*)arg;
        struct poar_dense dn;
        if(poar_dense_init(c->table, &dn) != OK){
                c->error = 1;
                return;
        }
        for(int i = start; i < end; i++){
                score_runs_row(c->table, &dn, i, c->numseq, c->tbl,
                               c->row_scores + (size_t)i * (size_t)c->table->n_present,
                               c->row_hist + (size_t)i * (POAR_MAX_RUNS + 1));
        }
        poar_dense_free(c->table, &dn);
}

int poar_run_is_alignment(const struct poar_table* table, const struct pos_matrix* pm, int k)
{
        if(table->pairs || table->off == NULL || k < 0 || k >= POAR_MAX_RUNS || table->col[k] == NULL){
                return 0;
        }
        if(pm->numseq != table->numseq || pm->alnlen != table->alnlen[k]){
                return 0;
        }
        for(int i = 0; i < pm->numseq; i++){
                const int32_t* col = table->col[k] + table->off[i];
                size_t li = table->off[i + 1] - table->off[i];
                size_t n = 0;
                for(int c = 0; c < pm->alnlen; c++){
                        if(pm->col_to_res[i][c] >= 0){
                                if(n >= li || col[n] != c){
                                        return 0;
                                }
                                n++;
                        }
                }
                if(n != li){
                        return 0;
                }
        }
        return 1;
}

int score_runs_poar(struct poar_table* table, int numseq, int n_alignments,
                    double* scores, threadpool_t* pool)
{
        double denom = (n_alignments > 1) ? (double)(n_alignments - 1) : 1.0;
        double tbl[POAR_MAX_RUNS + 1];
        double* row_scores = NULL;
        int64_t* row_hist = NULL;
        const int NR = table->n_present;

        ASSERT(table->pairs == NULL && table->off != NULL, "score_runs_poar needs an implicit POAR table");
        ASSERT(numseq == table->numseq, "score_runs_poar: %d sequences, table %d", numseq, table->numseq);
        ASSERT(NR >= 1 && NR <= POAR_MAX_RUNS, "score_runs_poar: table holds %d runs (1 to %d expected)", NR, POAR_MAX_RUNS);
        for(int sup = 0; sup <= POAR_MAX_RUNS; sup++){
                tbl[sup] = (double)(sup - 1) / denom;
        }
        /* the same pass counts every row's entries per support level; build_consensus reads the counts
           from the table instead of enumerating the pairs once more */
        MMALLOC(row_hist, sizeof(int64_t) * (size_t)numseq * (POAR_MAX_RUNS + 1));
        memset(row_hist, 0, sizeof(int64_t) * (size_t)numseq * (POAR_MAX_RUNS + 1));
        /* same association as score_alignment_poar_run: row sums added up with a pool, one running sum
           without */
        if(pool != NULL && numseq > 16){
                struct score_runs_ctx ctx;
                MMALLOC(row_scores, sizeof(double) * (size_t)numseq * (size_t)(NR > 0 ? NR : 1));
                for(size_t x = 0; x < (size_t)numseq * (size_t)NR; x++){
                        row_scores[x] = 0.0;
                }
                ctx.table = table;
                ctx.tbl = tbl;
                ctx.row_scores = row_scores;
                ctx.row_hist = row_hist;
                ctx.numseq = numseq;
                ctx.error = 0;
                tp_parallel_for(pool, 0, numseq - 1, score_runs_chunk, &ctx);
                if(ctx.error){
                        ERROR_MSG("POAR scoring failed");
                }
                for(int r = 0; r < NR; r++){
                        double total = 0.0;
                        for(int i = 0; i < numseq - 1; i++){
                                total += row_scores[(size_t)i * (size_t)NR + r];
                        }
                        scores[table->run_ids[r]] = total;
                }
                MFREE(row_scores);
        }else{
                struct poar_dense dn;
                double acc[POAR_MAX_RUNS];
                for(int r = 0; r < POAR_MAX_RUNS; r++){
                        acc[r] = 0.0;
                }
                RUN(poar_dense_init(table, &dn));
                for(int i = 0; i < numseq - 1; i++){
                        score_runs_row(table, &dn, i, numseq, tbl, acc,
                                       row_hist + (size_t)i * (POAR_MAX_RUNS + 1));
                }
                poar_dense_free(table, &dn);
                for(int r = 0; r < NR; r++){
                        scores[table->run_ids[r]] = acc[r];
                }
        }
        if(table->row_hist){
                MFREE(table->row_hist);
        }
        table->row_hist = row_hist;
        return OK;
ERROR:
        if(row_scores) MFREE(row_scores);
        if(row_hist) MFREE(row_hist);
        return FAIL;
}

struct score_poar_ctx {
        const struct poar_table* table;
        const struct pos_matrix* pm;
        const struct score_fast* sf;
        int numseq;
        double denom;
        double* row_scores;
        int error;
};

static void score_poar_chunk(int start, int end, void* arg)
{
        struct score_poar_ctx* c = (struct score_poar_ctx*)arg;
        struct poar_scratch sc;
        poar_scratch_init(&sc);
        for(int i = start; i < end; i++){
                double row_sum = 0.0;
                for(int j = i + 1; j < c->numseq; j++){
                        if(score_pair_dispatch(c->table, c->pm, c->sf, i, j, c->denom, &sc, &row_sum) != OK){
                                c->error = 1;
                                poar_scratch_free(&sc);
                                return;
                        }
                }
                c->row_scores[i] = row_sum;
        }
        poar_scratch_free(&sc);
}

static int score_alignment_poar_run(const struct poar_table* table,
                         const struct pos_matrix* pm,
                         const struct score_fast* sfp,
                         int numseq,
                         int n_alignments,
                         double* out_score
                         , threadpool_t* pool
                         )
{
        double total_score = 0.0;
        int i;
        double denom = (n_alignments > 1) ? (double)(n_alignments - 1) : 1.0;

        if(pool != NULL && numseq > 16){
                double* row_scores = NULL;
                MMALLOC(row_scores, sizeof(double) * numseq);
                for(i = 0; i < numseq; i++) row_scores[i] = 0.0;

                struct score_poar_ctx ctx = {
                        table, pm, sfp, numseq, denom, row_scores, 0
                };
                tp_parallel_for(pool, 0, numseq - 1, score_poar_chunk, &ctx);
                if(ctx.error){
                        ERROR_MSG("POAR scoring failed");
                }

                for(i = 0; i < numseq - 1; i++){
                        total_score += row_scores[i];
                }
                MFREE(row_scores);

                *out_score = total_score;
                return OK;
        ERROR:
                if(row_scores) MFREE(row_scores);
                return FAIL;
        }
        {
                struct poar_scratch sc;
                int j;
                poar_scratch_init(&sc);
                for(i = 0; i < numseq - 1; i++){
                        for(j = i + 1; j < numseq; j++){
                                if(score_pair_dispatch(table, pm, sfp, i, j, denom, &sc, &total_score) != OK){
                                        poar_scratch_free(&sc);
                                        return FAIL;
                                }
                        }
                }
                poar_scratch_free(&sc);
        }

        *out_score = total_score;
        return OK;
}

int score_alignment_poar(struct poar_table* table,
                         struct pos_matrix* pm,
                         int numseq,
                         int n_alignments,
                         double* out_score
                         , threadpool_t* pool
                         )
{
        double denom = (n_alignments > 1) ? (double)(n_alignments - 1) : 1.0;
        struct score_fast sf;
        const struct score_fast* sfp = NULL;
        int* rc = NULL;
        int i;
        int ret;

        /* residue -> column map of the scored alignment */
        if(!table->pairs && table->off != NULL && table->n_present <= POAR_MAX_RUNS &&
           pm->numseq == numseq){
                size_t total = table->off[numseq];
                int ok = 1;
                rc = malloc(sizeof(int) * (total > 0 ? total : 1));
                if(rc == NULL){
                        return FAIL;
                }
                for(i = 0; i < numseq && ok; i++){
                        size_t n = 0;
                        size_t li = table->off[i + 1] - table->off[i];
                        int* dst = rc + table->off[i];
                        for(int col = 0; col < pm->alnlen; col++){
                                if(pm->col_to_res[i][col] >= 0){
                                        if(n >= li){
                                                ok = 0;
                                                break;
                                        }
                                        dst[n++] = col;
                                }
                        }
                        if(n != li){
                                ok = 0;
                        }
                }
                if(ok){
                        sf.rc = rc;
                        for(int sup = 0; sup <= POAR_MAX_RUNS; sup++){
                                sf.tbl[sup] = (double)(sup - 1) / denom;
                        }
                        sfp = &sf;
                }
        }

        ret = score_alignment_poar_run(table, pm, sfp, numseq, n_alignments, out_score
                                       , pool
                                       );
        free(rc);
        return ret;
}

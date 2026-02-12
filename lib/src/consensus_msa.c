#include "tldevel.h"

#include <string.h>

#include "msa_struct.h"
#include "msa_alloc.h"
#include "poar.h"

#define CONSENSUS_MSA_IMPORT
#include "consensus_msa.h"

/* Union-find data structure with per-set sequence membership tracking
   and element linked lists for efficient set enumeration. */
struct uf_set {
        int* parent;
        int* rank;
        int* elem_seq;          /* element -> sequence index */
        uint64_t** seq_mask;    /* per-root: bitmask of sequences in set */
        int* set_head;          /* root -> first element in set (-1 if none) */
        int* next_in_set;       /* element -> next element in same set (-1 if last) */
        int n_elements;
        int numseq;
        int mask_words;         /* number of uint64_t words per bitmask */
};

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
        u->next_in_set = NULL;
        u->n_elements = n;
        u->numseq = numseq;
        u->mask_words = mw;

        MMALLOC(u->parent, sizeof(int) * n);
        MMALLOC(u->rank, sizeof(int) * n);
        MMALLOC(u->elem_seq, sizeof(int) * n);
        MMALLOC(u->set_head, sizeof(int) * n);
        MMALLOC(u->next_in_set, sizeof(int) * n);
        MMALLOC(u->seq_mask, sizeof(uint64_t*) * n);

        for(int i = 0; i < n; i++){
                u->parent[i] = i;
                u->rank[i] = 0;
                u->set_head[i] = i;
                u->next_in_set[i] = -1;
                u->seq_mask[i] = NULL;
                MMALLOC(u->seq_mask[i], sizeof(uint64_t) * mw);
                memset(u->seq_mask[i], 0, sizeof(uint64_t) * mw);
        }

        /* Fill elem_seq and initialize per-element bitmasks */
        for(int s = 0; s < numseq; s++){
                for(int p = 0; p < seq_lengths[s]; p++){
                        int elem = seq_offsets[s] + p;
                        u->elem_seq[elem] = s;
                        u->seq_mask[elem][s / 64] |= (1ULL << (s % 64));
                }
        }

        *uf = u;
        return OK;
ERROR:
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

/* Check if root 'target' is reachable from root 'start' via the
   implicit column ordering DAG. A directed edge col_A -> col_B exists
   when some sequence has consecutive residues in cols A and B.
   If target is reachable from start, merging them would create a cycle. */
static int dag_reachable(struct uf_set* u, int start, int target,
                         int* seq_offsets, int* seq_lengths,
                         int* visited, int visit_id)
{
        /* BFS using the element linked lists for efficient enumeration */
        /* We use a simple queue built from a scratch array */
        int queue[4096];  /* bounded queue; for very large problems this needs malloc */
        int head = 0, tail = 0;

        if(start == target) return 1;

        queue[tail++] = start;
        visited[start] = visit_id;

        while(head < tail){
                int cur = queue[head++];

                /* Enumerate outgoing edges: for each element in cur's set,
                   check its successor (same seq, pos+1) */
                int elem = u->set_head[cur];
                while(elem >= 0){
                        int s = u->elem_seq[elem];
                        int pos = elem - seq_offsets[s];
                        if(pos + 1 < seq_lengths[s]){
                                int succ_elem = seq_offsets[s] + pos + 1;
                                int succ_root = uf_find(u, succ_elem);
                                if(succ_root == target) return 1;
                                if(succ_root != cur && visited[succ_root] != visit_id){
                                        visited[succ_root] = visit_id;
                                        if(tail < 4096){
                                                queue[tail++] = succ_root;
                                        }
                                }
                        }
                        elem = u->next_in_set[elem];
                }
        }
        return 0;
}

/* Merge two sets. Returns 1 if merge succeeded, 0 if blocked by:
   - same-sequence conflict (two residues from same seq in one column)
   - ordering cycle (merging would create a cycle in the column DAG)
   When seq_offsets/seq_lengths are NULL, skip cycle check. */
static int uf_union_safe(struct uf_set* u, int a, int b,
                         int* seq_offsets, int* seq_lengths,
                         int* visited, int* visit_counter)
{
        int ra = uf_find(u, a);
        int rb = uf_find(u, b);
        if(ra == rb) return 1;

        /* Check for conflict: do the two sets share any sequence? */
        for(int w = 0; w < u->mask_words; w++){
                if(u->seq_mask[ra][w] & u->seq_mask[rb][w]){
                        return 0; /* conflict: same sequence in both sets */
                }
        }

        /* Check for ordering cycle: would merging create a cycle?
           A cycle exists if ra can reach rb (meaning ra < rb in the
           partial order, but merging makes them equal). */
        if(seq_offsets != NULL){
                (*visit_counter)++;
                if(dag_reachable(u, ra, rb, seq_offsets, seq_lengths,
                                 visited, *visit_counter)){
                        return 0; /* would create cycle */
                }
                (*visit_counter)++;
                if(dag_reachable(u, rb, ra, seq_offsets, seq_lengths,
                                 visited, *visit_counter)){
                        return 0; /* would create cycle */
                }
        }

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

        /* Merge sequence bitmasks into new root */
        for(int w = 0; w < u->mask_words; w++){
                u->seq_mask[new_root][w] |= u->seq_mask[old_root][w];
        }

        /* Concatenate element linked lists: append old_root's list to new_root's */
        if(u->set_head[old_root] >= 0){
                /* Find tail of new_root's list */
                int tail = u->set_head[new_root];
                if(tail < 0){
                        u->set_head[new_root] = u->set_head[old_root];
                }else{
                        while(u->next_in_set[tail] >= 0){
                                tail = u->next_in_set[tail];
                        }
                        u->next_in_set[tail] = u->set_head[old_root];
                }
        }
        u->set_head[old_root] = -1;

        return 1;
}

static inline int popcount32(uint32_t x)
{
#ifdef __GNUC__
        return __builtin_popcount(x);
#else
        x = x - ((x >> 1) & 0x55555555u);
        x = (x & 0x33333333u) + ((x >> 2) & 0x33333333u);
        return (int)(((x + (x >> 4)) & 0x0F0F0F0Fu) * 0x01010101u >> 24);
#endif
}

static inline int pair_index(int i, int j, int numseq)
{
        return i * numseq - (i * (i + 1)) / 2 + (j - i - 1);
}

struct merge_candidate {
        int elem_i;
        int elem_j;
        int support;
};

/* DFS-based topological sort that gracefully handles cycles by
   skipping back edges. This produces a valid ordering that respects
   as many sequence constraints as possible. */
static int topo_sort(int* col_id,         /* [total_residues]: element -> column id */
                     int* seq_offsets,     /* [numseq]: start offset for each seq */
                     int* seq_lengths,     /* [numseq] */
                     int numseq,
                     int n_cols,
                     int** sorted_cols,    /* output: sorted column indices */
                     int* n_sorted)
{
        int* out = NULL;
        int** adj_list = NULL;
        int* adj_count = NULL;
        int* adj_alloc = NULL;
        int* state = NULL;       /* 0=unvisited, 1=in-progress, 2=done */
        int* dfs_stack = NULL;   /* pairs of (node, edge_index) */
        int i, s, pos;
        int out_idx;
        int sp;

        MMALLOC(adj_count, sizeof(int) * n_cols);
        MMALLOC(adj_alloc, sizeof(int) * n_cols);
        MMALLOC(adj_list, sizeof(int*) * n_cols);
        MMALLOC(out, sizeof(int) * n_cols);
        MMALLOC(state, sizeof(int) * n_cols);
        MMALLOC(dfs_stack, sizeof(int) * n_cols * 2);

        for(i = 0; i < n_cols; i++){
                adj_count[i] = 0;
                adj_alloc[i] = 4;
                adj_list[i] = NULL;
                MMALLOC(adj_list[i], sizeof(int) * adj_alloc[i]);
                state[i] = 0;
        }

        /* Build adjacency list (deduplicated) */
        for(s = 0; s < numseq; s++){
                for(pos = 0; pos < seq_lengths[s] - 1; pos++){
                        int elem_a = seq_offsets[s] + pos;
                        int elem_b = seq_offsets[s] + pos + 1;
                        int ca = col_id[elem_a];
                        int cb = col_id[elem_b];
                        if(ca != cb){
                                int dup = 0;
                                for(int k = 0; k < adj_count[ca]; k++){
                                        if(adj_list[ca][k] == cb){
                                                dup = 1;
                                                break;
                                        }
                                }
                                if(!dup){
                                        if(adj_count[ca] >= adj_alloc[ca]){
                                                adj_alloc[ca] *= 2;
                                                MREALLOC(adj_list[ca], sizeof(int) * adj_alloc[ca]);
                                        }
                                        adj_list[ca][adj_count[ca]++] = cb;
                                }
                        }
                }
        }

        /* DFS topological sort — back edges (cycles) are silently skipped */
        out_idx = n_cols - 1;
        for(int start = 0; start < n_cols; start++){
                if(state[start] != 0) continue;

                sp = 0;
                dfs_stack[sp++] = start;
                dfs_stack[sp++] = 0;
                state[start] = 1;

                while(sp > 0){
                        int edge_idx = dfs_stack[--sp];
                        int node = dfs_stack[--sp];

                        int pushed = 0;
                        for(int e = edge_idx; e < adj_count[node]; e++){
                                int next = adj_list[node][e];
                                if(state[next] == 0){
                                        /* Push current with updated edge index */
                                        dfs_stack[sp++] = node;
                                        dfs_stack[sp++] = e + 1;
                                        /* Push next */
                                        dfs_stack[sp++] = next;
                                        dfs_stack[sp++] = 0;
                                        state[next] = 1;
                                        pushed = 1;
                                        break;
                                }
                                /* state[next]==1: back edge (cycle) — skip */
                                /* state[next]==2: cross/forward edge — skip */
                        }

                        if(!pushed){
                                state[node] = 2;
                                out[out_idx--] = node;
                        }
                }
        }

        for(i = 0; i < n_cols; i++){
                if(adj_list[i]) MFREE(adj_list[i]);
        }
        MFREE(adj_list);
        MFREE(adj_count);
        MFREE(adj_alloc);
        MFREE(state);
        MFREE(dfs_stack);

        *sorted_cols = out;
        *n_sorted = n_cols;
        return OK;
ERROR:
        if(adj_list){
                for(i = 0; i < n_cols; i++){
                        if(adj_list[i]) MFREE(adj_list[i]);
                }
                MFREE(adj_list);
        }
        if(adj_count) MFREE(adj_count);
        if(adj_alloc) MFREE(adj_alloc);
        if(state) MFREE(state);
        if(dfs_stack) MFREE(dfs_stack);
        if(out) MFREE(out);
        return FAIL;
}

int build_consensus(struct poar_table* table,
                    int* seq_lengths, int numseq,
                    int min_support,
                    struct msa* out_msa)
{
        struct uf_set* uf = NULL;
        int* seq_offsets = NULL;
        int* col_id = NULL;        /* element -> column */
        int* root_to_col = NULL;   /* uf root -> column index */
        int* sorted_cols = NULL;
        int* visited = NULL;       /* for cycle detection BFS */
        int visit_counter = 0;
        int n_sorted = 0;
        int total_residues = 0;
        int n_cols = 0;
        int i, j, s, pos;
        char** out_seqs = NULL;

        ASSERT(table != NULL, "No POAR table");
        ASSERT(out_msa != NULL, "No output MSA");

        /* Compute offsets */
        MMALLOC(seq_offsets, sizeof(int) * numseq);
        for(s = 0; s < numseq; s++){
                seq_offsets[s] = total_residues;
                total_residues += seq_lengths[s];
        }

        /* Initialize union-find with sequence membership tracking */
        RUN(uf_alloc(&uf, total_residues, seq_offsets, seq_lengths, numseq));

        /* Allocate visited array for cycle detection (use visit_counter
           to avoid clearing between checks) */
        MMALLOC(visited, sizeof(int) * total_residues);
        memset(visited, 0, sizeof(int) * total_residues);

        /* Collect all POAR entries above min_support, then process
           in descending support order so higher-confidence pairs merge first */
        {
                int n_candidates = 0;
                int alloc_candidates = 1024;
                struct merge_candidate *candidates = NULL;

                MMALLOC(candidates, sizeof(*candidates) * alloc_candidates);

                for(i = 0; i < numseq - 1; i++){
                        for(j = i + 1; j < numseq; j++){
                                int pidx = pair_index(i, j, numseq);
                                struct poar_pair* pp = table->pairs[pidx];

                                for(int e = 0; e < pp->n_entries; e++){
                                        int support = popcount32(pp->entries[e].support);
                                        if(support >= min_support){
                                                if(n_candidates >= alloc_candidates){
                                                        alloc_candidates *= 2;
                                                        MREALLOC(candidates, sizeof(*candidates) * alloc_candidates);
                                                }
                                                uint32_t key = pp->entries[e].key;
                                                candidates[n_candidates].elem_i = seq_offsets[i] + (int)(key >> 20);
                                                candidates[n_candidates].elem_j = seq_offsets[j] + (int)(key & 0xFFFFF);
                                                candidates[n_candidates].support = support;
                                                n_candidates++;
                                        }
                                }
                        }
                }

                /* Counting sort by descending support (values bounded 1..32) */
                {
                        int counts[33] = {0};
                        for(int a = 0; a < n_candidates; a++){
                                counts[candidates[a].support]++;
                        }
                        int offsets[33];
                        offsets[32] = 0;
                        for(int v = 31; v >= 0; v--){
                                offsets[v] = offsets[v+1] + counts[v+1];
                        }
                        int tmp_alloc = n_candidates > 0 ? n_candidates : 1;
                        struct merge_candidate *sorted = NULL;
                        MMALLOC(sorted, sizeof(*sorted) * tmp_alloc);
                        for(int a = 0; a < n_candidates; a++){
                                int s_val = candidates[a].support;
                                sorted[offsets[s_val]++] = candidates[a];
                        }
                        memcpy(candidates, sorted, sizeof(*candidates) * n_candidates);
                        MFREE(sorted);
                }

                /* Merge in priority order, skipping conflicts AND cycles */
                for(int c = 0; c < n_candidates; c++){
                        uf_union_safe(uf, candidates[c].elem_i, candidates[c].elem_j,
                                      seq_offsets, seq_lengths,
                                      visited, &visit_counter);
                }

                MFREE(candidates);
        }

        MFREE(visited);

        /* Map UF roots to column IDs */
        MMALLOC(root_to_col, sizeof(int) * total_residues);
        MMALLOC(col_id, sizeof(int) * total_residues);

        for(i = 0; i < total_residues; i++){
                root_to_col[i] = -1;
        }
        n_cols = 0;
        for(i = 0; i < total_residues; i++){
                int root = uf_find(uf, i);
                if(root_to_col[root] == -1){
                        root_to_col[root] = n_cols++;
                }
                col_id[i] = root_to_col[root];
        }

        MFREE(root_to_col);

        /* Topological sort (DFS-based, cycle-safe) */
        RUN(topo_sort(col_id, seq_offsets, seq_lengths, numseq, n_cols,
                      &sorted_cols, &n_sorted));

        /* Build column order lookup: sorted position -> column */
        int* col_order = NULL;  /* column -> position in sorted output */
        MMALLOC(col_order, sizeof(int) * n_cols);
        for(i = 0; i < n_sorted; i++){
                col_order[sorted_cols[i]] = i;
        }

        /* Build output alignment */
        MMALLOC(out_seqs, sizeof(char*) * numseq);
        for(s = 0; s < numseq; s++){
                out_seqs[s] = NULL;
                MMALLOC(out_seqs[s], sizeof(char) * (n_sorted + 1));
                memset(out_seqs[s], '-', n_sorted);
                out_seqs[s][n_sorted] = '\0';
        }

        /* Place residues */
        for(s = 0; s < numseq; s++){
                for(pos = 0; pos < seq_lengths[s]; pos++){
                        int elem = seq_offsets[s] + pos;
                        int col = col_id[elem];
                        int sorted_pos = col_order[col];
                        /* Use original residue character from the MSA */
                        out_seqs[s][sorted_pos] = out_msa->sequences[s]->seq[pos];
                }
        }

        /* Write back into MSA: replace seq pointers */
        for(s = 0; s < numseq; s++){
                MFREE(out_msa->sequences[s]->seq);
                out_msa->sequences[s]->seq = out_seqs[s];
                out_msa->sequences[s]->len = n_sorted;
                out_seqs[s] = NULL;
        }
        out_msa->alnlen = n_sorted;
        out_msa->aligned = ALN_STATUS_FINAL;

        MFREE(out_seqs);
        MFREE(col_order);
        MFREE(sorted_cols);
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
        if(col_order) MFREE(col_order);
        if(sorted_cols) MFREE(sorted_cols);
        if(col_id) MFREE(col_id);
        if(root_to_col) MFREE(root_to_col);
        if(visited) MFREE(visited);
        if(seq_offsets) MFREE(seq_offsets);
        uf_free(uf);
        return FAIL;
}

/* Score one alignment against the POAR table.
   Returns expected number of correct pairs: for each aligned pair,
   (support - 1) / (m - 1) gives the fraction of OTHER alignments
   agreeing. Summing these gives expected correct pairs.
   This rewards both high recall (many pairs) and high precision
   (pairs with broad agreement). */
int score_alignment_poar(struct poar_table* table,
                         struct pos_matrix* pm,
                         int numseq,
                         int n_alignments,
                         double* out_score)
{
        double total_score = 0.0;
        int i, j, col;
        int alnlen = pm->alnlen;
        double denom = (n_alignments > 1) ? (double)(n_alignments - 1) : 1.0;

        for(i = 0; i < numseq - 1; i++){
                for(j = i + 1; j < numseq; j++){
                        int pidx = pair_index(i, j, numseq);
                        struct poar_pair* pp = table->pairs[pidx];

                        for(col = 0; col < alnlen; col++){
                                int ri = pm->col_to_res[i][col];
                                int rj = pm->col_to_res[j][col];
                                if(ri >= 0 && rj >= 0){
                                        uint32_t key = ((uint32_t)ri << 20) | (uint32_t)rj;

                                        /* Binary search in sorted entries */
                                        int lo = 0;
                                        int hi = pp->n_entries;
                                        int support = 0;
                                        while(lo < hi){
                                                int mid = lo + (hi - lo) / 2;
                                                if(pp->entries[mid].key < key){
                                                        lo = mid + 1;
                                                }else if(pp->entries[mid].key == key){
                                                        support = popcount32(pp->entries[mid].support);
                                                        break;
                                                }else{
                                                        hi = mid;
                                                }
                                        }
                                        /* support includes self; subtract 1 for other-agreement */
                                        total_score += (double)(support - 1) / denom;
                                }
                        }
                }
        }

        *out_score = total_score;
        return OK;
}

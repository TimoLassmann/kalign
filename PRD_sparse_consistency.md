# PRD: Sparse Anchor Consistency Bonus Matrix

## 1. Problem Statement

The anchor consistency bonus matrix is currently stored as a dense `float[len_a * len_b]` array. This causes two critical failures for large DNA/RNA families:

**Integer overflow in DP indexing.** The DP inner loops compute `i * m->consistency_stride + j` using `int` arithmetic. When `len_a` and `len_b` both exceed ~46,340, the product overflows 32-bit `int` (46,341 × 46,341 = 2,147,488,281 > INT_MAX). The resulting negative index causes SIGBUS/segfault.

**Excessive memory.** A 50,000 × 50,000 dense matrix requires `50000 * 50000 * 4 = 10 GB` for a single merge step. With OpenMP parallelism, multiple merges may be in flight simultaneously.

**The matrix is extremely sparse.** Each anchor `k` contributes at most one entry per row `i`. With K anchors (typically 8), each row has at most K non-zero entries. For a 50k × 50k matrix: 400,000 non-zero entries out of 2.5 billion cells (0.016% density).

### Affected cases

MDSA DNA benchmark families with large profile lengths during progressive alignment:

| Family | Seqs | Max Length | N×MaxLen |
|--------|------|-----------|----------|
| RV60_sushi_ref6 | 384 | 5,178 | 1,988,352 |
| RV60_ank_ref6 | 210 | 9,297 | 1,952,370 |
| RV40_BB40047 | 54 | 13,644 | 736,776 |
| RV40_BB40023 | 22 | 23,769 | 522,918 |
| AAA (SMART) | 427 | 1,251 | 534,177 |

## 2. Proposed Solution

Replace the dense `float*` bonus matrix with a per-row sparse structure:

```c
struct sparse_bonus {
    int*   cols;    /* cols[i * K + k] = column index, or -1 if unused */
    float* vals;    /* vals[i * K + k] = bonus value */
    int    n_rows;  /* = len_a (number of DP rows) */
    int    K;       /* max entries per row (= n_anchors) */
};
```

**Memory comparison:**

| Scenario | Dense | Sparse (K=8) | Ratio |
|----------|-------|--------------|-------|
| Protein 500×500 | 1 MB | 32 KB | 32× |
| Protein 2000×2000 | 16 MB | 128 KB | 125× |
| DNA 50000×50000 | 10 GB | 3.2 MB | 3125× |

**DP lookup** via inline function scanning K=8 slots per cell:

```c
static inline float sparse_bonus_lookup(const struct sparse_bonus* sb, int i, int j)
{
    float bonus = 0.0f;
    const int base = i * sb->K;
    for(int k = 0; k < sb->K; k++){
        if(sb->cols[base + k] < 0) break;       /* early exit: unused slots at end */
        if(sb->cols[base + k] == j)
            bonus += sb->vals[base + k];
    }
    return bonus;
}
```

## 3. Detailed Design

### 3.1 New struct and functions in `anchor_consistency.h`

Add after `struct consistency_table`:

```c
struct sparse_bonus {
    int*   cols;
    float* vals;
    int    n_rows;
    int    K;
};

static inline float sparse_bonus_lookup(const struct sparse_bonus* sb, int i, int j)
{
    float bonus = 0.0f;
    const int base = i * sb->K;
    for(int k = 0; k < sb->K; k++){
        if(sb->cols[base + k] < 0) break;
        if(sb->cols[base + k] == j)
            bonus += sb->vals[base + k];
    }
    return bonus;
}

EXTERN void sparse_bonus_free(struct sparse_bonus* sb);
```

Update function signatures:

```c
EXTERN int anchor_consistency_get_bonus(struct consistency_table* ct,
                                        int seq_a, int len_a,
                                        int seq_b, int len_b,
                                        struct sparse_bonus** bonus_out);

EXTERN int anchor_consistency_get_bonus_profile(struct consistency_table* ct,
                                                struct msa* msa,
                                                int node_a, int len_a,
                                                int node_b, int len_b,
                                                struct sparse_bonus** bonus_out);
```

### 3.2 Changes to `aln_struct.h`

Add forward declaration:
```c
struct sparse_bonus;
```

**Remove:**
```c
float* consistency;
int consistency_stride;
```

**Replace with:**
```c
struct sparse_bonus* consistency;
```

### 3.3 Changes to `aln_mem.c`

In `alloc_aln_mem`, change:
```c
m->consistency = NULL;
m->consistency_stride = 0;
```
To:
```c
m->consistency = NULL;
```

### 3.4 Rewrite bonus construction in `anchor_consistency.c`

#### `sparse_bonus_free`:
```c
void sparse_bonus_free(struct sparse_bonus* sb)
{
    if(sb){
        if(sb->cols) MFREE(sb->cols);
        if(sb->vals) MFREE(sb->vals);
        MFREE(sb);
    }
}
```

#### `anchor_consistency_get_bonus_profile` (and `_get_bonus`):

Replace the dense allocation:
```c
MMALLOC(bonus, sizeof(float) * len_a * len_b);
memset(bonus, 0, sizeof(float) * len_a * len_b);
```

With sparse allocation:
```c
struct sparse_bonus* sb = NULL;
MMALLOC(sb, sizeof(struct sparse_bonus));
sb->cols = NULL;
sb->vals = NULL;
sb->n_rows = len_a;
sb->K = K;

MMALLOC(sb->cols, sizeof(int) * len_a * K);
MMALLOC(sb->vals, sizeof(float) * len_a * K);
for(i = 0; i < len_a * K; i++){
    sb->cols[i] = -1;
    sb->vals[i] = 0.0f;
}
```

Replace the dense accumulation (current line 532):
```c
bonus[i * len_b + bj] += per_anchor_weight * conf_a[i] * inv_conf_b[ak_pos];
```

With sparse slot insertion:
```c
{
    float val = per_anchor_weight * conf_a[i] * inv_conf_b[ak_pos];
    int base = i * K;
    int slot = -1;
    for(int s = 0; s < K; s++){
        if(sb->cols[base + s] == bj){ slot = s; break; }  /* existing entry */
        if(sb->cols[base + s] < 0){ slot = s; break; }    /* empty slot */
    }
    if(slot >= 0){
        sb->vals[base + slot] += val;
        sb->cols[base + slot] = bj;
    }
}
```

Note: the slot search handles both cases — accumulating into an existing entry (two anchors map the same row to the same column) and claiming a new slot. The `cols = bj` write is idempotent for existing entries.

Return `*bonus_out = sb` instead of `*bonus_out = bonus`.

### 3.5 Changes to DP consumer files (12 access sites)

All 12 sites follow the same transformation. Current:
```c
if(m->consistency){
    pa += m->consistency[i * m->consistency_stride + j];
}
```

Becomes:
```c
if(m->consistency){
    pa += sparse_bonus_lookup(m->consistency, ROW, COL);
}
```

Where ROW and COL are the same expressions currently used for `i` and `j`.

**`aln_seqseq.c`** — 4 sites:
- Forward inner loop (line 83): ROW=`i`, COL=`j`
- Forward end-of-row (line 105): ROW=`i`, COL=`j`
- Backward inner loop (line 199): ROW=`starta + i`, COL=`j`
- Backward end-of-row (line 223): ROW=`starta + i`, COL=`j`

**`aln_seqprofile.c`** — 4 sites:
- Forward inner loop (line 82): ROW=`i`, COL=`j`
- Forward end-of-row (line 107): ROW=`i`, COL=`j`
- Backward inner loop (line 193): ROW=`m->starta_2 + i`, COL=`j`
- Backward end-of-row (line 216): ROW=`m->starta_2 + i`, COL=`j`

**`aln_profileprofile.c`** — 4 sites:
- Forward inner loop (line 108): ROW=`i`, COL=`j`
- Forward end-of-row (line 138): ROW=`i`, COL=`j`
- Backward inner loop (line 251): ROW=`m->starta_2 + i`, COL=`j`
- Backward end-of-row (line 280): ROW=`m->starta_2 + i`, COL=`j`

### 3.6 Changes to `aln_run.c`

**`do_align`** (2 locations):

Setup (lines 259-261, 290-294): Remove `m->consistency_stride = dp_cols;` — stride is embedded in the sparse struct.

Teardown (lines 401-405): Replace `MFREE(m->consistency)` with `sparse_bonus_free(m->consistency)`.

**`do_align_inline_refine`** (2 locations):

Setup (lines 566-567, 595-598): Same — remove stride assignment.

Teardown (lines 736-740): Same — use `sparse_bonus_free`.

## 4. Implementation Ordering

Steps 1-2 are additive (nothing breaks). Steps 3-7 must be atomic (single commit).

1. Add `struct sparse_bonus`, `sparse_bonus_lookup`, `sparse_bonus_free` to header/source
2. (Can be tested in isolation with a unit test)
3. Change `aln_struct.h` — replace `float* consistency` + `int consistency_stride` with `struct sparse_bonus*`
4. Update `aln_mem.c` — remove `consistency_stride` init
5. Update 12 DP access sites to use `sparse_bonus_lookup`
6. Rewrite `anchor_consistency_get_bonus` and `_get_bonus_profile` to produce sparse output
7. Update `aln_run.c` setup/teardown

## 5. Testing Strategy

### 5.1 Bit-exact regression (CRITICAL)

For all BAliBASE protein families where the current implementation works correctly:
- Run alignment with current (dense) implementation, save output
- Run alignment with sparse implementation, compare output
- **Must be byte-identical** — same bonus values → same DP scores → same alignment paths

```bash
# Save reference outputs with current code
for f in tests/data/*.tfa; do
    ./build/src/kalign -i "$f" -o "/tmp/dense_$(basename $f)" --consistency 8
done

# After sparse implementation, compare
for f in tests/data/*.tfa; do
    ./build/src/kalign -i "$f" -o "/tmp/sparse_$(basename $f)" --consistency 8
    diff "/tmp/dense_$(basename $f)" "/tmp/sparse_$(basename $f)"
done
```

Also run the full BAliBASE benchmark and compare F1/TC scores:
```bash
uv run python -m benchmarks.runner --dataset balibase
```

### 5.2 Large DNA families (new capability)

Verify that MDSA families that crashed with dense implementation now complete:
```bash
./build-asan/src/kalign -i benchmarks/data/downloads/mdsa/unaligned/smart/AAA.afa \
    -o /tmp/test.fa --consistency 8 --realign 2
```

Run the full MDSA DNA benchmark:
```bash
uv run python -m benchmarks.optimize_unified --dataset mdsa --pop-size 20 --n-gen 2 --n-workers 4
```

### 5.3 Existing test suite

All existing tests must pass unchanged:
```bash
cd build && make test                        # C tests
uv run pytest tests/python/ -v               # Python tests
```

### 5.4 Performance comparison

Time the BAliBASE benchmark with both implementations:
```bash
# Before (dense)
time uv run python -c "
from benchmarks.scoring import run_case
from benchmarks.datasets import balibase_cases
for c in balibase_cases():
    run_case(c, method='python_api', refine='none')
"

# After (sparse) — same command
```

Expectation: equal or faster (better cache behavior). Any slowdown > 5% warrants investigation.

### 5.5 Unit test for sparse_bonus

Add a C test that:
1. Creates a `sparse_bonus` with known values
2. Verifies `sparse_bonus_lookup` returns correct values for filled slots
3. Returns 0.0f for empty slots and out-of-band columns
4. Correctly accumulates when two anchors map to the same (row, col)
5. Handles `sb == NULL` (returns 0.0f)

## 6. Risk Analysis

### 6.1 Incorrect accumulation
**Risk**: Two anchors mapping the same row to the same column must accumulate (add), not overwrite.
**Mitigation**: Slot-finding logic checks `cols[slot] == bj` before writing. The `+=` on vals handles accumulation. Unit test covers this case.

### 6.2 Hirschberg sub-problem offsets
**Risk**: The divide-and-conquer DP uses `starta`, `starta_2` offsets for row indices.
**Mitigation**: The sparse lookup receives the same (row, col) coordinates as the dense indexing. No change in semantics. Bit-exact regression testing confirms correctness.

### 6.3 Thread safety
**Risk**: Each merge step's `aln_mem` gets its own `sparse_bonus*`. The struct is read-only during DP.
**Mitigation**: Same thread-safety model as current dense implementation. No concurrent writes.

### 6.4 Rollback strategy
The change is purely internal — no API, file format, or CLI changes. If any regression is found, revert the single commit to restore the dense implementation.

## 7. Future Considerations

- **Row-pointer hoisting**: Hoist `cols + i*K` and `vals + i*K` outside the j-loop for better codegen
- **SIMD lookup**: For K=8, a single AVX2 `_mm256_cmpeq_epi32` could find the matching slot in one instruction
- **Adaptive K**: If future anchor schemes use K > 16, switch to sorted columns with binary search

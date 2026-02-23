# PRD: Gap-Fraction Down-Weighting & Normalized Gap Embedding

## Problem

Kalign's alignment quality degrades for divergent sequences. At BAliBASE tree
depth 4.0, kalign achieves nRF ≈ 0.30 vs MAFFT's 0.16. Analysis of the
profile-profile scoring model reveals two architectural issues that compound at
high divergence:

### Issue 1: No column quality weighting

When aligning two profiles, every column contributes equally to the match score
regardless of how gappy it is. A column with 1 real residue and 99 gaps
contributes the same structural weight as a column with 100 residues. MUSCLE
addresses this with `(1-fgap)` weighting; MAFFT uses similar position-specific
gap logic. Kalign has no equivalent.

At high divergence, many columns accumulate gaps from early progressive
alignment errors. Without down-weighting, these noisy columns have equal
influence on later alignment steps, propagating errors up the guide tree.

### Issue 2: Unbounded gap embedding in profiles

In `update_n()`, when a gap is placed at a column, kalign subtracts
`gpe × nsip_gapped` from all 23 substitution score slots. This "gap embedding"
is an absolute penalty that grows linearly with profile size. For a column
gapped 50 times, the substitution scores are reduced by `50 × gpe` per slot.

At high divergence with large profiles, this creates "score craters" at gappy
positions — the substitution scores become very negative, making these columns
repulsive to match against in subsequent alignments. The effect is self-
reinforcing: once gaps accumulate, more gaps are attracted.

## Solution

Two complementary changes, each controlled independently:

### Feature A: (1-fgap) Column Down-Weighting

Multiply the match score at each DP cell by `(1 - fgap₁) × (1 - fgap₂)`,
where `fgap = 1 - n_residues / n_total_seqs`.

- A clean column (fgap=0) contributes full weight
- A 50%-gapped column pair contributes 25% weight
- A 90%-gapped column pair contributes 1% weight

This is MUSCLE's "log-expectation" correction.

**Implementation:** Store residue count in profile slot 26 (currently unused).
At DP time, compute `wt = (prof[26] / nsip) × ...` and multiply the dot-product
match score.

### Feature B: Normalized Gap Embedding

Instead of subtracting `gpe × nsip_gapped` from substitution scores in
`update_n()`, subtract `gpe × fgap` where `fgap = nsip_gapped / nsip_total`.

- Bounded by `gpe` regardless of profile size
- Gap embedding becomes a fraction (0..1) of gpe, not an unbounded integer multiple
- Substitution scores at gappy columns stay in a reasonable range

**Implementation:** In `update_n()`, divide the gap penalty by `(sipa + sipb)`
before subtracting from substitution scores.

## Parameters

Add to `struct aln_param`:
```c
int use_gap_weight;         /* 1 = (1-fgap) column weighting */
int normalize_embedding;    /* 1 = normalize gap embed by profile size */
```

Default: both ON for protein and DNA. RNA: both OFF (conservative).

Python API: `gap_weight=True/False`, `normalize_embedding=True/False`.

## Files to Modify

### 1. `lib/src/aln_param.h` — Add fields

```c
int use_gap_weight;         /* 1 = (1-fgap) column weighting */
int normalize_embedding;    /* 1 = normalize gap embed by profile size */
```

### 2. `lib/src/aln_param.c` — Set defaults

In `aln_param_init()`, after biotype selection:
- Protein/DNA: `use_gap_weight = 1; normalize_embedding = 1;`
- RNA: `use_gap_weight = 0; normalize_embedding = 0;`

### 3. `lib/src/aln_struct.h` — Add nsip fields to `struct aln_mem`

```c
int nsip_a;   /* total sequences in profile a */
int nsip_b;   /* total sequences in profile b */
```

### 4. `lib/src/aln_setup.c` — Profile slot 26 + normalized embedding

**`make_profile_n()`:** Set `prof[26] = 1.0` at each sequence position, `0.0`
at sentinels.

**`update_n()`:** When `normalize_embedding`:
- Replace `gp = ap->gpe * (float)sipa` with
  `gp = ap->gpe * (float)sipa / (float)(sipa + sipb)`
- Same for gpo and tgpe terms

Slot 26 needs no explicit handling — element-wise sum/copy already works:
- Match: `newp[26] = profa[26] + profb[26]` = total residue count
- Gap-in-A: `newp[26] = profb[26]` = B's residue count only
- Gap-in-B: `newp[26] = profa[26]` = A's residue count only

### 5. `lib/src/aln_run.c` — Set nsip_a/nsip_b

In `do_align()`, before calling `aln_runner()`:
```c
m->nsip_a = msa->nsip[a];
m->nsip_b = msa->nsip[b];
```

For seq-seq (nsip==1 on both sides), nsip_a=1, nsip_b=1 → fgap=0 → no effect.

### 6. `lib/src/aln_profileprofile.c` — Apply (1-fgap) weighting

In forward, backward, and meetup functions, after computing the dot-product
match score, multiply by `wt`:

```c
// After dot-product loop:
if(m->ap->use_gap_weight){
    float nres1 = prof1[26];
    float nres2 = prof2[26];  // prof2 at current position
    float wt = (nres1 * nres2) / ((float)m->nsip_a * (float)m->nsip_b);
    // Apply to the dot product only, not to the transition scores
}
```

Three functions × inner loop + boundary = ~6 sites.

### 7. `lib/src/aln_seqprofile.c` — Apply (1-fgap) weighting

For seq-profile alignment, the sequence side has fgap=0 (always 1 residue).
Only the profile side contributes:

```c
if(m->ap->use_gap_weight){
    float wt = prof1[26] / (float)m->nsip_a;
    // multiply match score by wt
}
```

Three functions × inner loop + boundary = ~6 sites.

### 8. `lib/src/aln_seqseq.c` — No DP changes needed

Seq-seq always has fgap=0 on both sides. No modification required.

### 9. Python bindings

**`python-kalign/__init__.py`:** Add `gap_weight: Optional[bool] = None` and
`normalize_embedding: Optional[bool] = None` to `align()`,
`align_from_file()`, `align_file_to_file()`.

**`python-kalign/_core.cpp`:** Pass through to C API. Sentinel: `None` → -1
(use default), `True` → 1, `False` → 0.

### 10. Benchmark script

**`benchmarks/gap_weight_bench.py`:** Sweep 4 configurations on BAliBASE:
1. Baseline (both off)
2. (1-fgap) only
3. Normalized embedding only
4. Both on

Report F1, recall, precision per config + per-case deltas.

## Verification

1. **C tests:** `cd build && cmake .. && make && make test` — 13/13 pass
2. **Python tests:** `uv run pytest tests/python/ -v` — 140/140 pass
3. **BAliBASE benchmark:**
   ```
   uv run python benchmarks/gap_weight_bench.py -j 12
   ```
   Baseline F1 ≈ 0.719. Target: improvement with gap weighting.
4. **Phylo benchmark in container** (if BAliBASE shows improvement):
   Rebuild container, run `--depths 4.0`. Target: reduced nRF gap to MAFFT.

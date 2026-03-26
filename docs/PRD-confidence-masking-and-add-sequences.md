# PRD: Confidence masking + Add sequences to existing alignment

## Feature 1: Confidence-filtered output

### What

Expose the per-column ensemble confidence scores to the user. Allow masking (lowercasing or gap-replacing) columns below a threshold. Only meaningful for ensemble modes (recall, accurate) which compute POAR confidence.

### Why

No other fast aligner provides per-column reliability. GUIDANCE2 does but is extremely slow (runs external aligners repeatedly). Kalign's ensemble already does the multiple runs — confidence is free. This turns kalign into a one-stop alignment + quality filter tool, directly useful for phylogenetics, positive selection, and structure prediction pipelines.

### Paper angle

**New figure or table:** Score BAliBASE alignments on ONLY the columns kalign considers confident (confidence >= threshold). Show that SP/TC on confident columns is dramatically higher than on all columns — proving the confidence scores are meaningful. Compare at several thresholds (0.3, 0.5, 0.7, 0.9).

This directly addresses the TC weakness: "kalign's overall TC is lower than competitors, but if you trust its confidence scores and filter, the remaining columns have TC comparable to or better than MAFFT/MUSCLE."

**Second comparison:** Run phylogenetic tree inference (IQ-TREE) on:
1. Full kalign accurate alignment
2. Kalign accurate filtered at confidence >= 0.5
3. MAFFT alignment
4. MUSCLE alignment

Compare Robinson-Foulds distance to the true tree (from INDELible simulations, which the manuscript already has). If filtered kalign gives better trees, that's a compelling result.

### C CLI interface

```
kalign -i seqs.fa -o aligned.fa --mode accurate --confidence-threshold 0.7

Options:
  --confidence-threshold FLOAT   Mask columns with confidence below this value.
                                  0.0 = no masking (default). Requires ensemble mode.
  --confidence-style STRING      "lowercase" (default) or "remove".
                                  lowercase: uncertain residues become lowercase
                                  remove: uncertain columns replaced with gaps
  --confidence-output FILE       Write per-column confidence values to a separate file
                                  (one float per line, one line per column)
```

The `--confidence-output` flag is useful for downstream tools that want the raw scores (e.g., custom trimming scripts).

### Python API

```python
# align() and align_from_file() already return confidence when using ensemble modes
result = kalign.align(sequences, mode="accurate")
# result.confidence is a list of per-column floats [0.0, 1.0]

# New: mask_alignment convenience function
masked = kalign.mask_alignment(result, threshold=0.7, style="lowercase")
# masked.sequences has lowercase residues in low-confidence columns

# New: filter_alignment removes low-confidence columns entirely
filtered = kalign.filter_alignment(result, threshold=0.7)
# filtered.sequences are shorter — only high-confidence columns remain

# align_file_to_file gains optional confidence args
kalign.align_file_to_file(
    "in.fa", "out.fa",
    mode="accurate",
    confidence_threshold=0.7,
    confidence_style="lowercase",  # or "remove"
)

# Write raw confidence scores
kalign.write_confidence("confidence.txt", result)
```

### Implementation

**C library changes:**

1. `lib/src/msa_io.c` — modify `write_fasta` / `write_clustal` / `write_msf` to apply masking:
   - Accept threshold + style parameters
   - Before writing each character: if `col_confidence[col] < threshold`, apply style
   - For "lowercase": `seq[col] = tolower(seq[col])` (only if not gap)
   - For "remove": `seq[col] = '-'`

2. `src/run_kalign.c` — add CLI flags:
   - `--confidence-threshold` → `float conf_threshold`
   - `--confidence-style` → `enum {CONF_LOWERCASE, CONF_REMOVE}`
   - `--confidence-output` → write `msa->col_confidence[]` to file
   - Validate: if threshold > 0 and mode is not ensemble, warn and ignore

3. `python-kalign/__init__.py` — add `mask_alignment()`, `filter_alignment()`, `write_confidence()`
4. `python-kalign/_core.cpp` — pass threshold/style to `kalign_write_msa` or post-process in Python

**Key constraint:** Confidence is only available after ensemble alignment (`col_confidence` is NULL for single-run modes). The CLI and Python API must handle this gracefully — warn if the user requests masking with fast/default mode.

### Testing

1. **Unit test:** Align BB11001 with accurate mode, verify `col_confidence` is populated, apply threshold=0.5, check that output has lowercase characters in the right columns

2. **Round-trip test:** Write masked alignment, read it back, verify non-masked residues are unchanged

3. **Quality test (for paper):** Score all 218 BAliBASE cases:
   - Compute SP and TC on ALL columns (standard)
   - Compute SP and TC on ONLY columns with confidence >= threshold
   - Show improvement as threshold increases
   - Script: `scripts/bench_confidence_filtering.py`

4. **Phylogenetic test (for paper):** Use the manuscript's INDELible simulations:
   - Align with kalign accurate
   - Filter at various thresholds
   - Run IQ-TREE
   - Compare RF distance to true tree
   - Already partially done in the manuscript pipeline

---

## Feature 2: Add sequences to existing alignment

### What

Align new sequences against an existing (fixed) alignment without modifying the existing sequences. Each new sequence is independently aligned to the consensus profile of the existing alignment.

### Why

This is one of the most requested features in MSA tools. MAFFT `--add` is heavily cited specifically for this. Use cases:
- Metagenomics: add new sample sequences to a reference alignment
- Phylogenetics: add new taxa to a growing tree alignment
- Viral surveillance: daily additions to reference alignments (SARS-CoV-2, influenza)
- Database maintenance: adding sequences to curated family alignments

### Paper angle

**Benchmark against MAFFT --add:**

1. Take BAliBASE reference alignments. For each case:
   - Hold out 20% of sequences as "new"
   - Use remaining 80% as the "existing alignment"
   - Add the held-out sequences with kalign --add and mafft --add
   - Score the added sequences against the full reference alignment
   - Measure time

2. Larger-scale test with simulated data:
   - INDELible simulation with 500 sequences
   - Use first 400 as reference alignment (align normally)
   - Add remaining 100 with --add
   - Compare to full 500-sequence alignment
   - Time comparison: kalign --add vs mafft --add

This gives both quality and speed comparisons. Kalign should be faster (SIMD Hirschberg). Quality depends on how well the seq-to-profile alignment places gaps.

### C CLI interface

```
kalign --add new_seqs.fa --existing aligned.fa -o combined.fa

Options:
  --add FILE           Unaligned sequences to add to an existing alignment
  --existing FILE      Existing alignment (FASTA/MSF/Clustal). These sequences
                       are NOT re-aligned — their gaps are preserved exactly.
  -o FILE              Output: existing sequences (unchanged) + new sequences
                       (with gaps inserted to fit the existing column structure)
  --nthreads N         Parallel seq-to-profile alignments for each new sequence
```

**Behavior:**
- Read existing alignment → build profile
- For each new sequence: align to profile, insert gaps
- Output: existing sequences verbatim + new sequences with gaps
- Existing sequences are NEVER modified
- Column count may increase if new sequences have insertions not present in the existing alignment (new gap columns inserted in ALL sequences at those positions)

### Python API

```python
# From files
kalign.add_to_alignment(
    existing="reference.fa",
    new_sequences="new.fa",
    output="combined.fa",
    n_threads=8,
)

# In-memory
existing = kalign.align(ref_sequences, mode="accurate")
combined = kalign.add_sequences(existing, new_sequences, n_threads=8)
# combined.sequences = existing (unchanged) + new (gapped)
# combined.names = existing names + new names
```

### Implementation

**C library:**

1. **New public API function in `kalign.h`:**
   ```c
   int kalign_add_sequences(struct msa* existing_aln,
                            struct msa* new_seqs,
                            int n_threads);
   ```
   - `existing_aln`: finalized alignment (sequences have gap chars, alnlen set)
   - `new_seqs`: unaligned sequences
   - After call: `existing_aln` contains original + new sequences, all aligned
   - Returns OK/FAIL

2. **Core algorithm in new file `lib/src/aln_add.c`:**

   ```
   kalign_add_sequences(existing, new_seqs, n_threads):
     1. Detect biotype, encode new sequences to internal alphabet
     2. Build consensus profile from existing alignment
        - Walk all columns of existing alignment
        - At each column: count residue frequencies (weighted by sequence count)
        - Store as float[64] per column (same format as progressive profiles)
     3. For each new sequence (parallelizable with tp_parallel_for):
        a. Run seq-to-profile Hirschberg alignment (aln_seqprofile)
        b. Extract gap positions from alignment path
        c. Insert gaps into new sequence to match existing column structure
        d. If new sequence has insertions not in existing alignment:
           - Record insertion positions and lengths
     4. If any insertions were found:
        - Insert new gap columns into ALL sequences (existing + new) at those positions
        - Update alnlen
     5. Append new sequences to existing MSA
   ```

3. **Profile building from existing alignment:**
   - Similar to what `make_profile_n` does during progressive alignment
   - But operates on finalized character sequences with gap chars
   - Convert back to internal representation for the DP
   - Or: build profile directly from character frequencies per column

4. **CLI in `src/run_kalign.c`:**
   - Parse `--add` and `--existing` flags
   - Read both files
   - Call `kalign_add_sequences`
   - Write combined output

5. **Python bindings in `_core.cpp`:**
   - New `add_sequences` function
   - Reads existing alignment file, reads new sequences file, calls C API

**Parallelism:** Each new sequence's alignment to the profile is independent → `tp_parallel_for` over new sequences. This is embarrassingly parallel and should scale well.

**Key design decision — handling insertions in new sequences:**

When a new sequence has residues that don't map to any existing column, we must insert new columns. This affects ALL sequences (existing ones get gaps at those positions). Two approaches:

- **Strict (MAFFT --add style):** Never insert new columns. New sequence insertions are forced into existing columns or dropped. Existing alignment column count is preserved exactly.
- **Flexible:** Insert new columns as needed. More accurate for the new sequences but modifies the column structure.

Recommend: **strict mode as default** (existing sequences completely untouched, even column count preserved), with a `--add-insertions` flag for flexible mode. The strict mode is what users expect from "add to existing alignment."

In strict mode, the seq-to-profile alignment path may indicate insertions in the new sequence. These residues are simply lowercase or dropped (user choice). This matches MAFFT --add behavior.

### Testing

1. **Identity test:** Add a sequence that's already in the alignment. Result should be identical to the original.

2. **Residue preservation test:** After adding, count non-gap characters in each new sequence — must equal original sequence length.

3. **Existing-unchanged test:** After adding, the existing sequences must be byte-identical to the input existing alignment.

4. **BAliBASE holdout test (for paper):**
   - For each BAliBASE case: hold out 20% sequences
   - Align 80% normally → existing alignment
   - Add held-out 20% with kalign --add
   - Score the full result against BAliBASE reference
   - Compare to: MAFFT --add, and full kalign alignment of all sequences

5. **Speed test (for paper):**
   - 1000-sequence reference alignment + add 100 new sequences
   - Time kalign --add vs mafft --add
   - Repeat at 5000 + 500, 10000 + 1000

6. **Simulation test (for paper):**
   - INDELible: generate true alignment of 500 sequences
   - Split: 400 reference + 100 to add
   - Align 400 with each tool
   - Add 100 with each tool's --add
   - Score against true alignment
   - Also run full 500-sequence alignment as upper bound

---

## Implementation order

```
Phase 1: Confidence masking (simpler, builds on existing infrastructure)
  1a. C library: masking function in msa_op.c
  1b. CLI: --confidence-threshold, --confidence-style, --confidence-output
  1c. Python: mask_alignment(), filter_alignment(), write_confidence()
  1d. Tests: unit + BAliBASE quality at thresholds
  1e. Paper benchmark: SP/TC on confident columns only

Phase 2: Add sequences (more complex, new alignment mode)
  2a. C library: kalign_add_sequences in aln_add.c
  2b. Profile building from existing alignment
  2c. Seq-to-profile alignment + gap insertion
  2d. CLI: --add, --existing
  2e. Python: add_to_alignment(), add_sequences()
  2f. Tests: identity + preservation + unchanged + holdout
  2g. Paper benchmark: BAliBASE holdout + speed vs MAFFT --add
```

## Estimated complexity

| Component | Lines of C | Difficulty |
|-----------|:----------:|:----------:|
| Confidence masking (msa_op.c) | ~50 | Easy |
| Confidence CLI flags | ~30 | Easy |
| Confidence Python API | ~80 | Easy |
| Confidence paper benchmark script | ~150 | Easy |
| Add-sequences core (aln_add.c) | ~300 | Medium |
| Add-sequences profile builder | ~150 | Medium |
| Add-sequences CLI | ~50 | Easy |
| Add-sequences Python API | ~100 | Easy |
| Add-sequences paper benchmark | ~200 | Medium |

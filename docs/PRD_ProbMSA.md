# ProbMSA: Probabilistic Multiple Sequence Alignment via Quality-Annotated Pair-HMM Profiles

## Product Requirements Document (PRD)

**Version:** 1.0
**Date:** February 2026

---

## 1. Executive Summary

ProbMSA is a progressive multiple sequence alignment algorithm that achieves the accuracy benefits of probabilistic/consistency-based methods (ProbCons, MSAProbs) at the speed of progressive aligners (FAMSA, MAFFT, Kalign). The central innovation is representing alignment profiles not as flat frequency vectors, but as **quality-annotated pseudo-sequences** — where each column carries a probability distribution over residues plus a confidence/quality score derived from alignment posteriors. Every progressive merge step is a **pair-HMM forward-backward alignment** that produces posterior probabilities, which then propagate as quality annotations into the merged profile. This eliminates the need for hand-tuned gap penalties: gap open/extend costs are transition probabilities in the pair-HMM, estimable from data via Baum-Welch EM.

The implementation extends an existing, fully functional pair-HMM aligner (EA) that already supports quality-aware emissions, forward-backward computation, posterior decoding, MEA alignment, and bit-parallel pre-filtering.

---

## 2. Motivation and Literature Context

### 2.1 The Gap Penalty Problem

Substitution matrices (BLOSUM, PAM, MIQS) are derived from principled probabilistic models of sequence evolution — log-odds scores from observed substitution frequencies in aligned homologous sequences. Gap penalties, by contrast, have no analogous principled derivation. The standard affine gap model (gap-open `g_o`, gap-extend `g_e`) has two free parameters that are chosen by:

- Empirical grid search on structural benchmarks (BAliBASE, SABmark, PREFAB)
- Heuristic rules of thumb (e.g., gap-open ≈ 10–12, gap-extend ≈ 1)
- Context-free application regardless of protein family, evolutionary distance, or structural context

This is unsatisfying because the "correct" gap penalty depends on:

- Structural context (loops tolerate indels; secondary structure elements do not)
- Evolutionary distance (divergent families need different penalties than close homologs)
- Position in the progressive alignment (late merges affect more sequences)
- Alignment confidence (a gap adjacent to a confidently aligned region should be penalized differently than one near an uncertain region)

### 2.2 Existing Approaches to Gap Penalties

#### 2.2.1 MUSCLE Position-Specific Gap Correction (Edgar, 2004)

MUSCLE recognizes that inserting a gap column into a profile doesn't open a new gap in every sequence — some sequences already have gaps at that position. The correction scales gap-open penalty at each position by `(1 - f_gap)`, the fraction of sequences without an existing gap. This is motivated by correct SP score optimization.

**Limitation:** Binary — a gap is either present or not. No notion of gap confidence or alignment quality. All existing gaps are treated equivalently regardless of how they were introduced.

#### 2.2.2 FAMSA Gap Cost Innovations (Deorowicz et al., 2016)

FAMSA (from Deorowicz, Debudaj-Grabysz, and Gudyś at Silesian University of Technology) extends MUSCLE's approach in two ways:

1. **Six-case gap type correction:** Elaborate bookkeeping to determine the correct gap type (terminal open/extend, internal open/extend) when inserting columns during profile alignment, accounting for the context of neighboring gaps in each sequence (cases S1–S6).

2. **Set-size scaling:** All gap costs are multiplied by `f(k) = g_l / (1 + log2(k / g_d))` where `k` is the number of input sequences and `g_l = 45`, `g_d = 7` are experimentally tuned constants. This counteracts progressive alignment drift that widens alignments of large families.

**Limitation:** The set-size scaling is global (same factor everywhere) and the formula is a two-parameter hand-tuned heuristic. The gap type correction is elaborate bookkeeping rather than a principled model.

FAMSA was preceded by Kalign-LCS (Deorowicz et al., 2014), a variant of Kalign2 replacing Wu-Manber string matching with longest common subsequence (LCS) for distance estimation. FAMSA inherits this, adds SIMD-accelerated bit-parallel LCS, single-linkage guide trees via SLINK, in-place profile alignment with a novel gapped sequence representation, MIQS substitution matrix, and iterative refinement for k ≤ 1000 sequences. FAMSA2 (2025, bioRxiv) extends this further with stochastic medoid trees for O(MN log N) guide tree construction, enabling alignment of 3M sequences in minutes.

#### 2.2.3 Probabilistic Indel Models (TKF91/92, PRANK, BAli-Phy)

The TKF91 and TKF92 models (Thorne, Kishino, Felsenstein) treat indels as birth-death processes with explicit rate parameters, putting gaps on the same probabilistic footing as substitutions. Tools like PRANK, StatAlign, and BAli-Phy implement phylogeny-aware alignment where indel parameters are part of a full probabilistic model, estimable via maximum likelihood or MCMC.

**Limitation:** Computationally expensive (hours to days for moderate datasets). Not scalable to thousands of sequences. Realistic indel length distributions (power-law) are hard to capture tractably.

#### 2.2.4 Consistency-Based Methods (ProbCons, MSAProbs)

ProbCons computes pairwise pair-HMM forward-backward posteriors, then uses a "consistency transformation" to propagate information across all pairs before progressive alignment. MSAProbs refines this with partition function posteriors. Both achieve high accuracy but at O(N²L²) cost for N sequences of length L.

**Limitation:** The consistency transformation is a separate preprocessing step. Posterior information is not carried through the progressive alignment as an intrinsic property of the profiles.

#### 2.2.5 Iterative Parameter Learning (EM-Style Bootstrapping)

The idea of starting with naive parameters, aligning, estimating better parameters from the alignment, and re-aligning is conceptually an EM loop. The main challenge is degeneracy: naive EM can collapse to no-gaps or run away to all-gaps. Stabilization requires regularization via structural priors, transitivity consistency, or column entropy constraints. BAli-Phy's MCMC approach is the theoretically correct solution but is computationally prohibitive at scale.

#### 2.2.6 Machine Learning / Protein Language Model Approaches

Deep learning methods (DeepMSA, ESM-based aligners) sidestep explicit gap penalties by learning alignment features end-to-end. The "cost" of a gap is implicitly captured in learned representations. These are powerful but opaque and require training data.

### 2.3 The Key Insight: Quality-Annotated Profiles

None of the existing approaches propagate alignment uncertainty through the progressive alignment as an intrinsic property of the profile representation. ProbMSA fills this gap:

- **Alignment posteriors become quality annotations** that affect subsequent alignment steps
- **Gap penalties dissolve into pair-HMM transition probabilities** — no separate parameterization needed
- **MUSCLE/FAMSA gap corrections emerge naturally** from the probabilistic model
- **MEA alignment comes free** at every progressive step
- **Computational cost is identical** to standard progressive alignment (O(L²) per merge)

---

## 3. Algorithm Overview

### 3.1 High-Level Pipeline

```
Input: Set of N sequences (DNA or protein)

Phase 1: Pairwise Distances
  - For each pair (i,j): compute edit distance via bit-parallel Myers (BPM)
  - Use SIMD-accelerated BPM (AVX2 for patterns ≤ 255, serial for ≤ 63)
  - Alternative: LCS-based distance à la FAMSA/Kalign-LCS

Phase 2: Guide Tree Construction
  - Build UPGMA or single-linkage tree from distance matrix
  - For very large sets (N > 10000): use mBed embedding or medoid tree heuristic

Phase 3: Progressive Alignment with Quality Propagation
  For each internal node of the tree (bottom-up):
    a. Retrieve two child profiles: P_left, P_right
       (each position: sparse probability vector + quality/concentration parameter)
    b. Run pair-HMM forward-backward with quality-aware emissions
    c. Compute MEA alignment from posterior match probabilities
    d. Merge aligned positions into new profile:
       - Combine probability vectors (Dirichlet posterior update)
       - Compute new quality from input qualities + alignment posterior
    e. Gap columns introduced during merge receive low quality

Phase 4: Parameter Re-estimation (Optional, Iterative)
  - Re-estimate pair-HMM transition parameters (δ, ε) from completed alignment via Baum-Welch
  - Re-align from Phase 3 with updated parameters
  - Repeat until convergence (typically 2–3 iterations)

Phase 5: Iterative Refinement (Optional)
  - Split alignment at random gappy columns, re-align sub-profiles
  - Accept if SP score or column entropy improves

Output: Multiple sequence alignment in FASTA/Stockholm format
```

### 3.2 Profile Representation: The Quality-Annotated Pseudo-Sequence

Each profile position is represented as:

```c
typedef struct {
    float prob[ALPHA_SIZE];    // Probability distribution over alphabet
                               // DNA: ALPHA_SIZE = 4 (A,C,G,T)
                               // Protein: ALPHA_SIZE = 20 (amino acids)
    float quality;             // Concentration/confidence parameter (0.0 = unknown, high = confident)
    uint32_t depth;            // Number of sequences contributing to this column
    float gap_fraction;        // Fraction of sequences with a gap at this position
} ProfileColumn;

typedef struct {
    ProfileColumn *columns;
    uint32_t length;           // Number of columns
    uint32_t n_sequences;      // Number of sequences represented
} Profile;
```

**For a single sequence (leaf node):**
- `prob[base] = 1.0` for the observed residue, `0.0` for others
- `quality` = Phred-like score (high for confident positions)
  - DNA with base qualities: use actual Phred scores
  - DNA/protein without qualities: use a default high value (e.g., 40)
- `depth = 1`
- `gap_fraction = 0.0`

**After merging two profiles at an aligned position:**
- `prob[]` = normalized product of the two input distributions (Dirichlet posterior update)
- `quality` = f(quality_left, quality_right, alignment_posterior_at_this_position)
- `depth` = depth_left + depth_right
- `gap_fraction` = weighted average of input gap fractions, adjusted for newly introduced gaps

**For a newly introduced gap column:**
- `prob[]` = copy from the non-gapped profile's corresponding column
- `quality` = low (reflecting that this column's alignment is uncertain)
- `gap_fraction` = updated to reflect the gap

### 3.3 Quality-Aware Pair-HMM Emissions

This is the core extension of the existing EA pair-HMM. The existing implementation handles quality on one side via the Frith et al. (2010) marginalization:

```
P(a_i, b_j | M) = Σ_k subm[a_i][k] * P(k | b_j, Q_j)
```

For profile-profile alignment, both sides have uncertainty:

```
P(col_A, col_B | M) = Σ_s Σ_t P(s | col_A) * P(t | col_B) * subm[s][t]
```

Where `P(s | col_A) = col_A.prob[s]` is the probability distribution at that profile position.

**Connection to existing code:** The existing `qsubscore` table precomputes emission scores for all (base, quality) combinations. The profile extension replaces this with a runtime computation that is mathematically identical in structure but uses the profile probability vector instead of the quality-derived error distribution.

**Optimization for DNA (ALPHA_SIZE = 4):**
The double sum `Σ_s Σ_t` is only 4×4 = 16 terms. This can be computed as a matrix-vector product: emission = p_A^T · S · p_B, where S is the 4×4 substitution matrix and p_A, p_B are the 4-dimensional probability vectors. This is extremely cheap.

**Optimization for protein (ALPHA_SIZE = 20):**
The double sum is 20×20 = 400 terms. With sparse representations (top-k residues), this reduces to k² terms where k ≈ 3–4 in practice. Store profiles in sparse format and only iterate over nonzero entries.

### 3.4 Gap States: Emissions in Insert States

In the existing EA implementation, gap state emissions use uniform background:
```
P(a_i | X) = back[a_i] = 0.25  (DNA)
```

For profiles, the gap state emission when profile A emits into a gap (state X) is:
```
P(col_A | X) = Σ_s P(s | col_A) * back[s]
```

This weights the background by the column's distribution. For a well-conserved column (concentrated distribution), this is approximately `back[consensus]`. For a uniform column, this is approximately `Σ back[s]/ALPHA_SIZE = 1/ALPHA_SIZE`.

**Key insight:** A high-quality (conserved) column emitting into a gap state has low probability — the model naturally resists gapping conserved positions. A low-quality (uncertain) column emitting into a gap state has higher probability — the model is more willing to gap uncertain positions. This is the position-specific gap penalty behavior we want, emerging naturally from the probabilistic model.

### 3.5 Transition Probabilities as Gap Parameters

The existing EA pair-HMM uses:
```
           To M          To X          To Y
From M:    1 - g         g/2           g/2
From X:    1 - g         g             0
From Y:    1 - g         0             g

where g = base_error * indel_freq
```

For MSA, we parameterize transitions with:
- `δ` (delta): gap-open probability (M→X or M→Y)
- `ε` (epsilon): gap-extend probability (X→X or Y→Y)

```
           To M          To X          To Y
From M:    1 - 2δ        δ             δ
From X:    1 - ε         ε             0
From Y:    1 - ε         0             ε
```

**Phase 4 re-estimation:** After a complete progressive pass, `δ` and `ε` can be re-estimated from the alignment using Baum-Welch (the standard EM algorithm for HMMs). The E-step uses the forward-backward posteriors already computed during alignment. The M-step updates transition probabilities:

```
δ_new = (expected number of M→X transitions) / (expected number of transitions from M)
ε_new = (expected number of X→X transitions) / (expected number of transitions from X)
```

These expectations are computed from the posterior state probabilities summed over all progressive merge steps. This is the principled iterative bootstrapping loop — EM on a well-defined probabilistic model — that replaces hand-tuned gap penalties.

**Degeneracy prevention:** Unlike naive EM on raw alignments, the pair-HMM constrains parameters to be valid probabilities (0 < δ, ε < 1). The quality annotations provide further regularization: low-quality positions contribute less to parameter estimation. A weak Dirichlet prior on transitions (e.g., Beta(2, 100) for δ, Beta(2, 20) for ε) prevents collapse.

### 3.6 Profile Merging After Alignment

After the MEA alignment is computed, aligned positions from the two profiles are merged:

#### 3.6.1 Matched Positions (both profiles emit)

Given aligned positions `col_A` and `col_B`:

**Probability vector update (Dirichlet posterior product):**
```
merged.prob[s] ∝ col_A.prob[s]^(α_A) * col_B.prob[s]^(α_B)
```
where `α_A` and `α_B` are effective sample sizes derived from quality and depth. In practice, a simpler weighted average works well:
```
w_A = col_A.depth * col_A.quality
w_B = col_B.depth * col_B.quality
merged.prob[s] = (w_A * col_A.prob[s] + w_B * col_B.prob[s]) / (w_A + w_B)
```

**Quality update:**
```
alignment_posterior = P(match at this position | sequences)  // from forward-backward
merged.quality = combine(col_A.quality, col_B.quality, alignment_posterior)
```

A suitable combination function:
```
merged.quality = alignment_posterior * (col_A.quality + col_B.quality) / 2
```
This ensures that uncertain alignments (low posterior) produce low-quality merged columns, even if both input columns were high quality. This is the key error-propagation mechanism: early alignment mistakes get low quality and can be overwritten by later, more confident alignments.

**Depth and gap fraction:**
```
merged.depth = col_A.depth + col_B.depth
merged.gap_fraction = (col_A.depth * col_A.gap_fraction + col_B.depth * col_B.gap_fraction)
                      / merged.depth
```

#### 3.6.2 Gap Positions (one profile emits, the other gets a gap column)

When profile A emits at position i but profile B receives a gap:
```
merged.prob[s] = col_A.prob[s]    // from the emitting profile
merged.quality = col_A.quality * gap_posterior_penalty
merged.depth = col_A.depth + col_B.n_sequences
merged.gap_fraction = (col_A.depth * col_A.gap_fraction + col_B.n_sequences)
                      / merged.depth
```

The `gap_posterior_penalty` reflects that columns where one entire profile is gapped are inherently less certain. This is where the MUSCLE/FAMSA gap corrections emerge naturally: the gap_fraction tracks what fraction of sequences have gaps, and the quality is modulated by the alignment posterior.

### 3.7 Alignment Modes

The existing EA code supports semi-global alignment (free leading/trailing gaps in the database sequence). For MSA, we need:

- **Global alignment** for internal profile merges (both profiles should be fully aligned)
- **Semi-global alignment** for terminal gap handling (allow free terminal gaps when merging sequences of different lengths)

The existing forward/backward boundary conditions already implement the semi-global case. Global alignment requires penalizing all boundary gaps, which is a straightforward modification to the initialization conditions.

---

## 4. Mapping to Existing EA Codebase

### 4.1 Components to Reuse Directly

| EA Component | MSA Usage |
|---|---|
| `bpm()` / `bpm_256()` | Pairwise distance estimation for guide tree. Use edit distance between all sequence pairs. The AVX2 256-bit implementation handles sequences up to 255 residues; longer sequences use the serial 64-bit version. |
| `forward()` | Forward algorithm for profile-profile alignment. Core recurrence is identical; emission function is replaced with profile-aware version. |
| `backward()` | Backward algorithm. Same structural changes as forward. |
| `get_posterior()` | Posterior decoding. Unchanged — combines forward/backward matrices and divides by total score. |
| `max_posterior_alignment()` | MEA traceback. Unchanged — finds alignment maximizing sum of posterior match probabilities. |
| `random_score()` | Null model for log-odds scoring. Generalize to profile emissions for quality assessment of individual merge steps. |
| Log-space arithmetic (`logsum`, `logsum3`) | Used throughout. No changes needed. |
| Forward/backward consistency check | Reuse as correctness validation during development. |

### 4.2 Components Requiring Extension

#### 4.2.1 Emission Model (Critical Change)

**Current:** `get_qsubscore(quality, key)` — lookup from precomputed 3D table indexed by `(quality, base_a, base_b)`.

**New:** `profile_emission_match(col_A, col_B, submat)` — computes:
```c
float profile_emission_match(ProfileColumn *A, ProfileColumn *B, float submat[][ALPHA_SIZE]) {
    float score = 0.0;
    for (int s = 0; s < ALPHA_SIZE; s++) {
        if (A->prob[s] < MIN_PROB) continue;  // sparse optimization
        for (int t = 0; t < ALPHA_SIZE; t++) {
            if (B->prob[t] < MIN_PROB) continue;
            score += A->prob[s] * B->prob[t] * submat[s][t];
        }
    }
    return log(score);  // return in log space
}
```

**For DNA:** This is a 4×4 inner loop — negligible overhead.
**For protein with sparse profiles:** Typically 3–4 × 3–4 = 9–16 multiplications.

**Compatibility with existing single-sequence mode:** When both profiles are single sequences with quality 1.0, the emission reduces to `log(submat[a][b])`, which is the standard substitution score. When one profile is a single sequence with a Phred quality score, it reduces to the existing quality-aware emission. The generalization is strictly backward-compatible.

#### 4.2.2 Gap State Emissions

**Current:** `back[seq_a[i]]` — uniform 0.25 for DNA.

**New:** `profile_emission_gap(col, background)`:
```c
float profile_emission_gap(ProfileColumn *col, float background[]) {
    float score = 0.0;
    for (int s = 0; s < ALPHA_SIZE; s++) {
        score += col->prob[s] * background[s];
    }
    return log(score);
}
```

This naturally implements position-specific gap costs: conserved columns (concentrated prob) yield scores close to `log(back[consensus])`, while uncertain columns yield scores closer to `log(1/ALPHA_SIZE)`.

#### 4.2.3 Transition Probability Parameterization

**Current:** Parameterized by `base_error` and `indel_freq`, with `g = base_error * indel_freq`.

**New:** Parameterized directly by `delta` (gap-open) and `epsilon` (gap-extend):
```c
typedef struct {
    float delta;     // gap-open probability
    float epsilon;   // gap-extend probability
    float MM, MX, MY, XM, XX, XY, YM, YX, YY;  // log transition probabilities
} TransitionParams;

void init_transitions(TransitionParams *tp, float delta, float epsilon) {
    tp->delta = delta;
    tp->epsilon = epsilon;
    tp->MM = log(1.0 - 2.0 * delta);
    tp->MX = log(delta);
    tp->MY = log(delta);
    tp->XM = log(1.0 - epsilon);
    tp->XX = log(epsilon);
    tp->XY = LOG_ZERO;  // no X->Y transitions
    tp->YM = log(1.0 - epsilon);
    tp->YX = LOG_ZERO;  // no Y->X transitions
    tp->YY = log(epsilon);
}
```

**Default initialization:** `delta = 0.01`, `epsilon = 0.2` (reasonable starting point, equivalent to moderate gap-open penalty and geometric gap length distribution with mean 1/(1-ε) ≈ 1.25).

#### 4.2.4 DP Matrix Adaptation

The forward/backward recurrences remain structurally identical. The changes are:

1. Replace `get_qsubscore(qual_a[i], (seq_a[i] << 2) | seq_b[j])` with `profile_emission_match(&prof_A->columns[i], &prof_B->columns[j], submat)`
2. Replace `back[seq_a[i]]` with `profile_emission_gap(&prof_A->columns[i], background)`
3. Replace hardcoded transition log-probabilities with `TransitionParams` struct values

The DP matrix dimensions change from `(len_a+1) × (len_b+1)` for sequences to `(prof_A->length+1) × (prof_B->length+1)` for profiles — structurally identical.

### 4.3 New Components to Build

#### 4.3.1 Profile Data Structure and I/O

- `Profile` struct with array of `ProfileColumn`
- `profile_from_sequence()`: Create a leaf profile from a raw sequence (+ optional quality scores)
- `profile_merge()`: Merge two profiles given an MEA alignment path
- `profile_to_msa()`: Convert final root profile back to a gapped MSA for output
- `profile_column_entropy()`: Compute Shannon entropy per column (for diagnostics/regularization)

#### 4.3.2 Guide Tree Construction

- `distance_matrix_bpm()`: Compute all-pairs edit distances using existing BPM functions
- `upgma_tree()` or `slink_tree()`: Build guide tree from distance matrix
- `Tree` struct: binary tree with leaf→sequence and internal→merge-order mapping

#### 4.3.3 Progressive Alignment Driver

- `progressive_align()`: Walk the guide tree bottom-up, at each node:
  1. Retrieve child profiles
  2. Call forward-backward pair-HMM
  3. Call MEA traceback
  4. Call profile merge
  5. Store merged profile at the node

#### 4.3.4 Baum-Welch Parameter Re-estimation

- `accumulate_transition_counts()`: From forward-backward posteriors at each merge step, accumulate expected transition counts
- `reestimate_transitions()`: M-step — compute new δ, ε from accumulated counts
- `em_iterate()`: Outer loop — progressive align → accumulate → reestimate → repeat

#### 4.3.5 Substitution Matrix Handling

- `load_substitution_matrix()`: Read BLOSUM62, MIQS, or identity matrix
- `submat_to_joint_prob()`: Convert log-odds substitution matrix to joint probability matrix (needed for pair-HMM emissions)
- For DNA: simple 4×4 matrix parameterized by match/mismatch rates
- For protein: load standard matrices, convert from log-odds to probability space

#### 4.3.6 Sequence Tracking for Output

The final output must be a full MSA, not just a root profile. We need to track which gaps were inserted at each step so the final alignment can be reconstructed:

- **Edit string (e-string) approach (MUSCLE-style):** At each merge step, record a compact representation of which positions received gap insertions. The final alignment is reconstructed by composing e-strings along the tree.
- **Gapped sequence representation (FAMSA-style):** Store each sequence as (symbols, gap_counts) pairs, updating in-place as gaps are inserted. More memory-efficient for large families.

The FAMSA-style representation is recommended for scalability. The profile's `ProfileColumn` array represents the consensus view; a parallel per-sequence gap tracking array enables output.

---

## 5. Computational Complexity Analysis

### 5.1 Per-Merge Step

| Operation | Cost | Notes |
|---|---|---|
| Forward algorithm | O(L_A × L_B) | Identical to standard DP |
| Backward algorithm | O(L_A × L_B) | Identical to standard DP |
| Posterior decoding | O(L_A × L_B) | Element-wise combination |
| MEA traceback | O(L_A × L_B) | Second DP pass over posteriors |
| Profile merge | O(L_merged) | Linear scan of alignment |
| **Total per merge** | **O(L_A × L_B)** | **Same as standard profile alignment** |

**Constant factor overhead vs. standard DP:**
- DNA: Emission computation is ~16 multiplications instead of 1 table lookup. Offset by avoiding the quality table indirection.
- Protein (sparse, k=4): ~16 multiplications instead of 1. Still negligible compared to DP overhead.
- Forward-backward is ~2× Viterbi cost (running both directions). But we gain MEA alignment quality.

### 5.2 Full Pipeline

| Phase | Cost | Notes |
|---|---|---|
| Pairwise distances | O(N² × L / w) | BPM with word size w=64 or w=256 (AVX2) |
| Guide tree | O(N²) for UPGMA/SLINK | O(N log N) with medoid tree heuristic |
| Progressive alignment | O(N × L²) | N-1 merge steps, each O(L²) |
| Parameter re-estimation | O(N × L²) per iteration | Same cost as progressive alignment |
| **Total** | **O(N² × L/w + N × L²)** | **Same asymptotic as FAMSA/Kalign** |

### 5.3 Memory

| Component | Memory | Notes |
|---|---|---|
| Distance matrix | O(N²) or O(N) with SLINK | SLINK processes incrementally |
| DP matrices (forward/backward) | O(L²) × 3 states × 2 directions | Can be reduced to O(L) with checkpointing |
| Profile at each tree node | O(L × ALPHA_SIZE) | Profiles freed after merge |
| Sequence tracking | O(N × L_final) | For output reconstruction |
| **Peak memory** | **O(N × L + L²)** | **Competitive with FAMSA** |

---

## 6. Protein Support

### 6.1 Alphabet Extension

The existing implementation is DNA-only (ALPHA_SIZE = 4). Protein support requires:

- `ALPHA_SIZE = 20` (standard amino acids)
- `ALPHA_WITH_GAP = 21` (for profile representations including gap as a character)
- Encoding map: A=0, R=1, N=2, D=3, C=4, Q=5, E=6, G=7, H=8, I=9, L=10, K=11, M=12, F=13, P=14, S=15, T=16, W=17, Y=18, V=19
- Unknown residue X: uniform distribution `prob[s] = 1/20` for all s

### 6.2 Substitution Matrices

For DNA:
- Simple parameterized matrix: `match = (1-3μ)`, `mismatch = μ` for substitution rate μ
- Or the existing EA `subm[4][4]` construction

For protein:
- BLOSUM62 (standard), MIQS (recommended by FAMSA for distantly related sequences)
- Convert from log-odds to probability space: `P(s,t) = p_s * p_t * exp(S(s,t) * λ)` where λ is the scaling factor
- Background frequencies from the chosen matrix (e.g., Robinson-Robinson for BLOSUM)

### 6.3 Sparse Profile Optimization

For protein, most columns have at most 3–4 significant residues. Store profiles in sparse format:

```c
typedef struct {
    uint8_t residue;
    float prob;
} SparseEntry;

typedef struct {
    SparseEntry entries[MAX_SPARSE];  // MAX_SPARSE = 5 typically
    uint8_t n_entries;
    float quality;
    uint32_t depth;
    float gap_fraction;
} SparseProfileColumn;
```

Emission computation with sparse profiles:
```c
float sparse_emission_match(SparseProfileColumn *A, SparseProfileColumn *B, float submat[][20]) {
    float score = 0.0;
    for (int i = 0; i < A->n_entries; i++)
        for (int j = 0; j < B->n_entries; j++)
            score += A->entries[i].prob * B->entries[j].prob
                     * submat[A->entries[i].residue][B->entries[j].residue];
    return log(score);
}
```

With k=4 entries on each side, this is 16 multiply-adds — very fast.

---

## 7. Implementation Plan

### Phase 1: Core Profile Pair-HMM (MVP)

**Goal:** Align two profiles using the quality-annotated pair-HMM.

1. Define `ProfileColumn` and `Profile` data structures
2. Implement `profile_from_sequence()` for leaf profiles (DNA first)
3. Implement `profile_emission_match()` and `profile_emission_gap()`
4. Modify `forward()` to accept `Profile*` inputs instead of raw sequences
5. Modify `backward()` similarly
6. Implement `profile_merge()` using MEA alignment path
7. **Validation:** Align two sequences via profiles and verify that the result matches the existing EA pairwise alignment (when profiles are single sequences, results must be identical)

### Phase 2: Progressive MSA Pipeline

**Goal:** Align N sequences progressively using a guide tree.

1. Implement `distance_matrix_bpm()` using existing BPM functions
2. Implement UPGMA or SLINK guide tree construction
3. Implement `progressive_align()` tree-walking driver
4. Implement sequence tracking for MSA output reconstruction (e-string or FAMSA-style gap tracking)
5. Implement FASTA/A2M output
6. **Validation:** Benchmark on BAliBASE RV11/RV12 (closely related sequences). Compare SP and TC scores against MUSCLE, MAFFT, Kalign3

### Phase 3: Parameter Re-estimation (EM)

**Goal:** Learn gap parameters from data.

1. Implement transition count accumulation during progressive alignment
2. Implement Baum-Welch M-step for δ, ε
3. Implement EM outer loop with convergence detection
4. Add Dirichlet prior on transition parameters for regularization
5. **Validation:** Show that EM-estimated parameters improve alignment quality vs. default parameters on diverse benchmarks. Show convergence in 2–3 iterations.

### Phase 4: Protein Support

**Goal:** Extend to protein sequences.

1. Generalize alphabet handling (compile-time or runtime ALPHA_SIZE)
2. Implement substitution matrix loading and log-odds → probability conversion
3. Implement sparse profile representation for protein
4. **Validation:** Benchmark on full BAliBASE, PREFAB, SABmark, OXBench

### Phase 5: Scalability

**Goal:** Handle large families (10K+ sequences).

1. Implement SLINK for O(N) memory guide tree construction
2. Implement FAMSA-style gapped sequence representation for O(1) gap insertion
3. Parallelize independent subtree merges (OpenMP)
4. Parallelize pairwise distance computation (existing BPM already supports this conceptually)
5. Implement DP memory optimization (linear-space forward-backward if needed)
6. **Validation:** Benchmark on HomFam/extHomFam large families. Compare speed and accuracy against FAMSA, Clustal Omega, MAFFT

### Phase 6: Advanced Features

1. Position-specific transition probabilities (modulate δ, ε by local quality/conservation)
2. Column entropy regularization for EM
3. Transitivity consistency checking
4. Integration of predicted secondary structure as quality modulation (for protein)

---

## 8. Benchmarking Strategy

### 8.1 Accuracy Benchmarks

| Benchmark | Size Range | Purpose |
|---|---|---|
| BAliBASE 3.0 | 4–100 sequences | Standard accuracy (6 categories) |
| PREFAB 4.0 | ~1000 pairs | Pairwise accuracy across identity ranges |
| SABmark | Twilight zone | Low-identity accuracy |
| OXBench | ~400 families | General accuracy |
| HomFam / extHomFam | 200–415K sequences | Large family accuracy and scalability |

**Metrics:** Sum-of-pairs (SP) score, Total Column (TC) score via QSCORE.

### 8.2 Comparisons

| Method | Category | Expected Comparison |
|---|---|---|
| FAMSA / FAMSA2 | Fast progressive | Similar speed, better accuracy (quality propagation) |
| MAFFT (default, L-INS-i) | Progressive + iterative | Similar accuracy, better speed than L-INS-i |
| Clustal Omega | Progressive (HMM) | Better accuracy on large families |
| Kalign3 | Fast progressive | Similar speed, better accuracy |
| ProbCons / MSAProbs | Probabilistic consistency | Similar accuracy, much better speed |
| MUSCLE5 | Progressive | Better accuracy (quality propagation + MEA) |

### 8.3 Key Claims to Validate

1. **Quality propagation improves accuracy:** ProbMSA with quality propagation vs. ProbMSA with flat profiles (ablation study)
2. **EM parameter estimation works:** ProbMSA with EM vs. ProbMSA with hand-tuned defaults, across diverse benchmarks
3. **No hand-tuning needed:** Single set of EM-estimated parameters competitive with per-benchmark-optimized parameters of other tools
4. **Speed:** Within 2–3× of FAMSA on large families; orders of magnitude faster than ProbCons/MSAProbs
5. **The MUSCLE/FAMSA corrections are subsumed:** Show that ProbMSA without any explicit gap corrections matches or exceeds tools that use them

---

## 9. Key Design Decisions and Rationale

### 9.1 Why Pair-HMM Instead of Standard DP

Standard DP with substitution scores and affine gap penalties is a special case of the 3-state pair-HMM with Viterbi decoding. The pair-HMM framework gives us:

- **Forward-backward posteriors** at no additional asymptotic cost (2× constant factor)
- **MEA alignment** instead of Viterbi — empirically better for MSA (ProbCons demonstrated this)
- **Natural parameterization** of gap penalties as transition probabilities (constrained to valid probability space)
- **EM estimation** of parameters via standard Baum-Welch
- **Quality-aware emissions** via marginalization — already implemented in EA

There is no downside to the pair-HMM framework. It strictly generalizes standard DP while enabling the probabilistic features we need.

### 9.2 Why One-Sided Quality Extends to Two-Sided Naturally

The existing EA implementation marginalizes over base uncertainty on one side:
```
P(a, b|M) = Σ_k submat[a][k] * P(k|b, Q_b)
```

The two-sided generalization is:
```
P(col_A, col_B|M) = Σ_s Σ_t P(s|col_A) * P(t|col_B) * submat[s][t]
```

The first is a special case of the second where `P(s|col_A) = δ(s, a)` (delta function at the observed residue). The implementation change is minimal: replace the table lookup with a small nested loop.

### 9.3 Why Dirichlet Posterior for Merging (Not Simple Averaging)

Simple frequency averaging treats all sequences equally regardless of alignment quality. The Dirichlet posterior product weights by effective sample size (quality × depth), so:

- Confidently aligned columns contribute more to the merged distribution
- Columns with alignment ambiguity contribute less
- A single high-quality sequence can outweigh a large but poorly aligned profile

In practice, the weighted average approximation `merged = (w_A * col_A + w_B * col_B) / (w_A + w_B)` captures the key behavior and is computationally simpler. Use the full Dirichlet product only if the approximation proves insufficient in benchmarks.

### 9.4 Why This Subsumes MUSCLE/FAMSA Gap Corrections

MUSCLE's position-specific gap correction (scale gap-open by fraction of non-gap sequences) approximates the true SP score contribution of gap penalties. In ProbMSA:

- **Existing gaps in the profile** manifest as `gap_fraction > 0` in the ProfileColumn
- **The emission model** naturally down-weights gappy columns (uncertain distributions → higher gap state probability relative to match)
- **The quality annotation** is low at positions where previous alignment was ambiguous, further reducing the effective gap penalty

FAMSA's set-size scaling addresses the tendency of large alignments to widen. In ProbMSA:

- **Profiles grow in depth** as more sequences are added
- **Quality tends to decrease** as more sequences are added (unless conservation is strong), naturally increasing gap tolerance only where warranted
- **The EM re-estimation** adapts transition parameters to the actual gap frequency in the growing alignment, automatically adjusting for set size

These mechanisms are not ad hoc corrections bolted onto a score-based aligner — they emerge from the probabilistic model.

---

## 10. Risks and Mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| Profile emission computation is too slow for protein (20×20 inner loop in O(L²) DP cells) | Performance degradation | Sparse profile representation limits inner loop to ~16 multiplications. Profile emission caching for repeated column patterns. |
| EM does not converge or converges to poor parameters | Accuracy regression | Dirichlet prior on transitions. Bounded parameter space (δ ∈ [0.001, 0.1], ε ∈ [0.05, 0.5]). Fall back to defaults if EM diverges. |
| Quality propagation amplifies early errors instead of correcting them | Accuracy regression | Quality is multiplicatively attenuated at each step (alignment posterior < 1.0). Early errors get progressively lower quality. Validate with ablation study. |
| Memory overhead of profile representation for very large families | OOM on large inputs | FAMSA-style gapped sequence representation. Profiles freed after merge. Only O(active_nodes) profiles in memory at any time. |
| Forward-backward 2× overhead makes Phase 3 progressive alignment 2× slower than competitors | Speed inferiority | Inner loop is memory-bound (DP matrix access), not compute-bound. The 2× compute overhead may be partially hidden. MEA accuracy gain justifies modest speed cost. |

---

## 11. File Structure

```
probmsa/
├── src/
│   ├── main.c                    # CLI entry point
│   ├── pairhmm.c / pairhmm.h    # Core pair-HMM (forward, backward, posterior, MEA)
│   │                              # Extended from EA with profile emissions
│   ├── profile.c / profile.h    # Profile data structure, merge, I/O
│   ├── emission.c / emission.h  # Emission probability computation
│   │                              # (profile-profile match, gap emissions)
│   ├── transitions.c / .h       # Transition parameterization, Baum-Welch re-estimation
│   ├── guide_tree.c / .h        # Distance matrix, UPGMA/SLINK tree construction
│   ├── progressive.c / .h       # Progressive alignment driver
│   ├── sequence.c / .h          # Sequence I/O (FASTA parsing, MSA output)
│   ├── bpm.c / bpm.h            # Bit-parallel Myers (reuse from EA)
│   ├── submat.c / submat.h      # Substitution matrix loading and conversion
│   ├── logmath.c / logmath.h    # Log-space arithmetic (reuse from EA)
│   └── utils.c / utils.h        # Memory allocation, error handling
├── matrices/
│   ├── blosum62.txt
│   ├── miqs.txt
│   └── dna_default.txt
├── tests/
│   ├── test_profile_emission.c   # Unit tests for emission model
│   ├── test_pairhmm_profile.c    # Verify profile pair-HMM matches EA on single sequences
│   ├── test_progressive.c        # Small MSA correctness tests
│   └── test_baum_welch.c         # EM convergence tests
├── benchmarks/
│   ├── run_balibase.sh
│   ├── run_prefab.sh
│   └── run_homfam.sh
├── Makefile
└── README.md
```

---

## 12. Success Criteria

### 12.1 Phase 1 (Core Profile Pair-HMM)

- [ ] Profile-profile forward-backward produces identical results to EA when profiles are single sequences
- [ ] Forward and backward scores agree within tolerance (existing consistency check)
- [ ] MEA alignment of two profiles produces biologically sensible result on simple test cases

### 12.2 Phase 2 (Progressive MSA)

- [ ] BAliBASE RV11 SP score ≥ 80% (comparable to MUSCLE/MAFFT default)
- [ ] Correct MSA output reconstruction (all sequences recoverable from gapped alignment)
- [ ] Processes 1000 sequences of length 300 in < 60 seconds

### 12.3 Phase 3 (EM Parameter Estimation)

- [ ] EM converges in ≤ 5 iterations on 90% of test families
- [ ] EM-estimated parameters match or exceed hand-tuned defaults on BAliBASE aggregate score
- [ ] No parameter divergence (δ and ε stay within valid bounds) on any test case

### 12.4 Phase 4 (Protein)

- [ ] BAliBASE aggregate SP score within 2% of MAFFT-default
- [ ] PREFAB SP score within 2% of MUSCLE
- [ ] SABmark SP score competitive with FAMSA or better

### 12.5 Phase 5 (Scalability)

- [ ] Process 10,000 sequences in < 10 minutes (comparable to FAMSA)
- [ ] Process 100,000 sequences in < 2 hours
- [ ] Peak memory < 16 GB for 100,000 sequences of average length 300

---

## 13. Appendix: Mathematical Details

### 13.1 Pair-HMM Forward Recurrence with Profile Emissions

Let `A[1..m]` and `B[1..n]` be two profiles. Let `e_M(i,j)` be the log match emission for aligning column i of A to column j of B. Let `e_X(i)` and `e_Y(j)` be the log gap emissions for A column i and B column j respectively.

```
F_M[i][j] = logsum3(
    F_M[i-1][j-1] + log(1 - 2δ),
    F_X[i-1][j-1] + log(1 - ε),
    F_Y[i-1][j-1] + log(1 - ε)
) + e_M(i, j)

F_X[i][j] = logsum(
    F_M[i-1][j] + log(δ),
    F_X[i-1][j] + log(ε)
) + e_X(i)

F_Y[i][j] = logsum(
    F_M[i][j-1] + log(δ),
    F_Y[i][j-1] + log(ε)
) + e_Y(j)
```

Where:
```
e_M(i, j) = log( Σ_s Σ_t A[i].prob[s] * B[j].prob[t] * P_joint(s,t) )
e_X(i) = log( Σ_s A[i].prob[s] * background[s] )
e_Y(j) = log( Σ_t B[j].prob[t] * background[t] )
```

And `P_joint(s,t)` is the joint probability of aligned residues s and t (derived from the substitution matrix).

### 13.2 Posterior Match Probability

```
P(M at i,j | A, B) = exp(F_M[i][j] + B_M[i][j] - e_M(i,j) - total_score)
```

Where `total_score = logsum3(F_M[m][n], F_X[m][n], F_Y[m][n])` and `B_M` is the backward match matrix.

### 13.3 Baum-Welch Transition Re-estimation

Expected transition counts from a single merge step:
```
E[M→X] = Σ_i Σ_j exp(F_M[i][j] + log(δ) + e_X(i+1) + B_X[i+1][j] - total_score)
E[M→M] = Σ_i Σ_j exp(F_M[i][j] + log(1-2δ) + e_M(i+1,j+1) + B_M[i+1][j+1] - total_score)
```

Accumulated over all N-1 merge steps in the progressive alignment, then:
```
δ_new = (Σ E[M→X] + Σ E[M→Y]) / (Σ E[M→M] + Σ E[M→X] + Σ E[M→Y] + prior_counts)
ε_new = (Σ E[X→X] + Σ E[Y→Y]) / (Σ E[X→X] + Σ E[X→M] + Σ E[Y→Y] + Σ E[Y→M] + prior_counts)
```

### 13.4 Quality Update Function

For matched positions with alignment posterior `p`:
```
q_merged = p * (w_A * q_A + w_B * q_B) / (w_A + w_B)

where w_A = depth_A, w_B = depth_B
```

This satisfies the desiderata:
- If `p ≈ 1` (confident alignment): quality is the depth-weighted average of inputs
- If `p ≈ 0` (uncertain alignment): quality drops to near zero regardless of inputs
- If both inputs are high quality: merged quality is high (good)
- If one input is low quality: merged quality is diluted (good)

For gap positions:
```
q_gap = q_emitting * min(p_gap, 0.5)
```

Where `p_gap` is the posterior probability of the gap state at that position. The `min(0.5)` cap prevents gap positions from ever being high quality.

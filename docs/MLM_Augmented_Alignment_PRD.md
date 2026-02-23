# PRD: MLM-Augmented Sequence Alignment

**Using Masked Language Model Predictions to Build Context-Aware Substitution Scoring for Protein Alignment**

Target Integration: Kalign Multiple Sequence Aligner
February 2026 · Draft v1.0

| | |
|---|---|
| **Status** | Draft – For Review |
| **Classification** | Research & Development |
| **Target Conference** | Bioinformatics / NAR / ISMB |

---

## 1. Executive Summary

This document describes a novel approach to protein multiple sequence alignment (MSA) that augments traditional substitution matrices (BLOSUM/PAM) with **masked language model (MLM) probability distributions** from protein language models such as ESM-2. Unlike existing PLM-based alignment methods which rely on high-dimensional embedding vectors (1024–1280 dimensions), our approach uses the **20-dimensional amino acid prediction distributions** produced by the masked token prediction head. This signal is qualitatively different from embeddings, approximately 60× cheaper to compare, and captures what each position *expects* evolutionarily rather than what it *is* in context.

The method combines three complementary signals into a single scoring function:

1. Traditional BLOSUM substitution scores for amino acid similarity
2. Jensen-Shannon divergence between MLM distributions for positional evolutionary constraint matching
3. Information-content weighting so that the MLM signal dominates where it is informative and defers to BLOSUM where it is not

As a byproduct, MLM entropy provides position-specific gap penalties at no additional computational cost.

The target integration platform is Kalign, a fast progressive MSA tool, with the goal of exceeding MAFFT accuracy on divergent sequences while maintaining competitive speed.

---

## 2. Problem Statement

### 2.1 The Limitation of Static Substitution Matrices

All traditional MSA tools (MAFFT, MUSCLE, Clustal Omega, Kalign) score residue pairs using fixed substitution matrices like BLOSUM62. These matrices encode the average log-odds of observing one amino acid substituted for another across a broad set of aligned protein families. The fundamental limitation is that they are **context-free**: BLOSUM62 assigns the same score to a Leucine→Isoleucine substitution regardless of whether it occurs in a hydrophobic core, at a catalytic site, or in a disordered loop.

This context-blindness leads to two failure modes:

- **False positives:** Two residues with identical amino acids receive a high score even when they play completely different structural/functional roles and should not be aligned.
- **False negatives:** Two residues with different amino acids receive a low score even when they play equivalent roles (e.g., a catalytic histidine in one protein and a catalytic serine in another).

These failures are most acute in the twilight zone (20–40% sequence identity), where MAFFT's L-INS-i mode is considered state of the art but still produces many incorrect column assignments.

### 2.2 What Existing PLM Methods Do

Several recent methods replace substitution matrices with PLM embedding similarity (PLMAlign, EBA, pLM-BLAST, DEDAL, vcMSA, learnMSA2, ARIES). All of these use the high-dimensional per-residue embedding vectors (1024–1280 dimensions) and compare them using cosine similarity, Euclidean distance, or learned bilinear forms. While effective, this approach has practical drawbacks:

- **Computational cost:** Comparing two 1280-dimensional vectors for every cell of the DP matrix is expensive. For two sequences of length 300, this is ~115 million floating-point operations just for scoring.
- **Memory:** Storing embeddings for 1000 sequences of average length 300 requires ~1.5 GB.
- **Integration difficulty:** Embedding-based scoring does not map naturally onto existing DP-based alignment frameworks, which are optimized for scalar substitution lookups.

Critically, **no existing method uses the MLM probability distributions** (the 20-dimensional output of the masked token prediction head) as a substitution scoring signal for alignment.

---

## 3. Proposed Solution: MLM-Augmented Scoring

### 3.1 Core Concept

When a protein language model like ESM-2 masks a single position and predicts a probability distribution over the 20 standard amino acids, that distribution encodes the evolutionary constraints at that position given the full sequence context. This distribution answers the question: "What does evolution expect at this position?"

We propose comparing these distributions between positions in different sequences to determine whether two positions share the same evolutionary role. Two positions with similar MLM distributions are under the same evolutionary constraints and should therefore be aligned together, even if the actual amino acids at those positions differ.

### 3.2 The Three-Component Scoring Function

The core of the method is a hybrid scoring function that combines traditional and AI-derived signals:

```
S(i, j) = α · BLOSUM62[ seq₁[i], seq₂[j] ]  +  γ · (1 − JS( MLM₁[i], MLM₂[j] )) · C(i, j)

where:
    C(i, j) = ( IC(i) + IC(j) ) / 2
    IC(i)   = log₂(20) − H( MLM₁[i] )
```

Where the components are:

| Component | Description |
|---|---|
| **BLOSUM62[a, b]** | Standard substitution matrix score. Captures amino acid physicochemical similarity. A 20×20 lookup table — essentially free to compute. |
| **JS(P, Q)** | Jensen-Shannon divergence between two 20-dimensional MLM probability distributions. Bounded [0, 1], symmetric, and information-theoretically principled. A value of 0 means the two positions predict identical amino acid distributions; 1 means completely disjoint predictions. |
| **C(i, j)** | Information-content weighting. The average information content at both positions. Ranges from 0 (both positions are unconstrained, nearly uniform distributions) to log₂(20) ≈ 4.32 (both positions have deterministic predictions). This ensures the MLM signal is strongest where it matters most. |
| **α, γ** | Mixing weights, tunable on benchmark data (BAliBASE, SABmark, PREFAB). Expected optimal values: α near 1.0, γ in the range 0.5–2.0 depending on sequence divergence. |

### 3.3 Why Each Component Is Necessary

Each of the three components captures a qualitatively different biological signal:

| Scenario | BLOSUM | MLM Similarity | Correct Action |
|---|---|---|---|
| Same AA, same role (Leu in core ↔ Leu in core) | ✅ High (+4) | ✅ High (low JS) | Align — both signals agree |
| Diff AA, same role (His catalytic ↔ Ser catalytic) | ❌ Low (−1) | ✅ High (both predict catalytic residues) | Align — MLM rescues this case |
| Same AA, diff role (Gly in turn ↔ Gly in loop) | ❌ High (+6) | ✅ Low (sharp vs flat distribution) | Don't align — MLM prevents false positive |
| Diff AA, diff role (unconstrained positions) | Low | Weak signal (IC weighting → ~0) | Use BLOSUM as fallback |

The second row is the **key advantage over all existing methods**: detecting equivalent functional positions across divergent sequences where the amino acids are completely different. The third row is equally important: **preventing false positive alignments** where identical residues are playing unrelated roles.

### 3.4 Information-Content Weighting: The Critical Innovation

Not all MLM predictions are equally informative. Consider two positions that both produce nearly uniform distributions over all 20 amino acids. Their JS divergence will be close to zero (they look similar), but this similarity is meaningless — both positions are simply unconstrained.

The information-content weighting term C(i, j) solves this by modulating the MLM contribution based on how constrained the positions are:

- **Two highly constrained positions that agree** (e.g., both strongly predict histidine) → High IC, low JS divergence → strong positive MLM score. The model is confident these positions play the same role.
- **Two highly constrained positions that disagree** (e.g., one predicts histidine, the other predicts glycine) → High IC, high JS divergence → strong negative MLM score. The model is confident these positions should NOT be aligned.
- **Two unconstrained positions** → Low IC → MLM score approaches zero regardless of JS divergence. BLOSUM handles these positions, as it should.

This weighting creates a natural and principled handoff between the two scoring systems: the MLM signal takes over from BLOSUM exactly when it has the most to offer (constrained positions in divergent sequences) and yields to BLOSUM exactly when it has nothing to say (unconstrained positions).

### 3.5 Position-Specific Gap Penalties (Free Byproduct)

The MLM entropy at each position directly provides biologically motivated gap penalties at no additional computational cost:

```
gap_open(i) = base_gap_open × (1 + λ · IC(i))
```

Positions with high information content (constrained, low entropy) receive higher gap penalties because inserting a gap at a conserved site is biologically costly. Positions with low information content (unconstrained, high entropy, typically loops and disordered regions) receive lower gap penalties, reflecting that indels in these regions are common and tolerated.

No existing MSA tool derives position-specific gap penalties from a language model in this way.

---

## 4. Why This Approach Is Uniquely Advantageous

### 4.1 MLM Distributions vs. Embeddings: A Fundamental Distinction

Every existing PLM-based alignment method uses embeddings (the hidden-state vectors from intermediate or final transformer layers). Our approach uses a different, complementary output: the masked-token probability distributions from the model's prediction head. These are fundamentally different signals:

| Property | Embeddings (what others use) | MLM Distributions (our approach) |
|---|---|---|
| **Dimensionality** | 1024–1280 per residue | 20 per residue |
| **Comparison cost** | ~1280 FLOPs per cell | ~20 FLOPs per cell (~60× cheaper) |
| **Memory for 1000 seqs** | ~1.5 GB | ~24 MB |
| **Semantic meaning** | What the residue IS in context | What the position EXPECTS evolutionarily |
| **Interpretability** | Opaque (high-dimensional) | Directly interpretable (amino acid probabilities) |
| **Yields gap penalties?** | Not directly | Yes, from entropy |
| **Natural BLOSUM hybrid?** | Difficult (scale mismatch) | Natural (both are residue-pair scores) |

### 4.2 Breaking the Chicken-and-Egg Problem

Traditional MSA tools face a circular dependency: you need a good alignment to build a good profile, but you need a good profile to build a good alignment. Tools like MAFFT and MUSCLE address this through iterative refinement, but the initial alignment still starts from nothing.

MLM distributions **break this cycle entirely**. The language model provides per-position evolutionary profiles from single sequences — no alignment required. Each sequence is processed independently by ESM-2, which has internalized evolutionary patterns from training on hundreds of millions of sequences. This gives you something analogous to a profile-based prior *before a single pair of sequences has been aligned*. The very first progressive alignment is already informed by deep evolutionary knowledge from the broader protein universe.

### 4.3 Practical Integration Advantage

Because MLM distributions produce a scalar similarity score for each position pair, they integrate directly into existing DP-based alignment frameworks. The implementation requires changing only the scoring function — the alignment algorithm, guide tree construction, and iterative refinement can remain unchanged. Furthermore, ESM-2's MIT license and ONNX exportability enable a **pure-C deployment**: the model runs via ONNX Runtime's C API with no Python dependency at runtime (see §8.7). This is in contrast to methods like vcMSA (which abandons DP entirely) or learnMSA2 (which requires a differentiable HMM framework and A100 GPUs).

---

## 5. Novelty Assessment

A thorough literature review (as of February 2026) reveals the following landscape of PLM-aided alignment methods and confirms the novelty of our specific approach.

### 5.1 Existing Work (What Is Not Novel)

- **Embedding-based pairwise alignment:** PLMAlign (Liu et al., Nat Commun 2024), EBA (Pantolini et al., Bioinformatics 2024), pLM-BLAST (Kaminski et al., 2023), DEDAL (Llinares-López et al., 2023)
- **Embedding-based MSA:** vcMSA (McWhite & Singh, Genome Research 2023), learnMSA2 (Becker & Stanke, Bioinformatics 2024), ARIES (bioRxiv, January 2026)
- **MLM masked marginals for variant effect prediction:** ESM-1v (Meier et al., 2021), VESPA (Marquet et al., 2022)

### 5.2 What Is Novel In Our Approach

| Component | Novel? | Notes |
|---|---|---|
| JS divergence between MLM distributions as a substitution score | ✅ Yes | No prior work found |
| Information-content weighting for adaptive BLOSUM/MLM handoff | ✅ Yes | Novel formulation |
| BLOSUM + MLM hybrid scoring function | ✅ Yes | Unique combination |
| MLM entropy → position-specific gap penalties | ✅ Yes | Not in any alignment paper |
| 20d vs 1280d computational efficiency argument | ✅ Yes | Strong practical advantage |
| Pure-C deployment via ONNX Runtime (no Python) | ✅ Yes | No existing PLM alignment tool ships as a standalone C binary |
| PLM-aided alignment in general | ❌ No | Active field, many papers |
| Embeddings as substitution score | ❌ No | PLMAlign, EBA, DEDAL, etc. |

### 5.3 Publication Angle

> *"Everyone uses the embeddings from protein language models for alignment. We use a different, overlooked output: the masked-token predictions. These 20-dimensional probability distributions capture what each position expects evolutionarily — a qualitatively different signal from what embeddings capture. They are 60× cheaper to compare, naturally combine with BLOSUM matrices, provide information-content weighting for free, and yield position-specific gap penalties as a byproduct."*

---

## 6. System Architecture

### 6.1 Pipeline Overview

The full alignment pipeline consists of the following stages:

1. **Input:** N unaligned protein sequences in FASTA format.
2. **ESM-2 Forward Pass:** Run each sequence through ESM-2 once (without masking, using the pseudo-likelihood approximation). This single forward pass per sequence yields both the per-residue MLM probability distributions (20 values each) and the mean sequence embedding (for guide tree construction). In production, ESM-2 runs natively in C via ONNX Runtime (see §8.7); during prototyping, PyTorch is used.
3. **Pairwise Distance Estimation:** Compute pairwise distances between all sequences using mean embedding cosine distance (or fast k-mer distance for very large sets). Build a UPGMA/NJ guide tree.
4. **Progressive Alignment:** Align sequences/profiles progressively following the guide tree, using the three-component scoring function S(i, j) with MLM-derived position-specific gap penalties.
5. **Iterative Refinement:** Optionally re-score each column by its average MLM consensus, identify low-confidence regions, and re-align those regions with adjusted parameters.
6. **Output:** Multiple sequence alignment in standard formats (FASTA, CLUSTAL, PHYLIP).

### 6.2 Pseudo-Likelihood Approximation

A naive implementation would mask each position individually and run a forward pass for each, requiring L forward passes per sequence of length L. This is impractical.

Instead, we use the **pseudo-likelihood approximation**: run the model once on the unmasked sequence and take the output logits at each position as an approximation of the masked prediction. This is widely used for mutation effect scoring (ESM-1v, VESPA) and is empirically very close to true masked marginals for the purposes of generating probability distributions. **One forward pass per sequence gives us everything we need.**

### 6.3 Computational Complexity

| Stage | Cost | Notes |
|---|---|---|
| **ESM-2 inference** | O(N) forward passes | ~100 seq/sec on RTX 3090 for 650M model; ~1000 seq/sec for 150M model on CPU |
| **MLM storage** | 20 × L floats per sequence | ~24 MB for 1000 seqs of length 300 |
| **Scoring each DP cell** | ~40 FLOPs (BLOSUM lookup + JS on 20d) | vs. ~1280 FLOPs for embedding cosine |
| **Total overhead vs. standard Kalign** | ESM-2 inference + ~2× scoring cost | ESM-2 is one-time cost; DP scoring is modest |

---

## 7. Detailed: How MLM Scores Complement Substitution Matrices

### 7.1 What the MLM Distribution Encodes

When ESM-2 processes a sequence and outputs logits at each position (via the pseudo-likelihood approach), the softmax of those logits gives a 20-dimensional probability distribution. Different types of positions produce characteristically different distributions:

- **Highly conserved catalytic site:** Very sharp, low-entropy distribution. Example: P(His) = 0.85, P(Asn) = 0.05, P(Asp) = 0.04, all others near zero. The model is highly confident about what belongs here.
- **Hydrophobic core position:** Moderate entropy, concentrated on hydrophobic residues. Example: P(Leu) = 0.25, P(Ile) = 0.22, P(Val) = 0.20, P(Phe) = 0.12. The model knows this must be hydrophobic but allows variation.
- **Solvent-exposed loop:** High entropy, nearly uniform. Each amino acid gets roughly 0.04–0.07 probability. The model recognizes this position is unconstrained.
- **Structural glycine (tight turn):** Sharp distribution strongly favoring glycine. P(Gly) = 0.90. The model knows only glycine fits the structural constraint.

### 7.2 JS Divergence as a Similarity Measure

The Jensen-Shannon divergence between two probability distributions P and Q is defined as:

```
JS(P, Q) = ½ KL(P || M) + ½ KL(Q || M)    where M = ½(P + Q)
```

Key properties that make JS ideal for this application:

- **Bounded [0, 1]:** Scores are on a known scale, simplifying the combination with BLOSUM.
- **Symmetric:** JS(P, Q) = JS(Q, P), which is required for alignment scoring.
- **Information-theoretically principled:** Measures how much information you lose by representing both distributions with their mixture. This is a natural measure of evolutionary constraint similarity.
- **Efficient to compute:** For 20-dimensional distributions, requires ~60 multiplications and ~40 logarithms per cell.

### 7.3 The Adaptive Weighting in Practice

The information-content weighting creates a natural adaptive behavior across the alignment:

| Scenario | BLOSUM | MLM | Effective Behavior |
|---|---|---|---|
| Close homologs, conserved site | Dominant | Reinforcing | Both agree; high confidence alignment |
| Close homologs, variable site | Dominant | Weak (low IC) | BLOSUM leads; MLM contributes little |
| Distant homologs, conserved site | Weak (AAs differ) | **Dominant** | MLM rescues alignment where BLOSUM fails |
| Distant homologs, variable site | Weak | Weak | Low confidence; handled by gap model |

The third row is the critical case: **distant homologs at conserved sites**. This is precisely the twilight zone scenario where MAFFT and all traditional tools struggle, and where our method provides its greatest advantage. The MLM distributions at two equivalent constrained positions will be similar even when the actual amino acids are completely different, because the evolutionary constraints (encoded in the language model from training on millions of sequences) are the same.

---

## 8. Models, Software, and Resources

### 8.1 Recommended Language Model: ESM-2

ESM-2 is the recommended protein language model for this project due to its **MIT license** (fully open for commercial and non-commercial use), strong performance, and range of model sizes. Both the original GitHub repository and the HuggingFace weights carry the MIT license.

| Model | Params | VRAM | Speed | Quality | Recommended For |
|---|---|---|---|---|---|
| `esm2_t6_8M` | 8M | CPU OK | ~ms/seq | Decent | Testing only |
| `esm2_t12_35M` | 35M | ~2 GB | ~ms/seq | Good | CPU deployment |
| `esm2_t30_150M` | 150M | ~4 GB | ~10ms/seq | Very good | **Default choice** |
| **`esm2_t33_650M`** | **650M** | **~8 GB** | **~50ms/seq** | **Excellent** | **Best accuracy** |
| `esm2_t36_3B` | 3B | ~16 GB | ~200ms/seq | Best | If GPU budget allows |

The **650M model (`esm2_t33_650M_UR50D`)** is recommended as the primary model for benchmarking and publication, with the 150M model as the default for production use.

### 8.2 Where to Obtain ESM-2

- **Official GitHub Repository:** https://github.com/facebookresearch/esm — Contains all model code, pretrained weights (auto-downloaded), and usage examples. Install with: `pip install fair-esm`
- **HuggingFace Hub:** https://huggingface.co/facebook/esm2_t33_650M_UR50D (and other sizes). Use with the HuggingFace Transformers library: `pip install transformers torch`
- **PyTorch Hub:** Models can also be loaded via `torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")`

### 8.3 Quick Start Code: Extracting MLM Distributions

```python
import torch
import esm

# Load model (weights auto-download on first use)
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()

# Prepare sequences
data = [("protein1", "MKTLLILAVL..."), ("protein2", "MKGFVLIAL...")]
batch_labels, batch_strs, batch_tokens = batch_converter(data)

# Single forward pass — pseudo-likelihood approximation
with torch.no_grad():
    logits = model(batch_tokens)["logits"]  # (batch, seq_len, vocab_size)

# Extract 20 standard amino acid probabilities
# ESM-2 vocab includes special tokens; extract only the 20 standard AAs
aa_indices = [alphabet.get_idx(aa) for aa in "ACDEFGHIKLMNPQRSTVWY"]
mlm_logits = logits[:, :, aa_indices]  # (batch, seq_len, 20)
mlm_probs = torch.softmax(mlm_logits, dim=-1)  # probability distributions

# mlm_probs[0] is now (seq_len, 20) — the MLM distribution for each position
```

### 8.4 Quick Start Code: Computing the Hybrid Score

```python
import numpy as np
from scipy.spatial.distance import jensenshannon
from Bio.Align import substitution_matrices

blosum62 = substitution_matrices.load("BLOSUM62")
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"

def information_content(prob_dist):
    """IC = log2(20) - entropy"""
    entropy = -np.sum(prob_dist * np.log2(prob_dist + 1e-10))
    return np.log2(20) - entropy

def hybrid_score(seq1, seq2, mlm1, mlm2, alpha=1.0, gamma=1.5):
    """
    Compute the hybrid scoring matrix for two sequences.

    Args:
        seq1, seq2: amino acid sequences (strings)
        mlm1, mlm2: MLM probability distributions (L1×20, L2×20 numpy arrays)
        alpha: BLOSUM weight
        gamma: MLM weight

    Returns:
        score_matrix: L1 × L2 numpy array
        gap_penalties_1, gap_penalties_2: position-specific gap open penalties
    """
    L1, L2 = len(seq1), len(seq2)
    score_matrix = np.zeros((L1, L2))

    # Precompute information content
    ic1 = np.array([information_content(mlm1[i]) for i in range(L1)])
    ic2 = np.array([information_content(mlm2[j]) for j in range(L2)])

    for i in range(L1):
        for j in range(L2):
            # Component 1: BLOSUM62
            blosum_score = blosum62[seq1[i]][seq2[j]]

            # Component 2: MLM similarity (1 - JS divergence) weighted by IC
            js_div = jensenshannon(mlm1[i], mlm2[j]) ** 2  # squared = divergence
            mlm_sim = (1.0 - js_div)
            confidence = (ic1[i] + ic2[j]) / 2.0

            # Combined score
            score_matrix[i, j] = alpha * blosum_score + gamma * mlm_sim * confidence

    # Position-specific gap penalties from MLM entropy
    base_gap_open = -11.0
    lambda_gap = 0.3
    gap_penalties_1 = base_gap_open * (1 + lambda_gap * ic1)
    gap_penalties_2 = base_gap_open * (1 + lambda_gap * ic2)

    return score_matrix, gap_penalties_1, gap_penalties_2
```

### 8.5 Alternative Models

| Model | Source | Advantage | Consideration |
|---|---|---|---|
| **ProtT5-XL** | HuggingFace: `Rostlab/prot_t5_xl_uniref50` | T5 architecture; strong embeddings; proven in learnMSA2 | Larger; encoder-decoder architecture (use encoder only) |
| **ProstT5** | HuggingFace: `Rostlab/ProstT5` | Can translate between sequence and structure tokens | Interesting for structure-aware gap penalties |
| **ANKH** | HuggingFace: `ElnaggarLab/ankh-base` | Optimized for efficiency; smaller models available | Less community adoption; fewer benchmarks |
| **ESM-C** | EvolutionaryScale (check licensing) | Newer model; potentially better predictions | Licensing may restrict use; verify before adoption |

### 8.6 Licensing Comparison

| Model | License | C/C++ Inference | Notes |
|---|---|---|---|
| **ESM-2** | **MIT** | **Via ONNX Runtime C API** | Export once in Python, run pure C at runtime |
| ESM-C | Cambrian (custom) | No | Python only, restrictive license |
| ProtTrans/ProtBERT | MIT | Also ONNX-exportable | Older, lower quality than ESM-2 |

ESM-2's MIT license is a significant practical advantage: there are no restrictions on redistribution, modification, or commercial use. This simplifies packaging, distribution, and integration into existing open-source tools like Kalign.

### 8.7 Running ESM-2 from C via ONNX Runtime

There is no native C implementation of ESM-2, but the model can be run from pure C/C++ at runtime via ONNX Runtime, eliminating any Python dependency in the final deployed tool. This is the recommended path for production integration with Kalign.

**Step 1: One-time ONNX export (Python)**

Convert the PyTorch model to ONNX format. This is done once during development, not at runtime:

- **Dedicated tool:** https://github.com/ashrafgt/esm-2-onnx — a project specifically for converting ESM-2 to ONNX with embedding and logit extraction support.
- **Alternative:** HuggingFace's `optimum` library provides generic ONNX export: `optimum-cli export onnx --model facebook/esm2_t33_650M_UR50D esm2_650M_onnx/`

The exported `.onnx` file contains the full model graph and weights in a single portable file.

**Step 2: C/C++ inference via ONNX Runtime**

ONNX Runtime (https://onnxruntime.ai/) provides a full C API (`onnxruntime_c_api.h`) for loading and running ONNX models:

- Load the `.onnx` file, feed tokenized sequences as input tensors, get logits back
- **No Python dependency at runtime** — the binary links only against the ONNX Runtime shared library
- Supports CPU, CUDA, CoreML/Metal, TensorRT, and other execution providers
- Documentation: https://onnxruntime.ai/docs/get-started/with-cpp.html

```c
// Pseudocode: ESM-2 inference from C via ONNX Runtime
#include <onnxruntime_c_api.h>

// One-time setup
OrtEnv* env;
OrtCreateEnv(ORT_LOGGING_LEVEL_WARNING, "esm2", &env);
OrtSession* session;
OrtCreateSession(env, "esm2_650M.onnx", opts, &session);

// Per-sequence inference
int64_t token_ids[] = {0, 15, 7, 4, ...};  // tokenized sequence + BOS/EOS
OrtValue* input_tensor;   // create from token_ids
OrtValue* output_tensor;  // logits: (1, seq_len, vocab_size)
OrtRun(session, NULL, input_names, &input_tensor, 1, output_names, 1, &output_tensor);

// Extract 20 standard AA logits → softmax → MLM probability distributions
float* logits = OrtGetTensorData(output_tensor);
// softmax over the 20 AA indices, store as mlm_probs[pos][20]
```

**Deployment architecture:**

The recommended production deployment is a single statically-linked binary:

```
┌─────────────────────────────────────────────────────┐
│              kalign-mlm (single binary)              │
├──────────────────────┬──────────────────────────────┤
│   Kalign C engine    │   ONNX Runtime C library     │
│  (alignment, DP,     │  (loads esm2_*.onnx,         │
│   guide tree, I/O)   │   runs inference, CPU/GPU)   │
├──────────────────────┴──────────────────────────────┤
│           Hybrid scoring function (C)                │
│  BLOSUM62 lookup + JS divergence on 20d MLM dists   │
│  + IC weighting + position-specific gap penalties    │
└─────────────────────────────────────────────────────┘
```

No Python is required at any point in the runtime pipeline. The only Python involvement is the one-time ONNX export during development, and the pre-exported `.onnx` model files can be distributed alongside the binary.

### 8.8 Benchmark Datasets

- **BAliBASE 4.0:** https://lbgi.fr/balibase/ — The gold standard for MSA benchmarking. Contains reference alignments at various divergence levels. Essential for parameter tuning (α, γ, λ).
- **SABmark 1.65:** http://bioinformatics.vub.ac.be/databases/SABmark/ — Structural alignments for superfamily and twilight zone evaluation.
- **PREFAB 4.0:** Reference alignments from structural superposition. Available from the MUSCLE repository.
- **HomFam:** Large-scale families for testing scalability. Available from the Clustal group.
- **QuanTest 2.0:** Recent benchmark with quantitative scoring. Used by ARIES and learnMSA2.

### 8.9 Key Software Dependencies

- **Kalign source:** https://github.com/TimoLassmann/kalign — The target alignment engine. Written in C, fast, SIMD-optimized. The scoring function modifications would be implemented here.
- **ONNX Runtime:** https://onnxruntime.ai/ — C API for running the exported ESM-2 model. Required for production integration. Supports CPU and GPU backends. Install C library from GitHub releases or build from source.
- **esm-2-onnx:** https://github.com/ashrafgt/esm-2-onnx — One-time Python tool for exporting ESM-2 to ONNX format.
- **PyTorch:** https://pytorch.org/ — Required for prototyping (Phases 1–2) and for the one-time ONNX export step. Not needed at runtime in production.
- **SciPy:** For JS divergence computation (`scipy.spatial.distance.jensenshannon`) during prototyping.
- **NumPy:** For efficient array operations on probability distributions.
- **Biopython:** For BLOSUM matrix loading and sequence I/O (`pip install biopython`).

---

## 9. Implementation Roadmap

### 9.1 Phase 1: Proof of Concept (2–4 weeks)

- Implement the three-component scoring function in Python.
- Run on 10–20 hand-picked twilight-zone pairs from BAliBASE.
- Visualize: BLOSUM-only score matrix vs. hybrid score matrix vs. reference alignment.
- Validate that MLM distributions carry a measurably different signal from embeddings.
- **Deliverable:** Jupyter notebook with compelling visual evidence.

### 9.2 Phase 2: Pairwise Aligner (4–6 weeks)

- Build a complete Needleman-Wunsch aligner using the hybrid scoring function.
- Implement position-specific gap penalties from MLM entropy.
- Benchmark on BAliBASE pairwise subsets against BLOSUM-only, PLMAlign, and EBA.
- Tune α, γ, λ on a training split.
- **Deliverable:** Pairwise aligner with quantitative comparison to existing methods.

### 9.3 Phase 3: Integration into Kalign (6–10 weeks)

- Export ESM-2 model(s) to ONNX format using esm-2-onnx or HuggingFace optimum.
- Integrate ONNX Runtime C API into Kalign for native ESM-2 inference — no Python dependency at runtime.
- Implement the hybrid scoring function, JS divergence, IC weighting, and position-specific gap penalties in C.
- Modify Kalign's alignment loop to consume the 20-dimensional MLM distributions from ONNX Runtime output.
- Add an adaptive mode flag: `--mlm` for MLM-augmented scoring (default for <40% identity).
- Full MSA benchmarking on BAliBASE, SABmark, PREFAB, HomFam, QuanTest 2.0.
- **Deliverable:** Single self-contained binary (kalign-mlm) with MLM mode, ready for publication. Ships with pre-exported `.onnx` model files.

### 9.4 Phase 4: Publication and Release (4–6 weeks)

- Write manuscript targeting Bioinformatics, NAR, or ISMB proceedings.
- Release open-source code: single C binary + pre-exported ONNX model files.
- Provide Docker image with ONNX Runtime and model files pre-loaded (no Python required).
- Publish pip-installable Python package for prototyping/scripting use cases.
- Submit to bioRxiv for community feedback before journal submission.

---

## 10. Expected Outcomes and Success Criteria

The project will be considered successful if:

1. **Accuracy on twilight-zone alignments:** Our method achieves higher SP (sum-of-pairs) and TC (total column) scores than MAFFT L-INS-i on BAliBASE references with <30% sequence identity.
2. **Competitive speed:** Total alignment time (including ESM-2 inference) is within 5× of standard Kalign for typical inputs (<1000 sequences).
3. **MLM signal validation:** Ablation study demonstrates that removing the MLM component degrades accuracy, confirming it provides information beyond BLOSUM alone.
4. **Gap penalty validation:** Position-specific gap penalties from MLM entropy improve alignment quality compared to uniform gap penalties.
5. **Complementarity to embeddings:** Analysis shows cases where MLM distributions outperform embedding similarity and vice versa, confirming they capture different signals.

---

## 11. Risks and Mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| Pseudo-likelihood approximation is insufficiently accurate | MLM distributions do not match true masked marginals, reducing scoring quality | Empirically validated for mutation scoring; can fall back to true masking for a subset of positions if needed |
| Marginal gain over BLOSUM for close homologs | Method only helps in twilight zone, limiting impact | Expected and acceptable. The adaptive mode automatically uses BLOSUM-only for high-identity pairs, so no regression |
| Reviewers consider it incremental | Paper rejected as "just another PLM alignment method" | Emphasize the MLM vs. embedding distinction, the 60× efficiency gain, gap penalties, and information-content weighting as novel contributions |
| ESM-2 inference is too slow for production | Method is impractical for large-scale use | Use 35M model for CPU-only environments; ONNX Runtime provides optimized CPU/GPU inference with no Python overhead; precompute and cache distributions; use only for <40% identity pairs |

---

## 12. Key References

**Language Models:**

- Lin, Z. et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* 379(6637).
- Rives, A. et al. (2021). Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences. *PNAS* 118(15).
- Elnaggar, A. et al. (2022). ProtTrans: Toward understanding the language of life through self-supervised learning. *IEEE TPAMI* 44(10).

**PLM-Based Alignment:**

- Liu, W. et al. (2024). PLMSearch: Protein language model powers accurate and fast sequence search for remote homology. *Nat Commun* 15, 2775.
- Pantolini, L. et al. (2024). Embedding-based alignment: combining protein language models with DP alignment. *Bioinformatics* 40(1).
- Becker, F. & Stanke, M. (2024). learnMSA2: deep protein multiple alignments with large language and hidden Markov models. *Bioinformatics* 40(S2).
- McWhite, C. & Singh, M. (2023). Leveraging protein language models for accurate multiple sequence alignments. *Genome Research* 33(7).

**Alignment Benchmarking and Methodology:**

- Katoh, K. & Standley, D.M. (2013). MAFFT multiple sequence alignment software version 7. *Mol Biol Evol* 30(4).
- Lassmann, T. (2020). Kalign 3: multiple sequence alignment of large datasets. *Bioinformatics* 36(6).
- Marquet, C. et al. (2022). Embeddings from protein language models predict conservation and variant effects. *Human Genetics*.

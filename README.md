[![CMake](https://github.com/TimoLassmann/kalign/actions/workflows/cmake.yml/badge.svg)](https://github.com/TimoLassmann/kalign/actions/workflows/cmake.yml)
[![Python](https://github.com/TimoLassmann/kalign/actions/workflows/python.yml/badge.svg)](https://github.com/TimoLassmann/kalign/actions/workflows/python.yml)
[![Build Python Wheels](https://github.com/TimoLassmann/kalign/actions/workflows/wheels.yml/badge.svg)](https://github.com/TimoLassmann/kalign/actions/workflows/wheels.yml)
![CodeQL](https://github.com/TimoLassmann/kalign/workflows/CodeQL/badge.svg)

# Kalign

Kalign is a fast multiple sequence alignment program for biological
sequences. It aligns protein, DNA, and RNA sequences using a progressive
alignment approach with multi-threading support.

## Installation

### From source

Prerequisites: C compiler (GCC or Clang), CMake 3.18+.

```bash
mkdir build && cd build
cmake ..
make
make test
make install
```

Kalign runs on Linux and macOS. It parallelises with a built-in thread pool (pthreads); with `-n 1` it runs serially.

On x86-64, Kalign is built with AVX2 when the compiler and the build machine support it and does not check the CPU at run time. The x86-64 Python wheels need a CPU with AVX2 (Intel Haswell, 2013, or later; AMD Excavator, 2015, or later). On older CPUs, build from source with `cmake -DENABLE_AVX2=OFF ..`.

### Zig build (alternative)

Tested with zig 0.17.

```bash
zig build
```

### Python

```bash
pip install kalign-python
```

See [README-python.md](README-python.md) for the full Python documentation.

## Usage

```
kalign -i <input> -o <output>
```

Kalign 4 has two modes, each with a protein and a DNA/RNA setting; the
sequence type is detected from the input, or set with `--type protein`,
`--type dna` or `--type rna` (DNA and RNA use the same setting):

| Mode | Flag | Description |
|------|------|-------------|
| fast | `--mode fast` | Single run, fastest. This is the default (`--mode default` is an alias). |
| accurate | `--mode accurate` | 5-member ensemble with per-column confidence, slower. At most 5,000 sequences (`--mode recall` is an alias). |

The four settings were fitted on reference alignments inside grouped nested
cross-validation: fast by multi-objective search (NSGA-III); accurate by an
NSGA-III ensemble search followed by a small grid over guide-tree noise and
consensus support. Each shipped setting is the final fit on all training data,
then scored on held-out data that the fit did not use. The fit code, the
records and the settings are in [presets/](presets/).

### Examples

```bash
# Align sequences (default mode)
kalign -i sequences.fa -o aligned.fa

# Fast mode
kalign --mode fast -i sequences.fa -o aligned.fa

# Accurate mode (ensemble)
kalign --mode accurate -i sequences.fa -o aligned.fa

# Read from stdin
cat input.fa | kalign -i - -o aligned.fa

# Combine multiple input files
kalign seqsA.fa seqsB.fa -o combined.fa
```

### Options

Options that change a single alignment (gap penalties and the expert options)
only work with `--mode fast`. Ensemble options (`--fill-*`) only work with
`--mode accurate`. Using one in the wrong mode is an error.

General:

```
--mode         Mode: fast (default), accurate; "default" = fast, "recall" = accurate.
--type         Sequence type: protein, dna, rna (overrides detection). [auto]
--format, -f   Output format: fasta, msf, clu. [fasta]
-n             Number of threads (integer >= 1). [cores - 1, at most 16]
-q             Quiet: suppress messages.
--version      Print the version.
```

Fast mode only:

```
--gpo          Gap open penalty (overrides the fitted setting).
--gpe          Gap extension penalty (overrides the fitted setting).
--tgpe         Terminal gap extension penalty (overrides the fitted setting).
```

Expert options, also fast mode only: `--vsm-mode`, `--vsm-kappa`,
`--vsm-dzero` and `--guide-dist`. See `kalign` without arguments for details.

Accurate mode only:

```
--fill-in 0-3        Place unaligned consensus residues (1 realign, 2 project, 3 attach).
--fill-anchor-frac F Fill-in anchor needs at least F*numseq residues (minimum 2).
--fill-support N     Attach fill-in: votes needed to join a column.
```

Confidence (used by `--mode accurate`; in fast mode they have no effect and
Kalign warns):

```
--confidence-threshold F  Mask columns below this confidence (0-1). [off]
--confidence-style   Masking style: lowercase or remove. [lowercase]
--confidence-output  Write per-column confidence to a file. [off]
```

Adding sequences to an existing alignment (cannot be combined with `--mode`,
the fast-only options or the accurate-only options):

```
--add FILE       New sequences to add.
--existing FILE  Existing alignment the new sequences are added to.
```

The command line uses (cores - 1) threads by default, at most 16. The Python
package and `kalign-py` default to 1 thread; ask for more with
`kalign.align(..., n_threads=8)`, `kalign.set_num_threads(8)` or `kalign-py -n 8`.

### Output formats

```bash
kalign -i input.fa -f msf -o output.msf
kalign -i input.fa -f clu -o output.clu
```

## C library

Link Kalign into your C/C++ project:

```cmake
find_package(kalign)
target_link_libraries(<target> kalign::kalign)
```

Or include directly:

```cmake
add_subdirectory(<path>/kalign EXCLUDE_FROM_ALL)
target_link_libraries(<target> kalign::kalign)
```

`kalign()` aligns an array of sequences with the fitted fast setting for the
detected (or given) type, as the command line and Python do. A negative
gap penalty keeps the fitted value.

```c
#include <kalign/kalign.h>

char *seq[] = {"ATCGATCGAT", "ATCGTCGAT", "ATCGATCGGAT"};
int len[] = {10, 9, 11};
char **aligned = NULL;
int aln_len = 0;
/* 1 thread, type detected, fitted gap penalties */
kalign(seq, len, 3, 1, KALIGN_TYPE_UNDEFINED, -1.0f, -1.0f, -1.0f, &aligned, &aln_len);
/* aligned[0..2] are aln_len characters each; free each, then aligned */
```

To run a mode on a sequence file, as the command line does:

```c
struct msa *msa = NULL;
kalign_read_input("in.fa", &msa, 1);
kalign_align_mode(msa, "accurate", KALIGN_TYPE_UNDEFINED, 4);
kalign_write_msa(msa, "out.fa", "fasta");
kalign_free_msa(msa);
```

## Benchmarks

The Kalign 4 settings and their benchmark against FAMSA 2, MAFFT, MUSCLE 5 and
Clustal Omega (protein and DNA/RNA reference alignments, simulations, alignment
size up to 50,000 sequences, and trees, positive selection and homology search
built on the alignments) are described in a preprint in preparation. The
benchmarks of earlier versions are in the ChangeLog.

## Citation

Lassmann, Timo. "Kalign 3: multiple sequence alignment of large data sets."
Bioinformatics (2019). [DOI](https://doi.org/10.1093/bioinformatics/btz795)

## License

Apache License, Version 2.0. See [COPYING](COPYING).

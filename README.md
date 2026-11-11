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

### Zig build (alternative)

Tested with zig 0.16.

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

```
--mode         Mode preset: fast (default), accurate; "default" = fast, "recall" = accurate. [fast]
--format       Output format: fasta, msf, clu. [fasta]
--type         Sequence type: protein, dna, rna (overrides detection). [auto]
--gpo          Gap open penalty (overrides preset). [auto]
--gpe          Gap extension penalty (overrides preset). [auto]
--tgpe         Terminal gap extension penalty (overrides preset). [auto]
-n             Number of threads. [auto]
```

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

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

Prerequisites: C compiler (GCC or Clang), CMake 3.18+, optionally OpenMP.

```bash
mkdir build && cd build
cmake ..
make
make test
make install
```

On macOS, `brew install libomp` for OpenMP support.

### Zig build (alternative)

Requires zig version 0.12.

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

Kalign v3.5 has three modes:

| Mode | Flag | Description |
|------|------|-------------|
| default | (none) | Best general-purpose. |
| fast | `--fast` | Fastest. Same as kalign v3.4. |
| precise | `--precise` | Highest accuracy, ~10x slower. |

### Examples

```bash
# Align sequences
kalign -i sequences.fa -o aligned.fa

# Fast mode
kalign --fast -i sequences.fa -o aligned.fa

# Precise mode (ensemble + realign)
kalign --precise -i sequences.fa -o aligned.fa

# Read from stdin
cat input.fa | kalign -i - -o aligned.fa

# Combine multiple input files
kalign seqsA.fa seqsB.fa -o combined.fa

# Save ensemble consensus for re-thresholding
kalign --precise -i seqs.fa -o out.fa --save-poar consensus.poar
kalign -i seqs.fa -o out2.fa --load-poar consensus.poar --min-support 3
```

### Options

```
--format       Output format: fasta, msf, clu. [fasta]
--type         Sequence type: protein, dna, rna, divergent. [auto]
--gpo          Gap open penalty. [auto]
--gpe          Gap extension penalty. [auto]
--tgpe         Terminal gap extension penalty. [auto]
--ensemble N   Run N ensemble alignments. [off]
--refine       Refinement: none, all, confident. [none]
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

### Balibase
![Balibase_scores](https://user-images.githubusercontent.com/8110320/198513840-0e08a634-bb41-4826-bd58-7fc66eae1054.jpeg)

### Bralibase
![Bralibase_scores](https://user-images.githubusercontent.com/8110320/198513850-00e5037f-355f-45ec-828f-ed8d47497272.jpeg)

## Citation

Lassmann, Timo. "Kalign 3: multiple sequence alignment of large data sets."
Bioinformatics (2019). [DOI](https://doi.org/10.1093/bioinformatics/btz795)

## License

Apache License, Version 2.0. See [COPYING](COPYING).

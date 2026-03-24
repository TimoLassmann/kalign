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

Kalign uses a built-in thread pool for parallelization (requires pthreads, available on all POSIX systems). If pthreads is not available, it falls back to serial execution. To use OpenMP instead:

```bash
cmake -DUSE_OPENMP=ON -DUSE_THREADPOOL=OFF ..
```

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

Kalign has four mode presets, optimized for protein and nucleotide sequences:

| Mode | Flag | Description |
|------|------|-------------|
| fast | `--mode fast` | Single run, fastest. |
| default | `--mode default` | Single run with consistency anchors (default). |
| recall | `--mode recall` | Ensemble, optimized for recall. |
| accurate | `--mode accurate` | Ensemble, highest precision. |

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
--mode         Mode preset: fast, default, recall, accurate. [default]
--format       Output format: fasta, msf, clu. [fasta]
--type         Sequence type: protein, dna, rna, divergent. [auto]
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

### Balibase
![Balibase_scores](https://user-images.githubusercontent.com/8110320/198513840-0e08a634-bb41-4826-bd58-7fc66eae1054.jpeg)

### Bralibase
![Bralibase_scores](https://user-images.githubusercontent.com/8110320/198513850-00e5037f-355f-45ec-828f-ed8d47497272.jpeg)

## Citation

Lassmann, Timo. "Kalign 3: multiple sequence alignment of large data sets."
Bioinformatics (2019). [DOI](https://doi.org/10.1093/bioinformatics/btz795)

## License

Apache License, Version 2.0. See [COPYING](COPYING).

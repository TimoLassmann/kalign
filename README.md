<!-- ![C/C++ CI](https://github.com/TimoLassmann/kalign/workflows/C/C++%20CI/badge.svg) -->
[![CMake](https://github.com/TimoLassmann/kalign/actions/workflows/cmake.yml/badge.svg)](https://github.com/TimoLassmann/kalign/actions/workflows/cmake.yml)
![CodeQL](https://github.com/TimoLassmann/kalign/workflows/CodeQL/badge.svg)

# Kalign

Kalign is a fast multiple sequence alignment program for biological sequences.

# Installation

## Release Tarball

Download tarball from [releases](https://github.com/TimoLassmann/kalign/releases). Then:

``` bash
tar -zxvf kalign-<version>.tar.gz
cd kalign-<version>
mkdir build 
cd build
cmake .. 
make 
make test 
make install
```

on macOS, install [brew](https://brew.sh/) then:

``` bash
brew install cmake 
git clone https://github.com/TimoLassmann/kalign.git
cd kalign
mkdir build
cd build 
cmake ..
make 
make test 
make install
```

# Usage

The command line interface of Kalign accepts the following options:

```bash
Usage: kalign  -i <seq file> -o <out aln> 

Options:

   --format           : Output format. [Fasta]
   --type             : Alignment type (rna, dna, internal). [rna]
                        Options: protein, divergent (protein) 
                                 rna, dna, internal (nuc). 
   --gpo              : Gap open penalty. []
   --gpe              : Gap extension penalty. []
   --tgpe             : Terminal gap extension penalty. []
   -n/--nthreads      : Number of threads. [4]
   --version (-V/-v)  : Prints version. [NA]


```


Kalign expects the input to be a set of unaligned sequences in fasta format or aligned sequences in aligned fasta, MSF or clustal format. If the sequences are already aligned, kalign will remove all gap characters and re-align the sequences. 

By default, Kalign automatically detects whether the input sequences are protein or DNA and selects appropriate alignment parameters. 

The `--type` option gives users more direct control over the alignment parameters. Currently there are five core options:

- `protein`  : uses a the CorBLOSUM66_13plus substituion matrix (default for protein sequence)
- `divergent`: uses the gonnet 250 substituion matrix 
- `dna`      : default DNA parameters
  +  5 match score 
  + -4 mismatch score
  + -8 gap open penalty
  + -6 gap extension penalty 
  +  0 terminal gap extension penalty
- `internal` : same as above but terminal gaps set to 8 to encourage gaps within the sequences. 
- `rna`      : parameters optimised for RNA alignments.

The `--gpo`, `--gpe` and `--tgpe` options can be used to further fine tune the parameters.

# Examples

Passing sequences via stdin:

```bash
cat input.fa | kalign -f fasta > out.afa
```

Combining multiple input files:

```bash
kalign seqsA.fa seqsB.fa seqsC.fa -f fasta > combined.afa
```

Align sequences and output the alignment in MSF format:

```bash
kalign -i BB11001.tfa -f msf  -o out.msf
```

Align sequences and output the alignment in clustal format:

```bash
kalign -i BB11001.tfa -f clu -o out.clu
```

Re-align sequences in an existing alignment:

```bash
kalign -i BB11001.msf  -o out.afa
```

Reformat existing alignment:

```bash
kalign -i BB11001.msf -r afa -o out.afa
```

# Kalign library 

To incorporate Kalign into your own projects you can link to the library like this: 

```cmake 
find_package(kalign)
target_link_libraries(<target> kalign::kalign)
```

Alternatively, you can include the kalign code directly in your project and link with:

```cmake
if (NOT TARGET kalign)
  add_subdirectory(<path_to_kalign>/kalign EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(<target> kalign::kalign)
```
# Benchmark results

Here are some benchmark results. The code to reproduce these figures can be found at [here](scripts/benchmark.org).

## Balibase

![Balibase_scores](https://user-images.githubusercontent.com/8110320/198513840-0e08a634-bb41-4826-bd58-7fc66eae1054.jpeg)

## Bralibase

![Bralibase_scores](https://user-images.githubusercontent.com/8110320/198513850-00e5037f-355f-45ec-828f-ed8d47497272.jpeg)

# Please cite:
1. Lassmann, Timo. _Kalign 3: multiple sequence alignment of large data sets._ **Bioinformatics** (2019). [pdf](https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btz795/30314127/btz795.pdf)

# Other papers:
1. Lassmann, Timo, Oliver Frings, and Erik LL Sonnhammer. _Kalign2: high-performance multiple alignment of protein and nucleotide sequences allowing external features._ **Nucleic acids research** 37.3 (2008): 858-865. [Pubmed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647288/)
2. Lassmann, Timo, and Erik LL Sonnhammer. _Kalign: an accurate and fast multiple sequence alignment algorithm._ **BMC bioinformatics** 6.1 (2005): 298. [Pubmed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1325270/)

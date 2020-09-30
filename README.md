![C/C++ CI](https://github.com/TimoLassmann/kalign/workflows/C/C++%20CI/badge.svg)
![CodeQL](https://github.com/TimoLassmann/kalign/workflows/CodeQL/badge.svg)

# Kalign

Kalign is a fast multiple sequence alignment program for biological sequences.

# Installation

## Release Tarball

Download tarball from [releases](https://github.com/TimoLassmann/kalign/releases). Then:

``` bash
tar -zxvf kalign-<version>.tar.gz
cd kalign-<version>
./autogen.sh
./configure
make
make check
make install
```

## Homebrew
``` bash
brew install brewsci/bio/kalign
```

## Developer version
``` bash
git clone https://github.com/TimoLassmann/kalign.git
cd kalign
./autogen.sh
./configure
make
make check
make install
```

on macOS, install [brew](https://brew.sh/) then:

``` bash
brew install libtool
brew install automake
git clone https://github.com/TimoLassmann/kalign.git
cd kalign
./autogen.sh
./configure
make
make check
make install
```

# Usage


``` bash
Usage: kalign  -i <seq file> -o <out aln>

Options:

   --format           : Output format. [Fasta]
   --reformat         : Reformat existing alignment. [NA]
   --version          : Print version and exit
```

Kalign expects the input to be a set of unaligned sequences in fasta format or aligned sequences in aligned fasta, MSF or clustal format. Kalign automatically detects whether the input sequences are protein, RNA or DNA.

Since version 3.2.0 kalign supports passing sequence in via stdin and support alignment of sequences from multiple files.

# Examples

Passing sequences via stdin:

```
cat input.fa | kalign -f fasta > out.afa
```

Combining multiple input files:

```
kalign seqsA.fa seqsB.fa seqsC.fa -f fasta > combined.afa
```

Align sequences and output the alignment in MSF format:

```
kalign -i BB11001.tfa -f msf  -o out.msf
```

Align sequences and output the alignment in clustal format:

```
kalign -i BB11001.tfa -f clu -o out.clu
```

Re-align sequences in an existing alignment:

```
kalign -i BB11001.msf  -o out.afa
```

Reformat existing alignment:

```
kalign -i BB11001.msf -r afa -o out.afa
```

# Benchmark results

Here are some benchmark results. The code to reproduce these figures can be found at [here](scripts/benchmark.org).

## Balibase

![Balibase_scores](https://user-images.githubusercontent.com/8110320/66697423-7ea3d000-eca3-11e9-919a-995ca8e9f7c1.jpeg)

## Bralibase

![Bralibase_scores](https://user-images.githubusercontent.com/8110320/66697424-86637480-eca3-11e9-90ea-238f82b0ac6b.jpeg)

## Homfam

![Homfam_scores](https://user-images.githubusercontent.com/8110320/66697425-895e6500-eca3-11e9-97e7-63f3a79133cf.jpeg)

## Quantest2

![Quantest2_scores](https://user-images.githubusercontent.com/8110320/66698153-6c2c9500-eca9-11e9-904c-3d6ea9a1c44d.jpeg)

# Please cite:
1. Lassmann, Timo. _Kalign 3: multiple sequence alignment of large data sets._ **Bioinformatics** (2019). [pdf](https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btz795/30314127/btz795.pdf)

# Other papers:
1. Lassmann, Timo, Oliver Frings, and Erik LL Sonnhammer. _Kalign2: high-performance multiple alignment of protein and nucleotide sequences allowing external features._ **Nucleic acids research** 37.3 (2008): 858-865. [Pubmed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647288/)
2. Lassmann, Timo, and Erik LL Sonnhammer. _Kalign: an accurate and fast multiple sequence alignment algorithm._ **BMC bioinformatics** 6.1 (2005): 298. [Pubmed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1325270/)

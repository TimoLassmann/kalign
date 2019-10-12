# Kalign

Kalign is a fast multiple sequence alignment program for biological sequences.

![Image of example output](https://user-images.githubusercontent.com/8110320/62088549-330d8880-b255-11e9-928d-e5cc8031da97.png)

# Install

``` bash
git clone https://github.com/TimoLassmann/kalign.git 
cd kalign
./autogen.sh
make 
make check 
make install 
```

# Usage


``` sh
Usage: kalign  -i <seq file> -o <out aln> 

Options:

   --format           : Output format. [Fasta]
   --reformat         : Reformat existing alignment. [NA]
```

Kalign expects the input to be a set of unaligned sequences in fasta format or aligned sequences in aligned fasta, MSF or clustal format. Kalign automatically detects whether the input sequences are proten, RNA or DNA.

# Examples

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

## Balibase
![Optional Text](../master/doc/images/Balibase_scores.png)
## Bralibase

![Optional Text](../master/myFolder/image.png)
## Homfam

![Optional Text](../master/myFolder/image.png)

# Please cite:

Lassmann, Timo, Oliver Frings, and Erik LL Sonnhammer.
Kalign2: high-performance multiple alignment of protein and
nucleotide sequences allowing external features.
Nucleic acids research 37.3 (2008): 858-865.
        
Lassmann, Timo, and Erik LL Sonnhammer. Kalign–an accurate and
fast multiple sequence alignment algorithm.

BMC bioinformatics 6.1 (2005): 298.


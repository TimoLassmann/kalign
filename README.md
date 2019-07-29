# Kalign

Kalign is a fast multiple alignment program for protein and nucleotide sequences.


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

kalign -i <input sequences> -o <output alignment file> 
```

# Please cite:

Lassmann, Timo, Oliver Frings, and Erik LL Sonnhammer.
Kalign2: high-performance multiple alignment of protein and
nucleotide sequences allowing external features.
Nucleic acids research 37.3 (2008): 858-865.
        
Lassmann, Timo, and Erik LL Sonnhammer. Kalign–an accurate and
fast multiple sequence alignment algorithm.
BMC bioinformatics 6.1 (2005): 298.


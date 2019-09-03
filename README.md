# ReSequenceR
More realistic simulator for genomic DNA sequences from Illumina machines that achieves a similar k-mer spectrum as the original sequences

## Requirements
Compiler supporting C++14 (Tested on gcc version 7.2.1 20171019)\
Boost C++ libraries (https://www.boost.org/ tested on version 1.67.0)\
CMake (https://cmake.org/install/ tested on version 3.5.1)\
SWIG (http://www.swig.org/ tested on version 3.0.8)\
ZLIB ("sudo apt-get install zlib1g-dev" on Ubuntu 18)\
BZip2 ("sudo apt-get install libbz2-dev" on Ubuntu 18)

## Installation
```
cd /where/you/want/to/install/ReSequenceR
git clone https://github.com/schmeing/ReSequenceR.git
cd ReSequenceR
mkdir build
cd build
cmake ..
make
```
The executable file will afterwards be `/where/you/want/to/install/ReSequenceR/ReSequenceR/build/bin/reseq` and can be either added to the PATH variable or copied to the desired place

## Included libraries
Googletest (https://github.com/google/googletest.git, BSD 3-Clause license)\
NLopt (https://github.com/stevengj/nlopt.git, MIT license)\
SeqAn (https://github.com/seqan/seqan, BSD 3-Clause license)\
skewer (https://github.com/relipmoc/skewer, MIT license, made slight adaptations to the code to be able to include it)

## Publication
In preparation

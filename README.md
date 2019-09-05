# ReSequenceR
More realistic simulator for genomic DNA sequences from Illumina machines that achieves a similar k-mer spectrum as the original sequences

## Requirements

| Requirement               | Ubuntu/Debian                      | CentOS | Manual installation | Comments | Tested version |
|---------------------------|------------------------------------|--------|---------------------|----------|----------------|
| Linux system              |
| Compiler supporting C++14 | `sudo apt install build-essential` | `sudo yum install gcc gcc-c++ glibc-devel make` | | On CentOs 7 the g++ compiler is too old to support C++14 so you need to additionally to the yum command install a newer version following for example this guide (https://linuxhostsupport.com/blog/how-to-install-gcc-on-centos-7/). Afterwards the `CXX` variable has to be set to `/usr/local/bin/c++` and `/usr/local/lib64/` has to be added to your `LD_LIBRARY_PATH` before installing boost if you used the standard install path. | 7.2.1 20171019|
| ZLIB | `sudo apt-get install zlib1g-dev` | `sudo yum install zlib-devel` | | | |
| BZip2 | `sudo apt-get install libbz2-dev` | `sudo yum install bzip2-devel` | | | |
| Python libraries | `sudo apt-get install python-dev` | `ssudo yum install python-devel` | | | |
| Git | `sudo apt-get install git` | `sudo yum install git` | | | |
| CMake | `sudo apt-get install cmake` | (too old version) | https://cmake.org/install/  | | 3.5.1 |
| Boost C++ libraries | `sudo apt-get install libboost-all-dev` | `sudo yum install boost-devel`(version not working) | https://www.boost.org/doc/libs/1_71_0/more/getting_started/unix-variants.html | Only download and extraction in section 1 and library builds in section 5 are strictly needed, if you set a prefix you need to add `prefix/lib/` to your `LD_LIBRARY_PATH` and set `BOOST_ROOT` to `prefix` before the installation process below or you will get boost library errors at the cmake and make step. If you manually installed g++ run `./b2` without sudo so the environment variable CXX is found. | 1.67.0 |
| SWIG | `sudo apt-get install swig` | (too old version) | http://www.swig.org/Doc4.0/Preface.html | If you set a prefix you need to add prefix/bin to your PATH variable | 3.0.8 |

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
The executable file will afterwards be `/where/you/want/to/install/ReSequenceR/ReSequenceR/build/bin/reseq` and can be either added to the PATH variable or copied to the desired place.

To test the installation run:
```
reseq test
```

Some useful python scripts can be found in `/where/you/want/to/install/ReSequenceR/ReSequenceR/python`.

## Quick start examples
To create simulated data similar to real data you first need to map the real data to a reference. For example with `bowtie2`:
```
bowtie2-build my_reference.fa my_reference
bowtie2 -p 32 -X 2000 -x my_reference -1 my_data_1.fq -2 my_data_2.fq | samtools sort -m 10G -@ 4 -T _tmp -o my_mappings.bam -
```

To run the full simulation pipeline (Stats creation, Probability estimation, Simulation) execute:
```
reseq illuminaPE -j 32 -r my_reference.fa -b my_mappings.bam -1 my_simulated_data_1.fq -2 my_simulated_data_2.fq
```

The same is done by the following three commands for the three different steps (So you can run for example only the simulation the second time you want to simulate from the same real data):
```
reseq illuminaPE -j 32 -r my_reference.fa -b my_mappings.bam --statsOnly
reseq illuminaPE -j 32 -s my_mappings.bam.reseq --stopAfterEstimation
reseq illuminaPE -j 32 -R my_reference.fa -s my_mappings.bam.reseq --ipfIterations 0 -1 my_simulated_data_1.fq -2 my_simulated_data_2.fq
```
In order to add variation (to simulate diploid genomes or populations) the parameter `-V` needs to be added:
```
reseq illuminaPE -j 32 -r my_reference.fa -b my_mappings.bam -V my_variation.vcf -1 my_simulated_data_1.fq -2 my_simulated_data_2.fq
```
or
```
reseq illuminaPE -j 32 -R my_reference.fa -s my_mappings.bam.reseq -V my_variation.vcf --ipfIterations 0 -1 my_simulated_data_1.fq -2 my_simulated_data_2.fq
```

To run a simulation with tiles the tile information needs to stay in the read names after the mapping. This means there must not be a space before it, like it is often the case for read archive data. To replace the space on the fly with an underscore the reseq-prepare-names.py script is provided. In this case run the mapping like this:
```
bowtie2 -p 32 -X 2000 -x my_reference -1 <(reseq-prepare-names.py my_data_1.fq my_data_2.fq) -2 <(reseq-prepare-names.py my_data_2.fq my_data_1.fq) | samtools sort -m 10G -@ 4 -T _tmp -o my_mappings.bam -
```

## Included libraries
Googletest (https://github.com/google/googletest.git, BSD 3-Clause license)\
NLopt (https://github.com/stevengj/nlopt.git, MIT license)\
SeqAn (https://github.com/seqan/seqan, BSD 3-Clause license)\
skewer (https://github.com/relipmoc/skewer, MIT license, made slight adaptations to the code to be able to include it)

## Publication
In preparation

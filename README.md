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
## File Formats
| File type             | Ending     | Information |
|-----------------------|------------|-------------|
| Reference file        | .fa        | Standard reference file in fasta format. Everything not an A,C,G,T is randomly replaced by one of those before simulating. To consistently simulate the same base over multiple simulations use the replaceN mode of reseq before starting the simulation to create a reference without N(or other ambiguous bases, which are all treated as N). |
| Mapping file          | .bam       | Standard mapping file in bam format. Reference information need to match the provided reference. Soft clipped bases are yet untested. To include tiles the tile information needs to stay in the QNAME field. |
| Adapter file          | .fa        | File in fasta format giving the sequences of possibly used adapters. The shorter this list the less missidentification can happen. Some files to use are in the adapter folder. The direction of the adapters is irrelevant as always the adapter given and its reverse complement are checked, due to sequencing-machine dependence on the direction. |
| Adapter matrix        | .mat       | 0/1 matrix stating if the adapters can occur in a read pair together. (0:no; 1:yes). The n-th row/column is the n-th entry in the adapter file. Rows represent the adapter in the first read and columns represent the adapter in the second read. Columns are consecutive digits in a row. Number of rows and columns need to match the number of entries in the adapter file. |
| Variant file          | .vcf       | Standard vcf format. The reference information in the header and in the reference column must match the given reference file. Except of this only the CHROM, POS, ALT columns and the genotype information are used. No filtering by quality etc. takes place. All genotypes in the file will be simulated. No distinction is made if genotypes are in a single sample or spread out over multiple samples. All genotype information is considered phased independent of what is encoded in the file. |
| Stats file            | .reseq     | Boost archive: Not recommended to modify or create by hand even though it is in ASCII format |
| Probability file      | .reseq.ipf | Boost archive: Not recommended to modify or create by hand even though it is in ASCII format |
| Systematic error file | .fq        | Standard fastq format. The sequence represents the dominant error and the quality the error rate in percent at that position. There are two entries per reference sequence. The order of the reference sequences must be kept. The length must match the length of the reference sequence. The first entry per reference sequence is the reverse strand and reverse complemented. So an A in the first position means that a systematic error towards a T is simulated for the last base of the reference sequence in reads on the reverse strand. The second entry is the forward strand taken as is, so not reverse complemented. The error rate in percent is encoded similar to quality values with an offset of 33. Since the fastq format is limited to 94 quality values odd percentages over 86 are omitted. This mean `~` encodes 100% and `}` 98%. |
| Reference bias file   | .txt       | One line per specified bias. The line starts by a unique identifier of the reference sequence (the part before the first space in the reference sequence name in the reference file). The identifier can be followed by a space and after it some arbitrary information. The line ends with a space or tab separating a floating point number representing the bias for this sequence. It must be positive. All reference sequences in the reference file must have a bias given. However, the order of the sequences doesn't need to be kept and the bias for additional reference sequences could be specified. |

## Included libraries
Googletest (https://github.com/google/googletest.git, BSD 3-Clause license)\
NLopt (https://github.com/stevengj/nlopt.git, MIT license)\
SeqAn (https://github.com/seqan/seqan, BSD 3-Clause license)\
skewer (https://github.com/relipmoc/skewer, MIT license, made slight adaptations to the code to be able to include it)

## Publication
In preparation

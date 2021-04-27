# ReSeq
More realistic simulator for genomic DNA sequences from Illumina machines that achieves a similar k-mer spectrum as the original sequences.

## Table of Contents

- [Abstract](#abstract)
- [Requirements](#requirements)
- [Installation](#installation)
- [Bioconda](#conda)
- [Quick start examples](#quickstart)
- [Apply errors and qualities directly to sequences](#errormodel)
- [Parameter](#parameter)
- [File Formats](#formats)
- [FAQ](#faq)
- [Included libraries](#libraries)
- [Publication](#publication)

## <a name="abstract"></a>Abstract
Even though sequencing biases and errors have been deeply researched to adequately account for them, comparison studies, e.g. for error correction, assembly or variant calling, face the problem that synthetic datasets resemble the real output of high-throughput sequencers only in very limited ways, resulting in optimistic estimated performance of programs run on simulated data compared to real data. Therefore, comparison studies are often based on real data. However, this approach has its own difficulties, since the ground truth is unknown and can only be estimated, which introduces its own biases and circularity towards easy solutions and the methods used.

**Re**al **Seq**uence Reproducer shortens the gap between simulated and real data evaluations by adequately reproducing key statistics of real data, like the coverage profile, systematic errors and the k-mer spectrum. When these characteristics are translated into new synthetic computational experiments (i.e. simulated data), the performance can be more accurately estimated. Combining our simulator and real data gives two valuable perspectives on the performance of tools to minimize biases.

## <a name="requirements"></a>Requirements

| Requirement               | Ubuntu/Debian                      | CentOS | Manual installation | Comments | Tested version |
|---------------------------|------------------------------------|--------|---------------------|----------|----------------|
| Linux system              |
| Compiler supporting C++14 | `sudo apt install build-essential` | `sudo yum install gcc gcc-c++ glibc-devel make` | | On CentOs 7 the g++ compiler is too old to support C++14 so you need to additionally to the yum command install a newer version following for example this guide (https://linuxhostsupport.com/blog/how-to-install-gcc-on-centos-7/). If the standard install path `usr/local/` was used, afterwards the `CXX` variable has to be set to `/usr/local/bin/c++`, the `CC` variable to `/usr/local/bin/gcc`, `/usr/local/lib64/` has to be added to your `LD_LIBRARY_PATH` and `/usr/local/bin` to your `PATH` before installing boost. | 7.2.1 20171019|
| ZLIB                      | `sudo apt-get install zlib1g-dev`  | `sudo yum install zlib-devel`   | | | |
| BZip2                     | `sudo apt-get install libbz2-dev`  | `sudo yum install bzip2-devel`  | | | |
| Python 3                  |                                    |                                 | | | 3.6.11 |
| Python libraries          | `sudo apt-get install python3-dev`  | `sudo yum install python3-devel.x86_64` | | | |
| Git                       | `sudo apt-get install git`         | `sudo yum install git`          | | | |
| CMake                     | `sudo apt-get install cmake`       | (too old version)               | https://cmake.org/install/  | The newest versions (starting 3.16) require `sudo apt-get install libssl-dev` or `sudo yum install openssl-devel` | 3.5.1 |
| Boost C++ libraries       | `sudo apt-get install libboost-all-dev` | (version not working) | https://www.boost.org/doc/libs/1_71_0/more/getting_started/unix-variants.html | Only download and extraction in section 1 and library builds in section 5 are strictly needed, if you set a prefix you need to set `BOOST_ROOT` to this `prefix` before the installation process below or you will get boost library errors at the cmake and make step. If you manually installed g++ run `./b2` without sudo so the environment variables `CXX` and `CC` are found. | 1.67.0 |
| SWIG                      | `sudo apt-get install swig`        | (too old version)               | http://www.swig.org/Doc4.0/Preface.html | If you set a prefix you need to add prefix/bin to your PATH variable | 3.0.8 |

## <a name="installation"></a>Installation
To install to the standard folder `usr/local` or to keep everything in the build folder:
```
cd /where/you/want/to/build/ReSeq
git clone https://github.com/schmeing/ReSeq.git
cd ReSeq
mkdir build
cd build
cmake ..
make
```

To install to a different folder the same steps apply but the `cmake ..` line has to be exchange with:
```
cmake -DCMAKE_INSTALL_PREFIX=/where/you/want/to/install/ReSeq/ ..
```

The executable file will afterwards be `/where/you/want/to/build/ReSeq/ReSeq/build/bin/reseq` and can be added to the PATH variable or copied to the desired place.

Alternatively ReSeq can be install to the standard folder `usr/local` or the previously defined folder by:
```
make install
```

To test the installation run:
```
reseq test
```

Some useful python scripts can be found in `/where/you/want/to/install/ReSeq/ReSeq/python` or after an installation in `usr/local/bin` or `/where/you/want/to/install/ReSeq/bin/`.

## <a name="conda"></a>Bioconda
ReSeq can also be installed in an automatic fashion via anaconda/miniconda(https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html) with the following command:
```
conda install -c bioconda -c conda-forge reseq
```
However, updates will not be as frequent and the option to switch to the devel branch to get the most recent bugfixes is missing.

## <a name="quickstart"></a>Quick start examples
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
For best results it is always advised to create your own profiles from a dataset very closely matching the desired sequencer, chemistry, fragmentation, adapters, PCR cycles, etc. Furthermore, training on the same or a closely related species, is best to be sure that the necessary profile space (e.g. range of GC content) is well populated. However, in many situations there is no specific case that should be simulated, but a wide variety of datasets is important. Under this condition finding good datasets is tedious work and recreating the same profile from a given dataset does not help to ensure the quality of the simulated data. Therefore, [this repository](https://github.com/schmeing/ReSeq-profiles) is designed to provide high-quality profiles with detailed information on the original datasets. Whenever possible, method benchmarks should be performed on the simulated and the original dataset to verify that the simulation created realistic conditions for the particular use case.

## <a name="errormodel"></a>Apply errors and qualities directly to sequences
In case you cannot use the coverage model, you can directly provide sequences to ReSeq, which will be converted to reads. This includes adding qualities as well as InDel and substitution errors and cutting the sequence to the read length. If sequences are shorter than the read length ReSeq automatically adds an adapter.
```
reseq seqToIllumina -j2 -i my_sequences.fa -o my_simulated_reads.fq -s my_stats_profile.reseq
```
For it to work, all necessary informations need to be provided to ReSeq's error and quality model. Therefore, each input sequence in the fasta file must have the following form:
```
>{sequence id} {template segment};{fragment length};{error tendencies};{error rates}
{sequence to convert}
```
`{sequence id}`: The desired sequence id. It can contain spaces. The final read description in the output fastq will be:
@{sequence id} {cigar} E{number of errors in read}`

`{template segment}`: 1 for first reads or 2 for second reads.

`{sequence to convert}`: Sequence to which errors and qualities will be added. It may only contain A, C, G or T. Ns are not permitted, since a conversion should be performed in a consistent manner for all reads stemming from a given position in the reference (see `reseq replaceN`).

`{error tendencies}` and `{error rates}`: Both must have the same length as `{sequence to convert}` and the nth position always corresponds to the nth position in the sequence. Their format is described in detail for the Systematic error file under [File Formats](#formats). All bases stemming from the same position and strand in the reference must have the same values to properly simulate systematic errors. For insertions that are specific to one read use `N` and `!` respectively. The systematic errors for a reference can be simulated by calling:

`reseq illuminaPE -r my_reference.fa -s my_stats_profile.reseq --stopAfterEstimation --writeSysError my_systematic_errors.fq`

This call creates a fastq file with two sequences per reference sequence (one for each strand with the reverse strand first). From this file the corresponding error tendencies and rates can be extracted. Note that the the sequence for the reverse strand is already reverse complemented.

## <a name="parameter"></a>Parameter
`reseq illuminaPE [options]`

| Parameter         | Default | Description |
|-------------------|---------|-------------|
| **General**       |
| `-h` `--help`     |         | Prints help information and exits |
| `-j` `--threads`  | 0       | Number of threads used (0=auto) |
| `--verbosity`     | 4       | Sets the level of verbosity (4=everything, 0=nothing) |
| `--version`       |         | Prints version info and exits |
| **Stats**         |
| `--adapterFile`   | (AutoDetect) | Fasta file with adapter sequences |
| `--adapterMatrix` | (AutoDetect) | 0/1 matrix with valid adapter pairing (first read in rows, second read in columns) |
| `-b` `--bamIn`    | None    | Position sorted bam/sam file with reads mapped to `--refIn` |
| `--binSizeBiasFit`| 100000000 | Reference sequences large then this are split for bias fitting to limit memory consumption |
| `--maxFragLen`    | 2000    | Maximum fragment length to include pairs into statistics |
| `--minMapQ`       | 10      | Minimum mapping quality to include pairs into statistics |
| `--noBias`        |         | Do not perform bias fit. Results in uniform coverage if simulated from |
| `--noTiles`       |         | Ignore tiles for the statistics [default] |
| `-r` `--refIn`    | None    | Reference sequences in fasta format (gz and bz2 supported) |
| `--statsOnly`     |         | Only generate the statistics |
| `-s` `--statsIn`  | None    | Skips statistics generation and reads directly from stats file |
| `-S` `--statsOut` | `--bamIn`.reseq | Stores the real data statistics for reuse in given file |
| `--tiles`         |         | Use tiles for the statistics |
| `-v` `--vcfIn`    | None    | Ignore all positions with a listed variant for stats generation |
| **Probabilities** |
| `--ipfIterations` | 200     | Maximum number of iterations for iterative proportional fitting |
| `--ipfPrecision`  | 5       | Iterative proportional fitting procedure stops after reaching this precision (%) |
| `-p` `--probabilitiesIn` | `--statsIn`.ipf | Loads last estimated probabilities and continues from there if precision is not met |
| `-P` `--probabilitiesOut` | `--probabilitiesIn` | Stores the probabilities estimated by iterative proportional fitting |
| `--stopAfterEstimation` |   | Stop after estimating the probabilities |
| **Simulation**    |
| `-1` `--firstReadsOut` | reseq-R1.fq | Writes the simulated first reads into this file |
| `-2` `--secondReadsOut` | reseq-R2.fq | Writes the simulated second reads into this file |
| `-c` `--coverage` | 0       | Approximate average read depth simulated (0 = Corrected original coverage) |
| `--errorMutliplier` | 1.0   | Divides the original probability of correct base calls(no substitution error) by this value and renormalizes |
| `--methylation`   |         | Extended bed graph file specifying methylation for regions. Multiple score columns for individual alleles are possible, but must match with vcfSim. C->T conversions for 1-specified value in region. |
| `--noInDelErrors` |         | Simulate reads without InDel errors |
| `--noSubstitutionErrors` |  | Simulate reads without substitution errors |
| `--numReads`      | 0       | Approximate number of read pairs simulated (0 = Use `--coverage`) |
| `--readSysError`  | None    | Read systematic errors from file in fastq format (seq=dominant error, qual=error percentage) |
| `--recordBaseIdentifier` | ReseqRead | Base Identifier for the simulated fastq records, followed by a count and other information about the read |
| `--refBias`       | keep/no | Way to select the reference biases for simulation (keep [from refIn] / no [biases] / draw [with replacement from original biases] / file) |
| `--refBiasFile`   | None    | File to read reference biases from: One sequence per file (identifier bias) |
| `-R` `--refSim`   | `--refIn` | Reference sequences in fasta format to simulate from |
| `--seed`          | None    | Seed used for simulation, if none is given random seed will be used |
| `-V` `--vcfSim`   | None    | Defines genotypes to simulate alleles or populations |
| `--writeSysError` | None    | Write the randomly drawn systematic errors to file in fastq format (seq=dominant error, qual=error percentage) |

 `reseq queryProfile [options]`
 
| Parameter         | Default | Description |
|-------------------|---------|-------------|
| **General**       |
| `-h` `--help`     |         | Prints help information and exits |
| `-j` `--threads`  | 0       | Number of threads used (0=auto) |
| `--verbosity`     | 4       | Sets the level of verbosity (4=everything, 0=nothing) |
| `--version`       |         | Prints version info and exits |
| **queryProfile**  |
| `--maxLenDeletion`| None    | Output lengths of longest detected deletion to stdout |
| `--maxReadLength` | None    | Output lengths of longest detected read to stdout |
| `-r` `--ref`      | None    | Reference sequences in fasta format (gz and bz2 supported) |
| `--refSeqBias`    | None    | Output reference sequence bias to file (tsv format; `-` for stdout) |
| `-s` `--stats`    | None    | Reseq statistics file to extract reference sequence bias |

`reseq replaceN [options]`

| Parameter         | Default | Description |
|-------------------|---------|-------------|
| **General**       |
| `-h` `--help`     |         | Prints help information and exits |
| `-j` `--threads`  | 0       | Number of threads used (0=auto) |
| `--verbosity`     | 4       | Sets the level of verbosity (4=everything, 0=nothing) |
| `--version`       |         | Prints version info and exits |
| **ReplaceN**      |
| `-r` `--refIn`    | None    | Reference sequences in fasta format (gz and bz2 supported) |
| `-R` `--refSim`   | None    | File to where reference sequences in fasta format with N's randomly replace should be written to |
| `--seed`          | None    | Seed used for replacing N, if none is given random seed will be used |

`reseq seqToIllumina [options]`

| Parameter         | Default | Description |
|-------------------|---------|-------------|
| **General**       |
| `-h` `--help`     |         | Prints help information and exits |
| `-j` `--threads`  | 0       | Number of threads used (0=auto) |
| `--verbosity`     | 4       | Sets the level of verbosity (4=everything, 0=nothing) |
| `--version`       |         | Prints version info and exits |
| **seqToIllumina** |
| `--errorMutliplier` | 1.0   | Divides the original probability of correct base calls(no substitution error) by this value and renormalizes |
| `-i` `--input`    | `stdin` | Input file (fasta format, gz and bz2 supported) |
| `--ipfIterations` | 200     | Maximum number of iterations for iterative proportional fitting |
| `--ipfPrecision`  | 5       | Iterative proportional fitting procedure stops after reaching this precision (%) |
| `--noInDelErrors` |         | Simulate reads without InDel errors |
| `--noSubstitutionErrors` |  | Simulate reads without substitution errors |
| `-o` `--output`   | `stdout`| Output file (fastq format, gz and bz2 supported) |
| `-p` `--probabilitiesIn` | `--statsIn`.ipf | Loads last estimated probabilities and continues from there if precision is not met |
| `-P` `--probabilitiesOut` | `--probabilitiesIn` | Stores the probabilities estimated by iterative proportional fitting |
| `--seed`          | None    | Seed used for simulation, if none is given random seed will be used |
| `-s` `--statsIn`  | None    | Profile file that contains the statistics used for simulation |

## <a name="formats"></a>File Formats
| File type             | Ending     | Information |
|-----------------------|------------|-------------|
| Reference file        | .fa        | Standard reference file in fasta format. Everything not an A,C,G,T is randomly replaced by one of those before simulating. To consistently simulate the same base over multiple simulations use the replaceN mode of reseq before starting the simulation to create a reference without N(or other ambiguous bases, which are all treated as N). A single reference sequence is only supported up to  a maximum length of 4294967295 bases. The complete reference can be much bigger. |
| Mapping file          | .bam       | Standard mapping file in bam format. Reference information need to match the provided reference. To include tiles the tile information needs to stay in the QNAME field. Only primary alignments are used. |
| Adapter file          | .fa        | File in fasta format giving the sequences of possibly used adapters. The shorter this list the less missidentification can happen. Some files to use are in the adapter folder. The direction of the adapters is irrelevant as always the adapter given and its reverse complement are checked, due to sequencing-machine dependence on the direction. |
| Adapter matrix        | .mat       | 0/1 matrix stating if the adapters can occur in a read pair together. (0:no; 1:yes). The n-th row/column is the n-th entry in the adapter file. Rows represent the adapter in the first read and columns represent the adapter in the second read. Columns are consecutive digits in a row. Number of rows and columns need to match the number of entries in the adapter file. |
| Variant file          | .vcf       | Standard vcf format. The reference information in the header and in the reference column must match the given reference file. Ambiguous bases (e.g. N's) are not supported in the reference and alternative column. Except of this only the CHROM, POS, ALT columns and the genotype information are used. No filtering by quality etc. takes place. All genotypes in the file will be simulated. No distinction is made if genotypes are in a single sample or spread out over multiple samples. All genotype information is considered phased independent of what is encoded in the file. At most 128 genotypes are supported by default (see FAQ). |
| Methylation file	| .bed	     | Extended bed graph format with possibility of multiple score columns for individual alleles. Number of columns must be either 1 or the same number of alleles specified in the variant file. The number of columns need to be the same within each reference sequence. Bisulfite sequencing is simulated, so C->T conversions are inserted with a probability of 1-methylation value specified in this file. |
| Stats file            | .reseq     | Boost archive: Not recommended to modify or create by hand even though it is in ASCII format |
| Probability file      | .reseq.ipf | Boost archive: Not recommended to modify or create by hand even though it is in ASCII format |
| Systematic error file | .fq        | Standard fastq format. The sequence represents the error tendency and the quality the error rate in percent at that position. There are two entries per reference sequence. The order of the reference sequences must be kept. The length must match the length of the reference sequence. The first entry per reference sequence is the reverse strand and reverse complemented. So an A in the first position means that a systematic error towards a T is simulated for the last base of the reference sequence in reads on the reverse strand. The second entry is the forward strand taken as is, so not reverse complemented. The error rate in percent is encoded similar to quality values with an offset of 33. Since the fastq format is limited to 94 quality values odd percentages over 86 are omitted. This mean `~` encodes 100% and `}` 98%. |
| Reference bias file   | .txt       | One line per specified bias. The line starts by a unique identifier of the reference sequence (the part before the first space in the reference sequence name in the reference file). The identifier can be followed by a space and after it some arbitrary information. The line ends with a space or tab separating a floating point number representing the bias for this sequence. It must be positive. All reference sequences in the reference file must have a bias given. However, the order of the sequences doesn't need to be kept and the bias for additional reference sequences could be specified. The biases will be automatically normalized and define the relative base coverage of the reference sequences. Simulated base coverage will differ from this due to other biases additionally taken into account. |

## <a name="faq"></a>FAQ
**Can I simulate more than the default 128 alleles?**
Yes. You can set `kMaxAlleles` in `reseq/Reference.h` to any multiple of 64. After recompilation the new maximum number of alleles will be your chosen value.

**Can I simulate exome sequencing?**
Yes. You need to use a reference that only contains the exons as individual scaffolds. Using `--refBiasFile` you can specify the coverage of individual exons. To simulate intron contamination you can add the whole reference to the reference containing the exons and strongly reduce the coverage for these scaffolds using `--refBiasFile`.

**Can I train on datasets without adapters?**
Generally, it is not advised to use trimmed datasets, because they result in worse perfomance. However, by specifying decoy adapters with `--adapterFile TruSeq_single` you can skip the automatic adapter detection, which otherwise will prevent you from training on datasets without adapters.

**When I train the model a large part of the genome is excluded, because the sequences are too short.**
Lowering the `--maxFragLen` parameter most likely helps in this situation, because sequences that are not at least 100 bases longer than this parameter are excluded in any case. However, you need to check that you are not truncating your fragment lengths distribution by setting `--maxFragLen` too low.

## <a name="libraries"></a>Included libraries
Googletest (https://github.com/google/googletest.git, BSD 3-Clause license)\
NLopt (https://github.com/stevengj/nlopt.git, MIT license)\
SeqAn (https://github.com/seqan/seqan, BSD 3-Clause license)\
skewer (https://github.com/relipmoc/skewer, MIT license, made slight adaptations to the code to be able to include it)

## <a name="publication"></a>Publication
Schmeing, S., Robinson, M.D. ReSeq simulates realistic Illumina high-throughput sequencing data. Genome Biol 22, 67 (2021). https://doi.org/10.1186/s13059-021-02265-7

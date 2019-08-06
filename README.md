# tamock
Targeted bacterial mock communities

tamock creates sample-specific mock communities of metagenomic samples
by classifying all sequence reads and replacing all bacterial reads
with corresponding in silico sequences sampled randomly from RefSeq
genomes.
Classification is implemented using Centrifuge (Kim *et al.* 2016) and
sequencing read simulation is done by ART (Huang *et al.* 2012).


## Installation
### Dependencies

* GNU Scientific Library (GSL) https://www.gnu.org/software/gsl/
* Perl >= v5.12.0

Debian/Ubuntu

```bash
 apt install libgsl-dev
 #older distributions
 apt install libgsl0-dev
```

CentOS
``` bash
 yum install gsl-devel
```

MacOS
```bash
 brew install gsl
```

Centrifuge indexes have to be downloaded separately.
Prebuild indexes are provided on the centrifuge homepage at http://www.ccb.jhu.edu/software/centrifuge/

Example:

```bash
 #for details on provided indexes, please see documentation at the centrifuge homepage
 mkdir /<tamock_install_dir>/centrifuge-index
 cd /<tamock_install_dir>/centrifuge-index
 wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p+h+v.tar.gz
 tar -xzf p+h+v.tar.gz

 #if there is already a present centrifuge installation with indexes, the index folder can also be linked like
 ln -s /path/to/centrifuge-indexes /<tamock_install_dir>/centrifuge-index
```

If a centrifuge p+h+v index is placed at the install directory under
<install_dir/centrifuge-index/p+h+v*>
this index is used by default if no other index is provided to tamock.

RefSeq reference genomes are all saved within one directory. If the directory /<tamock_install_dir>/refseq-genomes is created, this directory is used by default unless another directory has to be provided via the command line.

```bash
 #optional
 mkdir /<tamock_install_dir>/refseq-genomes
```

To speed up analysis, all bacterial RefSeq genomes can be downloaded before mock community creation (>100 GB) instead of only downloading required genomes during a single run.

```bash
 #optional
 tamock --download-refseq -R /path/to/refseq-genomes
```

### Install

```bash
 git clone https://github.com/gerners/tamock.git
 cd tamock
 tamock --install-deps
 tamock --install-test
```

## Quick start

Paired end simulation with 4 threads

```bash
 tamock -1 <paired_1.fastq> -2 <paired_2.fastq> -o <output directory> \
 -R <reference genome directory> -x <centrifuge-index/p+h+v> -t 4
```

Paired end simulation with additional unpaired reads

```bash
 tamock -1 <paired_1.fastq> -2 <paired_2.fastq> -U <single.fastq>  \
 -o <output directory> -R <reference genome directory> -x <centrifuge-index/p+h+v>
```

Single end simulation

```bash
 tamock -U <single.fastq> -o <output directory> -R <reference genome directory> \
 -x <centrifuge-index/p+h+v>
```

Classification is the most expensive task. Multiple threads speed up classification while RAM usage depends on the index used. 
This step can also be skipped by providing pre-computed centrifuge classification results.

CAVEAT: When using extern centrifuge results, ensure that centrifuge version 1.0.4 or higher is used and that the same index is provided to tamock.

```bash
 tamock -1 <paired_1.fastq> -2 <paired_2.fastq> --centrifuge-kreport -o <output directory> \
 -R <reference genome directory> -x <centrifuge-index/p+h+v>
```

### Options

#### Mandatory

* -1/--forward

Forward paired reads. Can be used together with unpaired reads using -U option.

* -2/--reverse

Reverse paired reads. Can be used together with unpaired reads using -U option.

* -U/--unpaired

Unpaired or single end reads. Required if no paired reads are present.

* -o/--outdir
	
Output directory

* -R/--refgenomes
	
Directory to store all RefSeq reference genomes.

#### Optional

* --gzip

Gzip all sequence files

* -x/--index

Centrifuge index (*.cf files). Only basename needed, e.g. /path/to/index.X.cf should be provided as '-x /path/to/index'.  
Default: p+h+v if files present at /installdir/centrifuge-index/p+h+v.X.cf

Index files can be downloaded from the centrifuge homepage at http://www.ccb.jhu.edu/software/centrifuge/

* -a/--assembly-summary

NCBI assembly summary table from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt.  
Defaults to /installdir/assembly_summary.txt (downloaded during installation).

* --ncbi-sum-update

Update NCBI assembly summary to current version and replace old version at /installdir/assembly_summary.txt

* --download-refseq

Download all reference genomes from provided assembly-summary table. Requires option -R (directory to safe reference genomes)

* -t/--threads

Number of threads used by centrifuge. Defaults to 1.

* -v/--verbose

Verbose mode, print detailed information to screen

==== Mocks ====

* --centrifuge-out

Centrifuge tabbed output file (centrifuge -S <file>). If centrifuge-kreport and centrifuge output files are provided, classification is skipped and these files are used instead. 
CAVEAT: The same index used for classification has to be provided to tamock as well via -x/--index.

Centrifuge example command:
centrifuge -x index -1 forwardpairs.fq -2 reversepairs.fq --out-fmt tab -S centrifuge.out --report-file centrifuge.report

* --centrifuge-kreport

Centrifuge kraken-like report. If centrifuge-kreport and centrifuge output files are provided, classification is skipped and these files are used instead. 
CAVEAT: The same index used for classification has to be provided to tamock as well via -x/--index.
CAVEAT: Centrifuge version 1.0.4. or higher required

Centrifuge example command:
centrifuge-kreport -x index centrifuge.out > centrifuge.kreport

* --min-abund-species

Only include species with at least x classified reads into the mock community. Values above 1 effectively filter out spurious
classifications, eliminating entries with less than x reads classified to species and relating strains. Defaults to 1.

* --no-reassign

No reassignment of read counts from strains without a reference genome to other classified strains/reference genomes of same species.
Classifications to strains with no reference genome might occur from non-matching strain IDs from the centrifuge index to the current
NCBI RefSeq collection. Default off.

==== ART sequence simulator ====

* -M/--illumina-model

Either 'custom' to calculate error profile of input reads or use precalculated ART Illumina error model. Available precalculated profiles are:
GA1 - GenomeAnalyzer I (36bp,44bp),	GA2 - GenomeAnalyzer II (50bp, 75bp),
HS10 - HiSeq 1000 (100bp),			HS20 - HiSeq 2000 (100bp),			HS25 - HiSeq 2500 (125bp, 150bp),
HSXn - HiSeqX PCR free (150bp),		HSXt - HiSeqX TruSeq (150bp),		MinS - MiniSeq TruSeq (50bp),
MSv1 - MiSeq v1 (250bp),			MSv3 - MiSeq v3 (250bp),			NS50 - NextSeq500 v2 (75bp)

Defaults to sample specific error profile determination (custom).

* --qprof1

Precalculated forward-read quality profile for custom error profile in ART, together with -M custom. Replaces calculation of error profiles of input sequences.

* --qprof2

Precalculated reverse-read quality profile for custom error profile in ART, together with -M custom. Replaces calculation of error profiles of input sequences.

* -l/--length

Read length of simulated reads. Should be matched with a realistic error profile with according read lengths (-M). If not provided, determined from first 1000 reads of input sequences.

* --mean-fragment-length

Mean size of frament length for paired-end simulations by ART. Defaults to 200.

* --sd-fragment-length

Standard deviation of fragment size for paired-end simulations. Defaults to 10.

==== Dependencies ==== 

* --install-deps

Install centrifuge and ART default versions into the installation directory of tamock. Depends on libgsl (See Installation).

* --install-test

Test if installation was successfull.
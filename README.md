# tamock
Benchmark data creation

Tamock simulates habitat-specific benchmark data for metagenomic samples.

Simulated metagenomic reads are widely used to benchmark software and workflows for metagenome interpretation. 
Scope and power of metagenomic benchmarks depend on the selection of their underlying communities. As a result, 
conclusions of benchmark studies are limited for distant communities towards the benchmark data used. Ideally,  
simulations are therefore based on genomes, which resemble metagenomic communities realistically. 

Tamock facilitates the simulation of metagenomic reads according to a microbial community, derived from real 
metagenomic data. Thus, Tamock simulations enable an assessment of computational methods, workflows and parameters 
specific for a microbial habitat. Tamock automatically determines taxonomic profiles from shotgun metagenomic data, 
selects reference genomes accordingly and uses them to simulate metagenomic reads. 

By default, the bacterial fraction of a community is simulated, however other domains such as Eukaryota, Archaea or Viruses
can be simulated as well (-d/--domains option).

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

### Install

```bash
 git clone https://github.com/gerners/tamock.git
 cd tamock
 
 #if ART (version 20160605) and Centrifuge (>= v1.0.4) are available on the system, this step is not needed
 tamock --install-deps
 
 tamock --install-test
```

Centrifuge indexes have to be downloaded separately.
Prebuild indexes are provided on the centrifuge homepage at http://www.ccb.jhu.edu/software/centrifuge/

Example:

```bash
 #for details on provided indexes, please see documentation at the centrifuge homepage (link above)
 mkdir /<tamock_install_dir>/centrifuge-index
 cd /<tamock_install_dir>/centrifuge-index
 wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p+h+v.tar.gz
 tar -xzf p+h+v.tar.gz

 #if there is already a present centrifuge installation with indexes, the index folder can also be linked as follows to be used by default
 ln -s /path/to/centrifuge-index /<tamock_install_dir>/centrifuge-index
```

If a centrifuge p+h+v index is placed at the install directory under
<install_dir/centrifuge-index/p+h+v*>
this index is used by default if no other index is provided to tamock.

RefSeq reference genomes are all saved within one directory. If the directory /<tamock_install_dir>/refseq-genomes is created, this directory is used by default unless another directory is provided via the command line.

**CAVEAT**: If Eukaryotic genomes are simulated, downloaded genomes could potentially use >100 GB of file space 

```bash
 #optional
 mkdir /<tamock_install_dir>/refseq-genomes
```

If a local copy of all RefSeq genomes is present on the system, softlinks to all genome files (e.g. GCF_000013425.1_ASM1342v1_genomic.fna.gz) can be placed in /<tamock_install_dir>/refseq-genomes to use the local mirror of RefSeq.

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

CAVEAT: When using external centrifuge results, ensure that centrifuge version 1.0.4 or higher is used and that the same index used is provided to tamock.

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
	
Directory to store all RefSeq reference genomes. Defaults to <install_dir/refseq-genomes> if present

* -x/--index

Centrifuge index (*.cf files). Only basename needed, e.g. /path/to/index.X.cf should be provided as '-x /path/to/index'.  
Defaults to <installdir/centrifuge-index/p+h+v> if present and no other option given

Index files can be downloaded from the centrifuge homepage at http://www.ccb.jhu.edu/software/centrifuge/

#### Optional

* -d/--domains

Select domains which should be simulated. Options are E(ukaryota), B(acteria), V(iruses), A(rchaea)
For multiple selections, provide comma separated values e.g. <-d E,B>. Defaults to Bacteria (B).

* --gzip

Gzip all sequence files

* -k/--keep
Keep all intermediate ART and sequence files (sequence fraction with (*repl_by_sim*)/without (*not_repl*) reference genomes, 
simulated sequences (simulated*.fq)). Temporary files are stored in <output_dir/tmpreads>	

* -a/--assembly-summary

NCBI assembly summary table from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt.  
Defaults to /installdir/assembly_summary.txt (downloaded during installation).

* --ncbi-sum-update

Update NCBI assembly summary to current version and replace old version at /installdir/assembly_summary_refseq.txt
Alternatively manually download from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt and provide via -a

* -t/--threads

Number of threads used by centrifuge. Defaults to 1.

* -v/--verbose

Verbose mode, print detailed information to screen

==== Benchmarks ====

* --centrifuge-kreport

Centrifuge kraken-like report. If centrifuge-kreport and centrifuge output files are provided, classification is skipped and these files are used instead. 
CAVEAT: The same index used for classification has to be provided to tamock via -x/--index.
CAVEAT: Centrifuge version 1.0.4. or higher required

Centrifuge example command:
centrifuge-kreport -x index centrifuge.out > centrifuge.kreport

* --centrifuge-out

Centrifuge tabbed output file (centrifuge -S <file>). If centrifuge-kreport and centrifuge output files are provided, classification is skipped and these files are used instead. 
CAVEAT: The same index used for classification has to be provided to tamock as well via -x/--index.

Centrifuge example command:
centrifuge -x index -1 forwardpairs.fq -2 reversepairs.fq --out-fmt tab -S centrifuge.out --report-file centrifuge.report

* --no-reassign

No reassignment of read counts from strains without a reference genome to other classified strains/reference genomes of same species.
Only simulate sequences classified directly to a reference. Default off.

==== ART sequence simulator ====

* -M/--illumina-model

Either 'custom' to calculate error profile of input reads or use precalculated ART Illumina error model. Available precalculated profiles are:
GA1 - GenomeAnalyzer I (36bp,44bp),	GA2 - GenomeAnalyzer II (50bp, 75bp),
HS10 - HiSeq 1000 (100bp),			HS20 - HiSeq 2000 (100bp),			HS25 - HiSeq 2500 (125bp, 150bp),
HSXn - HiSeqX PCR free (150bp),		HSXt - HiSeqX TruSeq (150bp),		MinS - MiniSeq TruSeq (50bp),
MSv1 - MiSeq v1 (250bp),			MSv3 - MiSeq v3 (250bp),			NS50 - NextSeq500 v2 (75bp)

Defaults to sample specific error profile, calculated on-the-fly (custom).

* --qprof1

Pre-calculated forward-read quality profile for custom error profile in ART. Replaces calculation of error profiles of input sequences.

* --qprof2

Pre-calculated reverse-read quality profile for custom error profile in ART. Replaces calculation of error profiles of input sequences.

* -l/--length

Read length of simulated reads. Should be matched with a realistic error profile with according read lengths (-M). 
If not provided, the longest read length from the first 250 reads of input sequences are used.

* --rn-sim

Number of reads for simulated sequence fraction. By default, the same number of sequences from the input sample will be kept, this option alters number of sequences!

* --mean-fragment-length

Mean size of fragment length for paired-end simulations by ART. Defaults to 200.

* --sd-fragment-length

Standard deviation of fragment size for paired-end simulations. Defaults to 10.

==== Dependencies ==== 

* --install-deps

Install centrifuge and ART default versions into the installation directory of tamock. Depends on libgsl (See Installation).

* --install-test

Test installation.

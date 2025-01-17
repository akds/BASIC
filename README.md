## BASIC: BCR and TCR assembly from single cell RNA-seq

[![Build Status](https://travis-ci.org/akds/BASIC.svg?branch=master)](https://travis-ci.org/akds/BASIC)

### Prerequisites
* BASIC is tested to work on Python 2.7 and 3.5+
* BASIC requires Bowtie2 to run.
* The `db/` folder must remain with the BASIC.py file

### Installation
BASIC can be installed with [bioconda](https://bioconda.github.io/#using-bioconda). Once conda and bioconda are set up, running the following will install BASIC and [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml):

```bash
conda install basic
```

Alternatively, if you would like to use the development version and have Bowtie2 installed, this repository can be cloned and BASIC can be run as shown below.

### Usage

```bash
BASIC.py [-h] [-i TYPE] [-p NUM_THREADS] [-n NAME] [-SE FASTQ]
                [-PE_1 LEFT] [-PE_2 RIGHT] [-g GENOME] [-b BOWTIE] [-t TMPDIR]
                [-o OUTPUT_LOCATION] [-a] [-v] [-d DATABASE_PATH]
                [--subsample SUBSAMPLE] [--version]
```

```
  -h, --help            show this help message and exit
  -i TYPE               Type of receptor. Choices: "BCR" or "TCR" (default:
                        BCR)
  -p NUM_THREADS        Launch p > 2 threads that will run on separate
                        processors/cores (default: 2)
  -n NAME               Name of output file. Note: a ".fasta" extension will
                        be added (default: result.fasta)
  -SE FASTQ             Single end FASTQ file (optionally gzipped). (example:
                        se.fastq)
  -PE_1 LEFT            Paired end (left) FASTQ file (optionally gzipped).
                        -PE_2 is required and pairs must match order.
                        (example: pe_1.fastq)
  -PE_2 RIGHT           Paired end (right) FASTQ file (optionally gzipped).
                        (example: pe_2.fastq)
  -g GENOME             Options: "human" or "mouse" (default: human). Note:
                        other species are possible by adding the appropriate
                        Bowtie2 indices and following the existing db/
                        directory structure
  -b BOWTIE             Absolute path to bowtie2 executable or directory
                        containing it
  -t TMPDIR             Path to directory for writing intermediate files.
                        (default: current working directory)
  -o OUTPUT_LOCATION    Output directory (default: current working directory)
  -a                    Allow for partial reconstruction i.e. do not terminate
                        if one or both chains do not assemble.
  -v                    Turns on verbosity for more details.
  -d DATABASE_PATH      Path to database directory. Defaults to <path of
                        BASIC.py>/db
  --subsample SUBSAMPLE
                        Use the first <int> number of reads of the input.
                        Default: no limit
  --version             show program's version number and exit
```

### Paired-end reads example

```bash
BASIC.py -b <path_to_Bowtie2> -PE_1 R1.fastq.gz -PE_2 R2.fastq.gz -g human -i BCR
```

### More information
http://ttic.uchicago.edu/~aakhan/BASIC/

### Version
1.5.0 (2019/07/11)

### Publication
BASIC: BCR assembly from single cells.  
Canzar S, Neu KE, Tang Q, Wilson PC, Khan AA  
Bioinformatics, Volume 33, Issue 3, 1 February 2017.

### Contact
* Please contact Aly Azeem Khan <aakhan@ttic.edu> for any questions or comments.
* More recent updates to BASIC have been contributed by Derek Croote <dcroote@stanford.edu>

### License
Software provided to academic users under MIT License

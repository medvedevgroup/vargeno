# VarGeno
Accelerating SNP genotyping from whole genome sequencing data for bedside diagnostics

# Prerequisite
- A modern, C++11 ready compiler such as `g++` version 4.9 or higher.
- The cmake build system. (*for SDSL library. If SDSL library already installed, cmake is not needed*)
- A 64-bit operating system. Either Mac OS X or Linux are currently supported.
- 63 GB Memory for human whole genome sequencing SNP genotyping

# Quick Install

## Install SDSL

To download and install the SDSL library use the following commands.

```
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh
```
This installs the sdsl library into the `include` and `lib` directories in your home directory.


## Install VarGeno

```
git clone https://github.com/medvedevgroup/vargeno.git
cd vargeno
make all
```

You should see `vargeno`, `vqv`, `gbf` in vargeno directory.

# Quick Usage

VarGeno takes as input:
1. `ref.fa` reference genome sequence in FASTA file format.
2. list of known SNPs: `snp.txt` in [UCSC text file format](http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=varRep&hgta_track=snp141Common&hgta_table=snp141Common&hgta_doSchema=describe+table+schema).
> *VCF format support coming soon*
3. `reads.fq` sequencing reads from donor genome in FASTQ file format. If you have multiple FASTQ files, please `cat` them into one file.

> *Note here `ref.fa`, `snp.txt` and `reads.fq` are examples of filename. Your files do not need to be renamed as these.*

## Generate Dictionaries and Bloom filters

1. generate dictionaries `ref.dict`, `snp.dict` and also chrlens file `ref.fa.chrlens`

```vargeno ucscd ref.fa snp.txt ref.dict snp.dict```

2. generate Bloom filters `ref.bf` and `snp.dict`

```gbf ucsc ref.fa snp.txt ref.bf snp.bf```

## Use VarGeno
> *Note: VarGeno-QV is 3x faster than VarGeno, with slight accuracy decrease(0.04% in our experiment). Unless you have strict requirement of accuracy, we recommend using VarGeno-QV instead of VarGeno for quick variant genotyping.*

```vargeno geno ref.dict snp.dict reads.fq ref.fa.chrlens ref.bf snp.bf vargeno.out```

Here `vargeno.out` contains the variant genotyping result

## Use VarGeno-QV

```vqv geno ref.dict snp.dict reads.fq ref.fa.chrlens ref.bf snp.bf vqv.out```

Here `vqv.out` contains the variant genotyping result

## Output format

VarGeno variant genotyping output files contains 4 fields splited by `tab`:

  1. chromosome id
  2. genome position, 1-based. Here the first two fields together can uniquely identify a SNP in `snp.vcf`
  3. genotypes: `0/0`, `0/1` or `1/1` 
  4. quality score in [0,1], higher quality score means more confident genotyping result

# Experiments

To repeat experiments in our paper, please follow the instructions in [experiment.md](https://github.com/medvedevgroup/vargeno/blob/master/experiment/experiment.md) .

# Test VarGeno

To run VarGeno on test dataset, please follow the instructions in [test.md](https://github.com/medvedevgroup/vargeno/blob/master/test/test.md) .



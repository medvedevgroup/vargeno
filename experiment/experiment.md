# Experiment Guide

This file contains information to repeat experiments in our paper using VarGeno.

# Code Preparation

Please follow the instructions [here](https://github.com/medvedevgroup/vargeno) to install VarGeno.

Please make sure VarGeno properly installed before continue.

Suppose VarGeno is installed in directory `$VARGENO`.

# Data Preparation
1. reference genome sequence in FASTA format
```
wget http://cb.csail.mit.edu/cb/lava/data/hg19.fa.gz
gunzip hg19.fa.gz	
```

2. reads file in FASTQ format
```
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/NA12878/sequence_read/SRR622461_1.filt.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/NA12878/sequence_read/SRR622461_2.filt.fastq.gz
gunzip SRR622461_1.filt.fastq.gz
gunzip SRR622461_2.filt.fastq.gz
cat SRR622461_1.filt.fastq SRR622461_2.filt.fastq > reads.fq
```

3. known SNP list file in UCSC format
```
wget http://cb.csail.mit.edu/cb/lava/data/SNPs142_hg19_Common.filt.txt
mv SNPs142_hg19_Common.filt.txt snp.txt
```

# Preprocessing: generate dictionaries and Bloom filters

1. generate dictionaries `ref.dict`, `snp.dict` and also chrlens file `ref.fa.chrlens`

```$VARGENO/vargeno ucscd hg19.fa snp.txt ref.dict snp.dict```

2. generate Bloom filters `ref.bf` and `snp.dict`

```$VARGENO/gbf ucsc hg19.fa snp.txt ref.bf snp.bf```

# Genotyping

## Using VarGeno

```$VARGENO/vargeno geno ref.dict snp.dict reads.fq ref.fa.chrlens ref.bf snp.bf vargeno.out```
> **Do not forget the command `geno` after `vargeno`**

## Using VarGeno-QV

```$VARGENO/vqv geno ref.dict snp.dict reads.fq ref.fa.chrlens ref.bf snp.bf vqv.out```
> **Do not forget the command `geno` after `vqv`**

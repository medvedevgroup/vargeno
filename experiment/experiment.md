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

2. 6X reads file in FASTQ format
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

1. generate index files `ref.dict`, `snp.dict`, `ref.bf` and `snp.bf` and also chrlens file `ref.fa.chrlens`

```$VARGENO/vargeno index hg19.fa snp.txt experiment```

# Genotyping

## Using VarGeno

```$VARGENO/vargeno geno experiment reads.fq ref.fa.chrlens genotype.out```
> **Do not forget the command `geno` after `vargeno`**

# Other data

1. Affy SNP list

You can download Affy SNP list with command:
```
wget http://cb.csail.mit.edu/cb/lava/data/Affymetrix_6_SNPs.txt
```

2. High coverage data

The high coverage read data (15X, 25X and 51X) in our experiments is downloaded from GIAB FTP(ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/).

Please read the README file in above GIAB FTP link before download.

You can download the read data in our experiments with commands :

[Warning] the data is very huge (~500 GB), make sure you have enough free disk space.

```
wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_005_BH814YADXX/

wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_006_AH81VLADXX/
```

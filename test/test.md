# Test VarGeno

If you want to test if VarGeno is properly installed, or if you want to test each function module of VarGeno, please follow this manual.

This test tries to genotype 100 SNPs on human chromosome 22 with a small subset of 1000 Genome Project Illumina sequencing reads.

The test process is expected to finish within one minute(including preprossing and genotyping, not counting VarGeno installation).

The test process requires 34 GB memory. 

# Code Preparation

Please follow the instructions [here](https://github.com/medvedevgroup/vargeno) to install VarGeno.

Please make sure VarGeno properly installed before continue.

Suppose VarGeno is installed in directory `$VARGENO`.

# Go to data directory

```cd $VARGENO/test```

# Data Preparation
1. download chromosome 22 genome sequence (11MB)
```
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr22.fa.gz
gunzip chr22.fa.gz
mv chr22.fa > ref.fa	
```

2. download reads file(60MB)
```
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/NA12878/sequence_read/SRR622461.filt.fastq.gz
gunzip SRR622461.filt.fastq.gz
mv SRR622461.filt.fastq reads.fq
```

3. known SNP list file in UCSC format
The SNP list file `snp.txt` is already in `$VARGENO/test`.

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

## Output format

VarGeno variant genotyping output files contains 4 fields splited by `tab`:

  1. chromosome id
  2. genome position, 1-based. Here the first two fields together can uniquely identify a SNP in `snp.vcf`
  3. genotypes: `0/0`, `0/1` or `1/1` 
  4. quality score in [0,1], higher quality score means more confident genotyping result
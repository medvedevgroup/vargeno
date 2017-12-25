# VarGeno
Accelerating SNP genotyping from whole genome sequencing data for bedside diagnostics

# Prerequisites
- A modern, C++11 ready compiler, such as `g++` version 4.9 or higher.
- The cmake build system (*only necessary to install SDSL library. If SDSL library already installed, cmake is not needed*)
- A 64-bit operating system. Either Mac OS X or Linux are currently supported.

# Quick Install

VarGeno requires the SDSL library. If you do not have it already installed, please use the following commands.

```
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh
```
This installs the sdsl library into the `include` and `lib` directories in your home directory.

To install VarGeno itself,

```
git clone https://github.com/medvedevgroup/vargeno.git
cd vargeno
make all
```

You should then see `vargeno`, `vqv`, `gbf` in vargeno directory. To verify that your installation is correct, you can run the toy example below. 


# Quick Usage

VarGeno takes as input:
1. A reference genome sequence in FASTA file format.
2. A list of SNPs to be genotyped, in [UCSC text file format](http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=varRep&hgta_track=snp141Common&hgta_table=snp141Common&hgta_doSchema=describe+table+schema). VCF format support coming soon
3. Sequencing reads from the donor genome in FASTQ file format. If you have multiple FASTQ files, please `cat` them into one file.

Before genotyping an individual, you must construct indices for the reference using the following commands:
```
vargeno ucscd ref.fa snp.txt ref.dict snp.dict
gbf ucsc ref.fa snp.txt ref.bf snp.bf
```
This constructs the reference dictionaries `ref.dict` and `snp.dict`, the reference Bloom filters `ref.bf` and `snp.bf`, and also a file with the chromosome lengths `ref.fa.chrlens`.

To perform the genotyping, you can use either VarGeno-QV:
```
vqv geno ref.dict snp.dict reads.fq ref.fa.chrlens ref.bf snp.bf result.out
```
or VarGeno:
```
vargeno geno ref.dict snp.dict reads.fq ref.fa.chrlens ref.bf snp.bf result.out
```

VarGeno-QV is 3x faster than VarGeno, with only a slight decrease in accuracy (0.04% in our experiments). Unless you have strict requirement for accuracy, we recommend using VarGeno-QV instead of VarGeno 


## Output format

VarGeno variant genotyping output files contains 4 tab-separated fields for each SNP: 

  1. chromosome id
  2. genome position (1-based): The first two fields together uniquely identify a SNP in the input SNP list.
  3. genotypes: `0/0`, `0/1` or `1/1` 
  4. quality score in [0,1]: higher quality score means more confident genotyping result

# Example

In this example, we genotype 100 SNPs on human chromosome 22 with a small subset of 1000 Genome Project Illumina sequencing reads. The whole process should finish in around a minute and requries 34 GB RAM. Suppose VarGeno is installed in directory `$VARGENO`.

1. go to VarGeno directory
```
cd $VARGENO/test
```

2. pre-process the reference and SNP list to generate indices:
```
$VARGENO/vargeno ucscd chr22.fa snp.txt ref.dict snp.dict
$VARGENO/gbf ucsc chr22.fa snp.txt ref.bf snp.bf
```

3. genotype variants:
```
$VARGENO/vqv geno ref.dict snp.dict reads.fq ref.fa.chrlens ref.bf snp.bf result.out
```


# Experiments in paper

To repeat the experiments in our paper, please follow the instructions in [experiment.md](https://github.com/medvedevgroup/vargeno/blob/master/experiment/experiment.md) .

# Citation

If you use VarGeno in your research, please cite
* Chen Sun and Paul Medvedev, Accelerating SNP genotyping from whole genome sequencing data for bedside diagnostics

VarGeno's algorithm is built on top of LAVA's. Its code is built on top of LAVA's and it reuses a lot of LAVA's code. It uses some code from the [AllSome project](https://github.com/medvedevgroup/bloomtree-allsome).
* Shajii A, Yorukoglu D, William Yu Y, Berger B, [Fast genotyping of known SNPs through approximate k-mer matching,](https://academic.oup.com/bioinformatics/article/32/17/i538/2450790) Bioinformatics. 2016 32(17):i538-i544. Code is available [here](https://github.com/arshajii/lava/).

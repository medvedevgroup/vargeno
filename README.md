# VarGeno
Fase SNP genotyping tool for whole genome sequencing data and large SNP database.

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
You should then see `vargeno`, `gbf` in vargeno directory. To verify that your installation is correct, you can run the toy example below. 

To install VarGeno Lite version:
```
cd vargeno/vargeno_lite
make all
```
You should then see `vargeno_lite`, `gbf_lite` in vargeno/vargeno_lite directory.


# Quick Usage

VarGeno takes as input:
1. A reference genome sequence in FASTA file format.
2. A list of SNPs to be genotyped. VCF file format and [UCSC text file format](http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=varRep&hgta_track=snp141Common&hgta_table=snp141Common&hgta_doSchema=describe+table+schema) are supported.
3. Sequencing reads from the donor genome in FASTQ file format.

Before genotyping an individual, you must construct indices for the reference and SNP list using the following commands:
```
vargeno index ref.fa snp.txt index_prefix
```

If your SNP list to be genotyped is in VCF file format, the file extension should be `.vcf`. If the SNP list is in UCSC text file format, the file extension should be `.txt`.

To perform the genotyping:
```
vargeno geno index_prefix reads.fq ref.fa.chrlens output_filename
```

Here `index_prefix` should be the same string as index generating.

The use of VarGeno-Lite is similar, for detail usage and example, please refer to the README in `vargeno/vargeno_lite`.

## Output format

VarGeno variant genotyping output files contains 4 tab-separated fields for each SNP: 

  1. chromosome id
  2. genome position (1-based): The first two fields together uniquely identify a SNP in the input SNP list.
  3. genotypes: `0/0`, `0/1` or `1/1` 
  4. quality score in [0,1]: higher quality score means more confident genotyping result

# Example

In this example, we genotype 100 SNPs on human chromosome 22 with a small subset of 1000 Genome Project Illumina sequencing reads. The whole process should finish in around a minute and requries 34 GB RAM. Suppose VarGeno is installed in directory `$VARGENO`.

1. go to test data directory
```
cd $VARGENO/test
```

2. pre-process the reference and SNP list to generate indices:
```
$VARGENO/vargeno index chr22.fa snp.txt test_prefix
```

3. genotype variants:
```
$VARGENO/vargeno geno test_prefix reads.fq chr22.fa.chrlens genotype.out
```

# Experiments in paper

To repeat the experiments in our paper, please follow the instructions in [experiment.md](https://github.com/medvedevgroup/vargeno/blob/master/experiment/experiment.md) .

# Citation

If you use VarGeno in your research, please cite
* Chen Sun and Paul Medvedev, Toward fast and accurate SNP genotyping from whole genome sequencing data for bedside diagnostics.

VarGeno's algorithm is built on top of LAVA's. Its code is built on top of LAVA's and it reuses a lot of LAVA's code. It uses some code from the [AllSome project](https://github.com/medvedevgroup/bloomtree-allsome).
* Shajii A, Yorukoglu D, William Yu Y, Berger B, [Fast genotyping of known SNPs through approximate k-mer matching,](https://academic.oup.com/bioinformatics/article/32/17/i538/2450790) Bioinformatics. 2016 32(17):i538-i544. Code is available [here](https://github.com/arshajii/lava/).

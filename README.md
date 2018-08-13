# VarGeno
Fase SNP genotyping tool for whole genome sequencing data and large SNP database.

# Prerequisites
- A modern, C++11 ready compiler, such as `g++` version 4.9 or higher.
- The cmake build system (*only necessary to install SDSL library. If SDSL library already installed, cmake is not needed*)
- A 64-bit operating system. Either Mac OS X or Linux are currently supported.

# Quick Install
```
git clone https://github.com/medvedevgroup/vargeno.git
cd vargeno
export PREFIX=$HOME
bash ./install.sh
```
You should then see `vargeno` in vargeno directory. To verify that your installation is correct, you can run the toy example below. 

# Quick Usage

VarGeno takes as input:
1. A reference genome sequence in FASTA file format.
2. A list of SNPs to be genotyped in VCF file format.
3. Sequencing reads from the donor genome in FASTQ file format.

Before genotyping an individual, you must construct indices for the reference and SNP list using the following commands:
```
vargeno index ref.fa snp.vcf index_prefix
```

To perform the genotyping:
```
vargeno geno index_prefix reads.fq snp.vcf output_filename
```

Here `index_prefix` should be the same string as index generating.

## Output format: VCF

VarGeno's genotyping results are in the "FORMAT" column of VCF file.

  1. genotypes: in "GT" field: `0/0`, `0/1` or `1/1`.
  2. genotype quality: in "GQ" field, encoded as a phred quality (Integer).

For details of "GT" and "GQ" fields, please refer to [The Variant Call Format(VCF) Version 4.2 Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

# Example

In this example, we genotype 100 SNPs on human chromosome 22 with a small subset of 1000 Genome Project Illumina sequencing reads. The whole process should finish in around a minute and requries 34 GB RAM. Suppose VarGeno is installed in directory `$VARGENO`.

1. go to test data directory
```
cd $VARGENO/test
```

2. pre-process the reference and SNP list to generate indices:
```
$VARGENO/vargeno index chr22.fa snp.vcf test_prefix
```

3. genotype variants:
```
$VARGENO/vargeno geno test_prefix reads.fq snp.vcf genotyped.vcf
```

# Memory Lite Version

The memory lite version of VarGeno (VarGeno-Lite) is maintained as an independent project in [https://github.com/medvedevgroup/vargeno_lite](https://github.com/medvedevgroup/vargeno_lite).

# Citation

If you use VarGeno in your research, please cite
* Chen Sun and Paul Medvedev, Toward fast and accurate SNP genotyping from whole genome sequencing data for bedside diagnostics.

VarGeno's algorithm is built on top of LAVA's. Its code is built on top of LAVA's and it reuses a lot of LAVA's code. It uses some code from the [AllSome project](https://github.com/medvedevgroup/bloomtree-allsome).
* Shajii A, Yorukoglu D, William Yu Y, Berger B, [Fast genotyping of known SNPs through approximate k-mer matching,](https://academic.oup.com/bioinformatics/article/32/17/i538/2450790) Bioinformatics. 2016 32(17):i538-i544. Code is available [here](https://github.com/arshajii/lava/).

# HiMBar:  Multi-kingdom Barcode Extractor 
![Good](https://img.shields.io/badge/latest%20version-v1.1.0-red) ![Active](https://www.repostatus.org/badges/latest/active.svg) ![GPL](https://img.shields.io/badge/license-GPLv3.0-blue)

## Introduction
### `HiMBar.sh` is an automated pipeline for extracting and annotating barcoding sequences (rDNA, CO1, rbcL) from PacBio HiFi metagenomic data. It performs:
1. CCS error correction (optional if FASTA input is provided)
2. Barcode extraction using HMMER (rDNA) and BLAST (CO1/rbcL)
3. VSEARCH clustering and taxonomy annotation (SILVA + PR2)
4. Phylogenetic tree construction (MAFFT + Gblocks + RAxML)


<p align="center">
  <img src="https://github.com/user-attachments/assets/c2fcb444-8f93-4072-9cd1-c7ded08bc650" width="600">
</p>




## Dependencies 
- bash
- Perl
- hmmer (v3.3+)
- blast+ (v2.11+)
- seqtk
- Python (optional)

## Install
We recommend using Conda:
```bash
conda create -n HIFI_barcoding hmmer blast seqkit vsearch mafft gblocks raxml -y
conda activate HIFI_barcoding
```

```
$ git clone https://github.com/Xchenkai2019/HIFI_barcoding.git
# give executable permission to all scripts in HiMBar scripts directory
$ chmod a+x HIFI_barcoding/HiMBar.sh
# add HiMBar scripts directory to your PATH environment variable
$ echo 'PATH=$(pwd)/HIFI_barcoding/:$PATH' >> ~/.bashrc
$ source ~/.bashrc
```
## Usage

```bash
bash HiMBar.sh \
  --input test_data/subreads.ccs.bam \
  --outdir results \
  --hmm-rdna db/all.rDNA.hmm \
  --db-co1   db/PF00115_full.txt.400up.550down.fa \
  --db-rbcl  db/rbcl_complete_cds.fasta \
  --silva    db/SILVA_138.1_SSURef_NR99_tax_silva.fasta \
  --pr2      db/pr2_version_4.14.0_SSU_UTAX.fasta \
  --threads 4 --min-passes 3 --min-rq 0.99

If the input is already in FASTA format (rather than BAM), the CCS step can be skipped：
bash HiMBar.sh \
  --input test_data/subreads.ccs.bam.fa --skip-ccs \
  --outdir results \
  --hmm-rdna db/all.rDNA.hmm \
  --db-co1   db/PF00115_full.txt.400up.550down.fa \
  --db-rbcl  db/rbcl_complete_cds.fasta \
  --silva    db/SILVA_138.1_SSURef_NR99_tax_silva.fasta \
  --pr2      db/pr2_version_4.14.0_SSU_UTAX.fasta

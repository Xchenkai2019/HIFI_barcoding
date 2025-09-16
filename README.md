# HiMBar:  Multi-kingdom Barcode Extractor 
![Good](https://img.shields.io/badge/latest%20version-v1.1.0-red) ![Active](https://www.repostatus.org/badges/latest/active.svg) ![GPL](https://img.shields.io/badge/license-GPLv3.0-blue)

## Introduction
### HiMBar is an automated pipeline to extract taxonomic barcodes from PacBio CCS reads, supporting multiple  kingdoms (e.g., microbial SSU, algal rbcl, protist markers).


<p align="center">
  <img src="https://github.com/user-attachments/assets/c2fcb444-8f93-4072-9cd1-c7ded08bc650" width="600">
</p>



## Features 
- Supports HMM scanning for rDNA barcoding extraction

- Supports marker gene alignment for Mitochondria and chloroplasts barcoding extractio (e.g., CO1, rbcL)

- Supports parameterized control (e.g., --min-passes, --min-rq)

## Dependencies 
- bash
- Perl
- hmmer (v3.3+)
- blast+ (v2.11+)
- seqtk
- Python (optional)

## Usage

```bash
bash HiMBar.sh \
  --input test_data/subreads.ccs.bam.fa \
  --hmm_db db/all.rDNA.hmm \
  --rbcl_db db/rbcl_complete_cds.fasta \
  --pfam_db db/PF00115_full.txt.400up.550down.fa \
  --min-passes 3 \
  --min-rq 0.99 \
  --num-threads 4 \
  --outdir output

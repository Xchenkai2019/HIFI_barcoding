# HiMBar:  Multi-kingdom Barcode Extractor 

## Introduction
![六界系统条形码](https://github.com/user-attachments/assets/c69d1558-af04-4647-829a-0ce937f7c2fc)

HiMBar is an automated pipeline to extract taxonomic barcodes from PacBio CCS reads, supporting multiple  kingdoms (e.g., microbial SSU, algal rbcl, protist markers).


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

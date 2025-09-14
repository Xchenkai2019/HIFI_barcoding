# HiMBar:  Multi-kingdom Barcode Extractor | 多生物界条形码提取工具

## Introduction | 简介

HiMBar is an automated pipeline to extract taxonomic barcodes from PacBio CCS reads, supporting multiple  kingdoms (e.g., microbial SSU, algal rbcl, protist markers).

HiMBar 是一个用于从三代测序数据（如 PacBio CCS）中提取多生物界条形码（如 SSU rRNA, rbcl 等）的自动化工具。适用于微生物、藻类、植物等不同界别的群落结构分析。

## Features | 功能 
- 支持标准HMM扫描提取rDNA
- 支持功能基因比对（rbcl, PF00115）
- 支持参数化控制 (--min-passes, --min-rq)
- 模块化结构，自动日志记录
- 支持全流程中文与英文注释

## Dependencies | 安装依赖 
- bash
- Perl
- hmmer (v3.3+)
- blast+ (v2.11+)
- seqtk
- Python (optional)

## Usage | 用法

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

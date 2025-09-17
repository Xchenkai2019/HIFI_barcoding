#!/usr/bin/env bash
# HiMBar.sh ¡ª Automated HiFi barcoding pipeline
# Steps:
# 1) CCS error-correction and FASTA export
# 2) rDNA / CO1 / rbcL barcode extraction
# 3) Taxonomy assignment (rDNA with SILVA+PR2; CO1/rbcL with NCBI or user DB)
# 4) Phylogeny (MAFFT -> Gblocks -> RAxML)
set -euo pipefail

###############
# Usage/Help  #
###############
usage() {
  cat <<'EOF'
Usage:
  bash HiMBar.sh \
    --input subreads.bam \
    --outdir results \
    --hmm-rdna db/all.rDNA.hmm \
    --db-co1 db/PF00115_full.txt.400up.550down.fa \
    --db-rbcl db/rbcl_complete_cds.fasta \
    --silva db/SILVA_138.1_SSURef_NR99_tax_silva.fasta \
    --pr2   db/pr2_version_4.14.0_SSU_UTAX.fasta \
    [--threads 16] [--min-passes 3] [--min-rq 0.99] [--skip-ccs] [--keep-temp]

Required:
  --input       Input PacBio subreads BAM (or FASTA if --skip-ccs is used)
  --outdir      Output directory
  --hmm-rdna    rDNA HMM database
  --db-co1      CO1 reference FASTA
  --db-rbcl     rbcL reference FASTA
  --silva       SILVA SSU reference FASTA
  --pr2         PR2 SSU reference FASTA

Optional:
  --threads     Number of CPU threads (default: 4)
  --min-passes  Minimum CCS passes (default: 3)
  --min-rq      Minimum CCS read quality (default: 0.99)
  --skip-ccs    Skip CCS if input is already FASTA
  --keep-temp   Keep intermediate files for debugging
  --help        Print this help message

Outputs (under --outdir):
	01_ccs/         CCS reads and FASTA files
	02_barcode/     Extracted barcodes for rDNA, CO1, and rbcL (FASTA)
	03_taxonomy/    Annotation results (BLAST/VSEARCH tables, top-hits, integrated table)
	04_phylogeny/   Alignments, filtering, and phylogenetic tree for rDNA representative sequences


EOF
  exit 1
}

#####################
# Default parameters #
#####################
threads=4
min_passes=3
min_rq=0.99
skip_ccs=0
keep_temp=0

input=""
outdir=""
hmm_rdna=""
db_co1=""
db_rbcl=""
db_silva=""
db_pr2=""

#####################
# Parse arguments   #
#####################
[[ $# -eq 0 ]] && usage
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)       input="$2"; shift 2;;
    --outdir)      outdir="$2"; shift 2;;
    --hmm-rdna)    hmm_rdna="$2"; shift 2;;
    --db-co1)      db_co1="$2"; shift 2;;
    --db-rbcl)     db_rbcl="$2"; shift 2;;
    --silva)       db_silva="$2"; shift 2;;
    --pr2)         db_pr2="$2"; shift 2;;
    --threads)     threads="$2"; shift 2;;
    --min-passes)  min_passes="$2"; shift 2;;
    --min-rq)      min_rq="$2"; shift 2;;
    --skip-ccs)    skip_ccs=1; shift 1;;
    --keep-temp)   keep_temp=1; shift 1;;
    --help)        usage;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

#####################
# Sanity checks     #
#####################
err() { echo "[ERROR] $*" >&2; exit 1; }

[[ -z "$input" || -z "$outdir" || -z "$hmm_rdna" || -z "$db_co1" || -z "$db_rbcl" || -z "$db_silva" || -z "$db_pr2" ]] && usage
[[ -f "$input" ]]     || err "Input not found: $input"
[[ -f "$hmm_rdna" ]]  || err "rDNA HMM not found: $hmm_rdna"
[[ -f "$db_co1" ]]    || err "CO1 DB not found: $db_co1"
[[ -f "$db_rbcl" ]]   || err "rbcL DB not found: $db_rbcl"
[[ -f "$db_silva" ]]  || err "SILVA DB not found: $db_silva"
[[ -f "$db_pr2" ]]    || err "PR2 DB not found: $db_pr2"

# tool deps
need_tools=(nhmmer hmmscan blastn seqkit vsearch mafft Gblocks)
# CCS and bam2fastx only needed when not skipping CCS
if [[ $skip_ccs -eq 0 ]]; then need_tools+=(ccs bam2fastx); fi

for t in "${need_tools[@]}"; do
  command -v "$t" >/dev/null 2>&1 || err "Tool not in PATH: $t"
done

############
# Logging  #
############
mkdir -p "$outdir"
log="$outdir/HiMBar_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$log") 2>&1

echo "[INFO] Start HiMBar pipeline"
echo "[INFO] threads=$threads min_passes=$min_passes min_rq=$min_rq skip_ccs=$skip_ccs keep_temp=$keep_temp"

# layout
d01="$outdir/01_ccs"
d02="$outdir/02_barcode"
d03="$outdir/03_taxonomy"
d04="$outdir/04_phylogeny"
mkdir -p "$d01" "$d02" "$d03" "$d04"

################################
# Step1. CCS & FASTA exporting #
################################
fasta="$d01/ccs_reads.fasta"
if [[ $skip_ccs -eq 0 ]]; then
  echo "[STEP1] CCS correction and export FASTA"
  ccs "$input" "$d01/subreads.ccs.bam" --min-passes "$min_passes" --min-rq "$min_rq" --num-threads "$threads"
  bam2fastx -a -A -o "$fasta" "$d01/subreads.ccs.bam"
else
  echo "[STEP1] Skip CCS; treat --input as FASTA"
  # link/copy to standard name
  cp -f "$input" "$fasta"
fi
[[ -s "$fasta" ]] || err "No sequences found in $fasta"

####################################################
# Step2. Extract barcodes (rDNA, CO1, rbcL)       #
####################################################
echo "[STEP2] Extracting barcodes"

# 2.1 rDNA using nhmmer (DNA HMMs) -> domains table -> regions -> subseq
rdna_tbl="$d02/rdna.dom.tbl"
nhmmer --cpu "$threads" --tblout "$rdna_tbl" "$hmm_rdna" "$fasta" >/dev/null

# helper function: from nhmmer --tblout, pick best domain per query per target-key (SSU/LSU/5.8S/5S)
pick_best_regions() {
  local key="$1" out_file="$2"
  # nhmmer --tblout columns are comment-heavy; skip '#' lines
  # We'll extract qname, strand, qstart, qend heuristically from the "alignment numbering" columns.
  # Many nhmmer table formats: safer approach is to re-run nhmmer with --dfamtblout? For now, simple filter by $2 (target name contains key)
  awk -v K="$key" '
    $0 !~ /^#/ && index($2,K)>0 {
      q=$1; s=$13; qs=$17; qe=$18;     # qs,qe are approximate: check your nhmmer version
      if(qs=="" || qe==""){next}
      if(s=="+"){print q"\t"qs"\t"qe}
      else      {print q"\t"qe"\t"qs}
    }' "$rdna_tbl" \
  | sort -k1,1 -k2,2n \
  | awk '!seen[$1]++' > "$out_file" || true
}

# pick SSU/LSU/5.8S/5S regions (best per query)
pick_best_regions "SSU"  "$d02/SSU.regions"
pick_best_regions "LSU"  "$d02/LSU.regions"
pick_best_regions "5.8S" "$d02/5_8S.regions"
pick_best_regions "5S"   "$d02/5S.regions"

# extract with seqkit subseq; reverse-complement handled by regions choice above
[[ -s "$d02/SSU.regions"  ]] && seqkit subseq "$fasta" --bed <(awk '{print $1"\t"$2-1"\t"$3}' "$d02/SSU.regions")  > "$d02/SSU.fasta"  || true
[[ -s "$d02/LSU.regions"  ]] && seqkit subseq "$fasta" --bed <(awk '{print $1"\t"$2-1"\t"$3}' "$d02/LSU.regions")  > "$d02/LSU.fasta"  || true
[[ -s "$d02/5_8S.regions" ]] && seqkit subseq "$fasta" --bed <(awk '{print $1"\t"$2-1"\t"$3}' "$d02/5_8S.regions") > "$d02/5.8S.fasta" || true
[[ -s "$d02/5S.regions"   ]] && seqkit subseq "$fasta" --bed <(awk '{print $1"\t"$2-1"\t"$3}' "$d02/5S.regions")   > "$d02/5S.fasta"   || true

# 2.2 CO1 via BLASTN to PF00115_full... (outfmt6 -> first best hit per query -> regions -> subseq)
co1_b6="$d02/CO1.blastn.tsv"
blastn -query "$fasta" -db "$db_co1" -out "$co1_b6" -outfmt 6 -num_threads "$threads" -evalue 0.01
# choose first best per query by bitscore desc
co1_plus="$d02/CO1.plus.regions"
co1_minus="$d02/CO1.minus.regions"
sort -k1,1 -k12,12gr "$co1_b6" | awk '!seen[$1]++{ if($7>$8){print $1"\t"$7"\t"$8 > "'"$co1_plus"'" } else {print $1"\t"$8"\t"$7 > "'"$co1_minus"'" }}'
[[ -s "$co1_plus"  ]] && seqkit subseq "$fasta" --bed <(awk '{print $1"\t"$2-1"\t"$3}' "$co1_plus")  > "$d02/CO1.fasta" || true
[[ -s "$co1_minus" ]] && seqkit subseq "$fasta" --bed <(awk '{print $1"\t"$2-1"\t"$3}' "$co1_minus") >> "$d02/CO1.fasta" || true

# 2.3 rbcL via BLASTN
rbcl_b6="$d02/rbcl.blastn.tsv"
blastn -query "$fasta" -db "$db_rbcl" -out "$rbcl_b6" -outfmt 6 -num_threads "$threads" -evalue 0.01
rbcl_plus="$d02/rbcl.plus.regions"
rbcl_minus="$d02/rbcl.minus.regions"
sort -k1,1 -k12,12gr "$rbcl_b6" | awk '!seen[$1]++{ if($7>$8){print $1"\t"$7"\t"$8 > "'"$rbcl_plus"'" } else {print $1"\t"$8"\t"$7 > "'"$rbcl_minus"'" }}'
[[ -s "$rbcl_plus"  ]] && seqkit subseq "$fasta" --bed <(awk '{print $1"\t"$2-1"\t"$3}' "$rbcl_plus")  > "$d02/rbcl.fasta" || true
[[ -s "$rbcl_minus" ]] && seqkit subseq "$fasta" --bed <(awk '{print $1"\t"$2-1"\t"$3}' "$rbcl_minus") >> "$d02/rbcl.fasta" || true

####################################################
# Step3. Taxonomy assignment                       #
####################################################
echo "[STEP3] Taxonomy (rDNA: VSEARCH cluster + BLAST to SILVA/PR2)"

# 3.1 rDNA (SSU as an example; LSU/5.8S can also be processed similarly if needed)
if [[ -s "$d02/SSU.fasta" ]]; then
  vsearch --cluster_size "$d02/SSU.fasta" --strand both --id 0.97 \
          --centroids "$d03/SSU.centroids97.fasta" \
          --otutabout "$d03/SSU.otu97.tsv" \
          --sizeout --threads "$threads" --uc "$d03/SSU.uc.97"
  # BLAST against SILVA & PR2
  blastn -query "$d03/SSU.centroids97.fasta" -db "$db_silva" -out "$d03/SSU.silva.b6" -outfmt 6 -evalue 0.01 -num_threads "$threads"
  blastn -query "$d03/SSU.centroids97.fasta" -db "$db_pr2"   -out "$d03/SSU.pr2.b6"   -outfmt 6 -evalue 0.01 -num_threads "$threads"

  # Get top-hit for each query (based on bitscore)
  awk 'BEGIN{OFS="\t"}{k=$1; if(!(k in s) || $12>s[k]){s[k]=$12; line[k]=$0}} END{for(k in line)print line[k]}' "$d03/SSU.silva.b6" > "$d03/SSU.silva.top"
  awk 'BEGIN{OFS="\t"}{k=$1; if(!(k in s) || $12>s[k]){s[k]=$12; line[k]=$0}} END{for(k in line)print line[k]}' "$d03/SSU.pr2.b6"   > "$d03/SSU.pr2.top"

  # Simplified merging (use sseqid from the top line as representative species name; 
  # if full taxonomy is required, add your mapping or parsing script here)
  cut -f1,2 "$d03/SSU.silva.top" > "$d03/SSU.silva.map"
  cut -f1,2 "$d03/SSU.pr2.top"   > "$d03/SSU.pr2.map"
  # Merge results: prioritize SILVA hits (use PR2 for missing entries)
  awk 'NR==FNR{a[$1]=$2; next} {if(a[$1]!=""){print $1"\t"a[$1]} else {print $1"\t"$2}}' "$d03/SSU.silva.map" "$d03/SSU.pr2.map" > "$d03/SSU.taxonomy.tsv"
fi

# 3.2 CO1 / rbcL annotation (example: using sseqid from top-hit as the label;
# if NCBI taxonomic hierarchy is needed, perform additional mapping for sseqid)
if [[ -s "$d02/CO1.fasta" ]]; then
  blastn -query "$d02/CO1.fasta" -db "$db_co1" -out "$d03/CO1.b6" -outfmt 6 -num_threads "$threads" -evalue 0.01
  awk 'BEGIN{OFS="\t"}{k=$1; if(!(k in s) || $12>s[k]){s[k]=$12; line[k]=$0}} END{for(k in line)print line[k]}' "$d03/CO1.b6" | cut -f1,2 > "$d03/CO1.taxonomy.tsv"
fi

if [[ -s "$d02/rbcl.fasta" ]]; then
  blastn -query "$d02/rbcl.fasta" -db "$db_rbcl" -out "$d03/rbcl.b6" -outfmt 6 -num_threads "$threads" -evalue 0.01
  awk 'BEGIN{OFS="\t"}{k=$1; if(!(k in s) || $12>s[k]){s[k]=$12; line[k]=$0}} END{for(k in line)print line[k]}' "$d03/rbcl.b6" | cut -f1,2 > "$d03/rbcl.taxonomy.tsv"
fi

####################################################
# Step4. Phylogeny (SSU representatives)           #
####################################################
echo "[STEP4] Phylogeny (SSU representatives)"
if [[ -s "$d03/SSU.centroids97.fasta" ]]; then
  mafft --thread "$threads" --auto "$d03/SSU.centroids97.fasta" > "$d04/SSU.centroids97.aln.fasta"
  Gblocks "$d04/SSU.centroids97.aln.fasta" -t=d -b2=0.85 -b3=8 -b4=2 -b5=a -e=.gb
  sed -i 's/ //g' "$d04/SSU.centroids97.aln.fasta.gb"
  if command -v raxmlHPC-PTHREADS >/dev/null 2>&1; then
    raxmlHPC-PTHREADS -x 1234567890 -p 1234567890 -f a -# 100 -T "$threads" \
      -m GTRGAMMA -s "$d04/SSU.centroids97.aln.fasta.gb" -n SSU.tree
  else
    raxmlHPC -x 1234567890 -p 1234567890 -f a -# 100 -m GTRGAMMA \
      -s "$d04/SSU.centroids97.aln.fasta.gb" -n SSU.tree
  fi
fi

####################
# Cleanup & Finish #
####################
if [[ $keep_temp -eq 0 ]]; then
  # Intermediate large files can be deleted here; by default, they are retained for traceability.
  :
fi

echo "[DONE] Pipeline completed. See outputs in: $outdir"

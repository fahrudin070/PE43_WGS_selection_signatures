#!/bin/bash
# Script: 08_snpeff_annotation.sh
# Purpose: Annotate filtered PE43 variants using snpEff against the ARS1 goat genome database.
# Input:
#   - filtered PE43 autosomal SNP VCF
#   - snpEff ARS1 database and config
# Output:
#   - annotated VCF (.vcf.gz)
#   - snpEff HTML summary report
# Main tools:
#   - snpEff
#   - bgzip
#   - tabix
# Notes:
#   - Final functional annotation step for PE43 variants
#SBATCH --job-name=snpeff_PE43_ARS1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=snpeff_PE43_ARS1.%j.out
#SBATCH --error=snpeff_PE43_ARS1.%j.err
set -euo pipefail

module load 2024

# snpEff + java
export PATH=/lustre/nobackup/WUR/ABGC/husen002/env_snpeff/bin:$PATH

# bgzip/tabix/bcftools (ambil dari env pipeline yang biasa kamu pakai)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /lustre/nobackup/WUR/ABGC/shared/pipelines_version2/envs/population-variant-calling-version2

WORK=/lustre/nobackup/WUR/ABGC/husen002/pe43_popvar
VCALL=$WORK/results/results/variant_calling
VCF=$VCALL/final_qc/PE43.Q20DP5.autosome.MAF05.vcf.gz

SNPEFF_BASE=$WORK/data/snpeff/ARS1
SNPEFF_DATA=$SNPEFF_BASE/data
GENOME=Capra_hircus_ARS1
CFG=$SNPEFF_BASE/snpEff.config

OUTDIR=$VCALL/annotation_final/snpeff_ARS1
mkdir -p $OUTDIR

echo "snpEff:"; which snpEff
echo "bgzip:"; which bgzip
echo "tabix:"; which tabix

# IMPORTANT: matikan download otomatis
snpEff -Xmx50g -c $CFG -dataDir $SNPEFF_DATA -v -nodownload \
  -stats $OUTDIR/PE43.Q20DP5.MAF05.snpeff.html \
  $GENOME $VCF | bgzip -c > $OUTDIR/PE43.Q20DP5.MAF05.snpeff.vcf.gz

tabix -f -p vcf $OUTDIR/PE43.Q20DP5.MAF05.snpeff.vcf.gz

#!/usr/bin/env bash
# Script: 07_beagle_imputation.sh
# Purpose: Perform chromosome-wise genotype imputation using Beagle.
# Input:
#   - chromosome-specific VCF files
#   - Beagle JAR file
# Output:
#   - imputed VCFs per chromosome
# Main tools:
#   - Beagle
#   - Java
#   - SLURM array jobs
# Notes:
#   - Imputation step for multi-breed ARS1-ready VCF dataset
#SBATCH --job-name=beagle29
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --array=1-29
#SBATCH --output=logs/beagle_%A_%a.out
#SBATCH --error=logs/beagle_%A_%a.err

set -eo pipefail

BASE=/lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/data/external_breeds/ars1_ready/final_vcfs
cd "$BASE"

BGL=/lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/data/external_breeds/pca/selection/tools/beagle.jar
IN=02_impute/by_chr/diploid_fixed_clean
OUT=02_impute/by_chr/imputed_fromSRC


mkdir -p logs "$OUT"

# ambil chr ke-N sesuai urutan sort -V
chr=$(ls -1 $IN/NC_*.vcf.gz | sort -V | sed -n "${SLURM_ARRAY_TASK_ID}p" | xargs -n1 basename | sed 's/\.vcf\.gz$//')

echo "CHR=$chr"
echo "INPUT=$IN/${chr}.vcf.gz"

java -Xmx55g -jar "$BGL" \
  gt="$IN/${chr}.vcf.gz" \
  out="$OUT/${chr}.imputed" \
  nthreads=8

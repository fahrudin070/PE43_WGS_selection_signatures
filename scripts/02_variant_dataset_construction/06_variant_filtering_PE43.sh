#!/bin/bash
# Script: 06_variant_filtering_PE43.sh
# Purpose: Filter raw PE43 VCF to retain high-quality biallelic SNPs using quality, depth, missingness, and MAF thresholds.
# Input:
#   - PE43_EBV.vcf.gz
# Output:
#   - filtered SNP VCF files
#   - indexed filtered VCF files
# Main tools:
#   - bcftools
# Notes:
#   - Final SNP filtering step for PE43 downstream analyses
#SBATCH --job-name=PE43_filter
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=PE43_filter_%j.out
#SBATCH --error=PE43_filter_%j.err

set -euo pipefail

# --- Activate conda environment (FULL PATH, jangan kepotong) ---
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate /lustre/nobackup/WUR/ABGC/shared/pipelines_version2/envs/population-variant-calling-version2

# --- Working directory (FULL PATH, jangan kepotong) ---
WORKDIR="/lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/results/results/variant_calling"
cd "$WORKDIR"

# --- Sanity checks (biar kalau salah path langsung ketahuan) ---
echo "PWD: $(pwd)"
echo "Using bcftools: $(which bcftools)"
ls -lh PE43_EBV.vcf.gz PE43_EBV.vcf.gz.tbi

echo "Step 1: Biallelic SNP + QUAL>=30 + INFO/DP>=5"
bcftools view -m2 -M2 -v snps -Ou PE43_EBV.vcf.gz \
| bcftools filter -i 'QUAL>=30 && INFO/DP>=5' -Ou \
| bcftools sort -Ou \
| bcftools view -Oz -o PE43_EBV.Q30.DP5.snp.vcf.gz

bcftools index -f PE43_EBV.Q30.DP5.snp.vcf.gz

echo "Step 2: Add tags (MAF, F_MISSING) then filter MAF>=0.05 & Missing<=10%"
bcftools +fill-tags PE43_EBV.Q30.DP5.snp.vcf.gz -Ou -- -t MAF,F_MISSING \
| bcftools filter -e 'F_MISSING > 0.10 || MAF < 0.05' -Oz -o PE43_EBV.Q30.DP5.MAF05.MISS10.snp.vcf.gz

bcftools index -f PE43_EBV.Q30.DP5.MAF05.MISS10.snp.vcf.gz

echo "Filtering finished successfully"
echo "Outputs:"
ls -lh PE43_EBV.Q30.DP5.snp.vcf.gz* PE43_EBV.Q30.DP5.MAF05.MISS10.snp.vcf.gz*


#!/bin/bash
#SBATCH --job-name=PCA370
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/PCA370_%j.out
#SBATCH --error=logs/PCA370_%j.err

set -eo pipefail   # <-- JANGAN pakai -u sebelum conda activate

cd /lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/data/external_breeds/ars1_ready/final_vcfs
mkdir -p logs

date
echo "HOST=$(hostname)"
echo "PWD=$(pwd)"
echo "JOBID=${SLURM_JOB_ID:-NA}"

source /home/WUR/husen002/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/nobackup/WUR/ABGC/shared/pipelines_version2/envs/population-variant-calling-version2

# setelah conda aktif, kalau kamu mau strict boleh aktifkan -u:
set -euo pipefail

VCF=ALL_BREEDS_ARS1_AUTOSOME_Q20_DP5.norm.vcf.gz

echo "Samples:"
bcftools query -l ${VCF} | wc -l

# SNP-only
bcftools view -v snps -Oz -o ALL_BREEDS.SNP.vcf.gz ${VCF}
bcftools index -t ALL_BREEDS.SNP.vcf.gz

# set ID unik
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz \
  -o ALL_BREEDS.SNP.setID.vcf.gz ALL_BREEDS.SNP.vcf.gz
bcftools index -t ALL_BREEDS.SNP.setID.vcf.gz

# PLINK bed
plink --vcf ALL_BREEDS.SNP.setID.vcf.gz \
  --double-id --allow-extra-chr --make-bed --out ALL_BREEDS

# remove duplicate variants if any
plink --bfile ALL_BREEDS --allow-extra-chr \
  --list-duplicate-vars ids-only suppress-first --out ALL_BREEDS_dup

if [ -s ALL_BREEDS_dup.dupvar ]; then
  plink --bfile ALL_BREEDS --allow-extra-chr \
    --exclude ALL_BREEDS_dup.dupvar --make-bed --out ALL_BREEDS_nodup
else
  cp ALL_BREEDS.bed ALL_BREEDS_nodup.bed
  cp ALL_BREEDS.bim ALL_BREEDS_nodup.bim
  cp ALL_BREEDS.fam ALL_BREEDS_nodup.fam
fi

# LD prune
plink --bfile ALL_BREEDS_nodup --allow-extra-chr \
  --indep-pairwise 50 5 0.2 --out ALL_BREEDS_LD

# PCA
plink --bfile ALL_BREEDS_nodup --allow-extra-chr \
  --extract ALL_BREEDS_LD.prune.in --pca 20 --out ALL_BREEDS_PCA

ls -lh ALL_BREEDS_PCA.eigenvec ALL_BREEDS_PCA.eigenval
date
echo "DONE"

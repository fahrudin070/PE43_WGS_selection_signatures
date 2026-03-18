#!/bin/bash
#SBATCH --job-name=ROH_ADMIXREADY
#SBATCH --output=ROH_ADMIXREADY.fixed.out
#SBATCH --error=ROH_ADMIXREADY.fixed.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -eo pipefail   # <- HAPUS -u supaya conda deactivate tidak bikin FAILED

EXT=/lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/data/external_breeds

# Activate conda safely
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/nobackup/WUR/ABGC/shared/pipelines_version2/envs/population-variant-calling-version2

# Run ROH
plink --bfile ${EXT}/pca/admixture/merged.ADMIXREADY \
  --allow-extra-chr --chr-set 29 \
  --homozyg \
  --out ${EXT}/pca/roh/ROH_ADMIXREADY


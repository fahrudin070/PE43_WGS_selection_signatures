#!/bin/bash
#SBATCH --job-name=smk_map_NUB
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=snakemake_master_%j.out
#SBATCH --error=snakemake_master_%j.err

# PENTING: jangan pakai -u karena conda deactivate bisa memanggil variable yg belum ada
set -eo pipefail

source /home/WUR/husen002/miniconda3/etc/profile.d/conda.sh

# env snakemake
conda activate /lustre/nobackup/WUR/ABGC/shared/pipelines_version2/envs/population-variant-calling-version2

# tambah bwa-mem2 dari env PE43
export PATH=/lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/conda_envs/env_mapping_pe43/bin:$PATH

cd /lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/runs/mapping_nubian/mapping_pipeline

mkdir -p /lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/runs/mapping_nubian/output/logs_slurm

snakemake -j 8 --cluster-config cluster.yaml \
  --cluster "sbatch --mem={cluster.mem} --time={cluster.time} --cpus-per-task={cluster.threads} \
  --job-name={cluster.name} --output={cluster.output} --error={cluster.error}" \
  --configfile config.yaml

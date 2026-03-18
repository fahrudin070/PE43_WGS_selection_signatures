#!/bin/bash
#SBATCH --job-name=smk_map_BBG
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=snakemake_master_%j.out
#SBATCH --error=snakemake_master_%j.err

set -eo pipefail   # <-- HILANGKAN -u

source /home/WUR/husen002/miniconda3/etc/profile.d/conda.sh

conda activate /lustre/nobackup/WUR/ABGC/shared/pipelines_version2/envs/population-variant-calling-version2

# bwa-mem2 dari env PE43
export PATH=/lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/conda_envs/env_mapping_pe43/bin:$PATH

cd /lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/runs/mapping_black_bengal/mapping_pipeline
mkdir -p /lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/runs/mapping_black_bengal/output/logs_slurm

snakemake -j 8 --cluster-config cluster.yaml \
  --cluster "sbatch --export=ALL --mem={cluster.mem} --time {cluster.time} --cpus-per-task {cluster.threads} \
  --job-name={cluster.name} --output={cluster.output} --error={cluster.error}" \
  --conda-frontend conda \
  --configfile config.yaml

#!/usr/bin/env bash
# Script: 02_map_one_sample.sh
# Purpose: Map paired-end reads from one sample to the ARS1 reference genome and generate sorted/indexed BAM files.
# Input:
#   - reference.fa
#   - sample_ids.txt
#   - paired FASTQ files per sample
# Output:
#   - sorted BAM
#   - BAM index
#   - flagstat report
#   - alignment statistics report
# Main tools:
#   - bwa
#   - samtools
#   - SLURM
# Notes:
#   - One-sample mapping script used in the PE43 workflow
#SBATCH --job-name=PE_map
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=/lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/logs/mapping/%x_%A_%a.out
#SBATCH --error=/lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/logs/mapping/%x_%A_%a.err

set -euo pipefail

BASE="/lustre/nobackup/WUR/ABGC/husen002/pe43_popvar"
REF="${BASE}/data/ref/ARS1/reference.fa"
SAMPLE_LIST="${BASE}/config/sample_ids.txt"
OUTDIR="${BASE}/results/bams"

mkdir -p "$OUTDIR"

# ---- Modules ----
source /etc/profile >/dev/null 2>&1 || true
module purge >/dev/null 2>&1 || true
module load legacy
module load bwa/gcc/64/0.7.17

# ---- Conda hook (disable nounset temporarily to avoid CONDA_BACKUP_* unbound) ----
set +u
source /home/WUR/husen002/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/nobackup/WUR/ABGC/shared/pipelines_version2/envs/population-variant-calling-version2
set -u

# ---- Sample ----
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
if [[ -z "${SAMPLE}" ]]; then
  echo "[ERROR] Empty sample for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 3
fi

R1="${BASE}/data/fastq/${SAMPLE}/${SAMPLE}_clean_1.fq.gz"
R2="${BASE}/data/fastq/${SAMPLE}/${SAMPLE}_clean_2.fq.gz"
OUTBAM="${OUTDIR}/${SAMPLE}.sorted.bam"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
  echo "[ERROR] FASTQ not found for ${SAMPLE}" >&2
  echo "[ERROR] R1=${R1}" >&2
  echo "[ERROR] R2=${R2}" >&2
  exit 4
fi

RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}\tPU:${SAMPLE}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "[INFO] Sample=${SAMPLE}"
echo "[INFO] bwa=$(which bwa)"
echo "[INFO] samtools=$(which samtools)"
echo "[INFO] Start mapping..."

bwa mem -t "$THREADS" -R "$RG" "$REF" "$R1" "$R2" \
  | samtools sort -@ "$THREADS" -o "$OUTBAM" -

samtools index -@ "$THREADS" "$OUTBAM"
samtools flagstat "$OUTBAM" > "${OUTBAM}.flagstat.txt"
samtools stats "$OUTBAM" > "${OUTBAM}.stats.txt"

echo "[INFO] Done ${SAMPLE}"

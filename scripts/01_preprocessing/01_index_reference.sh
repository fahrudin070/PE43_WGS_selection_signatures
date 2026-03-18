#!/usr/bin/env bash
# Script: 01_index_reference.sh
# Purpose: Build reference genome index files required for mapping and variant calling.
# Input:
#   - ARS1 reference genome FASTA
# Output:
#   - BWA index files
#   - samtools FASTA index (.fai)
# Main tools:
#   - bwa
#   - samtools
# Notes:
#   - Run before mapping and variant calling
set -euo pipefail

REF="/lustre/nobackup/WUR/ABGC/husen002/pe43_popvar/data/ref/ARS1/reference.fa"

echo "[INFO] Reference: $REF"
test -f "$REF" || { echo "[ERROR] Reference not found: $REF" >&2; exit 1; }

# Build BWA index (using bwa from module legacy)
echo "[INFO] Running bwa index..."
bwa index "$REF"

# Build FASTA index for tools like freebayes-parallel
echo "[INFO] Running samtools faidx..."
samtools faidx "$REF"

echo "[INFO] Done indexing reference"

#!/usr/bin/env bash
set -euo pipefail

REG_DIR=${1:-04_fst/regions}
GENE_BED=${2:-04_fst/ref/ARS1.genes.bed}
OUT_DIR=${3:-04_fst/genes}

mkdir -p "$OUT_DIR"
pairs=(PE_vs_Nubian PE_vs_Bengal PE_vs_TG PE_vs_GBG PE_vs_HBG)

for p in "${pairs[@]}"; do
  bed="$REG_DIR/${p}.top1pct.WEIGHTED.bed"
  [[ -f "$bed" ]] || { echo "Missing: $bed" >&2; exit 1; }

  out="$OUT_DIR/${p}.WEIGHTED.genes.tsv"

  bedtools intersect -a "$bed" -b "$GENE_BED" -wa -wb \
    | awk 'BEGIN{FS=OFS="\t"}{print $8}' \
    | sort -u > "$out"

  echo "OK: $p -> $out"
done

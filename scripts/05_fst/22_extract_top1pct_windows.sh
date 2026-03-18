#!/usr/bin/env bash
set -euo pipefail

RAW_DIR=${1:-04_fst/raw}
BOTH_TSV=${2:-04_fst/summary/FST_summary_table.BOTH.tsv}
OUT_DIR=${3:-04_fst/regions}

mkdir -p "$OUT_DIR"
pairs=(PE_vs_Nubian PE_vs_Bengal PE_vs_TG PE_vs_GBG PE_vs_HBG)

# figure out p99_weighted column index
p99_col=$(head -n1 "$BOTH_TSV" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="P99_weighted") print i}')
[[ -n "$p99_col" ]] || { echo "Cannot find P99_weighted column in $BOTH_TSV" >&2; exit 1; }

for p in "${pairs[@]}"; do
  f="$RAW_DIR/${p}.50kb.windowed.weir.fst"
  [[ -f "$f" ]] || { echo "Missing: $f" >&2; exit 1; }

  thr=$(awk -F'\t' -v P="$p" -v c="$p99_col" 'NR>1 && $1==P {print $c}' "$BOTH_TSV")
  [[ -n "$thr" ]] || { echo "Missing P99_weighted for $p in $BOTH_TSV" >&2; exit 1; }

  out="$OUT_DIR/${p}.top1pct.WEIGHTED.bed"

  hdr=$(head -n 1 "$f" | tr -d '\r')

  chrom_col=$(echo "$hdr" | awk 'BEGIN{IGNORECASE=1}{n=split($0,a,/[\t ]+/); for(i=1;i<=n;i++) if(tolower(a[i])=="chrom"||tolower(a[i])=="chr"){print i; exit}}')
  start_col=$(echo "$hdr" | awk 'BEGIN{IGNORECASE=1}{n=split($0,a,/[\t ]+/); for(i=1;i<=n;i++) if(tolower(a[i])=="bin_start"||tolower(a[i])=="start"){print i; exit}}')
  end_col=$(echo "$hdr"   | awk 'BEGIN{IGNORECASE=1}{n=split($0,a,/[\t ]+/); for(i=1;i<=n;i++) if(tolower(a[i])=="bin_end"||tolower(a[i])=="end"){print i; exit}}')
  wcol=$(echo "$hdr"      | awk 'BEGIN{IGNORECASE=1}{n=split($0,a,/[\t ]+/); for(i=1;i<=n;i++) if(tolower(a[i])=="weighted_fst"||tolower(a[i])=="weightedfst"){print i; exit}}')

  [[ -n "$chrom_col" && -n "$start_col" && -n "$end_col" && -n "$wcol" ]] || {
    echo "Cannot detect required columns in $f" >&2
    echo "Header: $hdr" >&2
    exit 1
  }

  awk -v FS='[ \t]+' -v OFS="\t" \
      -v C="$chrom_col" -v S="$start_col" -v E="$end_col" -v W="$wcol" \
      -v THR="$thr" -v PAIR="$p" \
      'NR>1{
        gsub(/\r/,"",$W);
        if($W!="NA" && $W!="" && ($W+0==$W) && $W>=THR){
          print $C, $S-1, $E, PAIR
        }
      }' "$f" | sort -k1,1 -k2,2n > "$out"

  echo "OK: $p  P99_weighted=$thr  -> $out"
done

#!/usr/bin/env Rscript
# Script: 20_plot_roh_landscape.R
# Purpose: Generate ROH frequency landscape plot and identify genomic ROH enrichment patterns.
# Input:
#   - ROH window/frequency summary tables
# Output:
#   - final ROH landscape figure
# Main packages:
#   - ggplot2
#   - data.table
# Notes:
#   - Final ROH visualization used in manuscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# =========================
# CONFIG
# =========================
infile  <- "tables/PE43_ROHfreq_50kb.tsv"
out_png <- "plots/Fig_ROHfreq_50kb_thr30_FINAL.png"

N_TOTAL <- 43
thr_pct <- 30
thr_n   <- ceiling(thr_pct/100 * N_TOTAL)   # 13
TOP_LABEL_N <- 10

COL_PINK <- "#ff4fb3"
COL_BLUE <- "#2f7bff"

# =========================
# READ (7 columns, no header)
# V1 chr | V2 start | V3 end | V4 winID | V5 count | V6 N | V7 freq(0-1)
# =========================
dt <- fread(infile, header = FALSE)
if (ncol(dt) < 7) stop("Input harus 7 kolom: chr start end winID count N freq")

setnames(dt, 1:7, c("chr","start","end","win_id","count","Ncol","freq"))

dt[, `:=`(
  chr   = as.character(chr),
  start = as.numeric(start),
  end   = as.numeric(end),
  count = as.numeric(count),
  Ncol  = as.numeric(Ncol),
  freq  = as.numeric(freq)
)]

dt <- dt[!is.na(start) & !is.na(end) & end > start]
dt[is.na(count), count := 0]
dt[is.na(Ncol) | Ncol <= 0, Ncol := N_TOTAL]

# % individuals with ROH in window
dt[, pct_indiv := 100 * freq]
bad <- is.na(pct_indiv) | pct_indiv < 0 | pct_indiv > 100
if (any(bad)) dt[bad, pct_indiv := 100 * count / Ncol]
dt <- dt[pct_indiv >= 0 & pct_indiv <= 100]

# =========================
# AUTOSOMES ONLY: NC_030808.1 ... NC_030836.1  (Chr1..Chr29)
# =========================
dt <- dt[grepl("^NC_0308(0[8-9]|[1-2][0-9]|3[0-6])\\.", chr)]
if (nrow(dt) == 0) stop("Setelah filter autosome, dt kosong. Cek kolom chr.")

# Convert chr to Chr1..Chr29
dt[, chr_base := sub("\\..*$", "", chr)]                     # NC_030808
dt[, chr_num  := as.integer(sub("^NC_0*", "", chr_base))]    # 30808
dt[, chr_idx  := chr_num - 30807L]                           # 1..29
dt <- dt[chr_idx >= 1 & chr_idx <= 29]
dt[, chr_short := paste0("Chr", chr_idx)]

setorder(dt, chr_idx, start)

# =========================
# Manhattan-like genome coordinate
# =========================
chr_len <- dt[, .(chr_len = max(end, na.rm=TRUE)), by=.(chr, chr_idx, chr_short)]
setorder(chr_len, chr_idx)
chr_len[, chr_start := shift(cumsum(chr_len), fill=0)]
chr_len[, chr_mid   := chr_start + chr_len/2]

dt <- merge(dt, chr_len[, .(chr, chr_start, chr_len, chr_mid, chr_idx, chr_short)],
            by="chr", all.x=TRUE)
dt[, pos := start + chr_start]

# alternating colors by chromosome
dt[, chr_group := ifelse(chr_idx %% 2 == 1, "odd", "even")]

# =========================
# Threshold 30%
# =========================
dt[, is_thr30 := pct_indiv >= thr_pct]
cand <- dt[is_thr30 == TRUE]
setorder(cand, -pct_indiv)

topN <- cand[1:min(TOP_LABEL_N, .N)]
if (nrow(topN) > 0) {
  topN[, label := sprintf("%s:%d-%d (%.1f%%; %d/%d)",
                          chr_short, start, end, pct_indiv,
                          ceiling(pct_indiv/100 * N_TOTAL), N_TOTAL)]
}

subtitle_txt <- sprintf(
  "ROH frequency per 50 kb window (PE43). Threshold %d%% (≥%d/%d individuals). Top %d peaks labeled.",
  thr_pct, thr_n, N_TOTAL, min(TOP_LABEL_N, nrow(topN))
)

# =========================
# PLOT
# =========================
p <- ggplot(dt, aes(x=pos, y=pct_indiv)) +
  geom_point(aes(color=chr_group, alpha=is_thr30), size=0.9) +
  scale_alpha_manual(values=c("FALSE"=0.20, "TRUE"=0.95), guide="none") +
  scale_color_manual(values=c("odd"=COL_PINK, "even"=COL_BLUE), guide="none") +
  geom_hline(yintercept=thr_pct, linetype="dashed", linewidth=0.7, color="black") +
  annotate("text",
           x=max(dt$pos, na.rm=TRUE), y=thr_pct,
           label=sprintf("%d%% (≥%d/%d)", thr_pct, thr_n, N_TOTAL),
           hjust=1.02, vjust=-0.6, size=3.2) +
  scale_x_continuous(breaks=chr_len$chr_mid, labels=chr_len$chr_short,
                     expand=expansion(mult=c(0.002, 0.01))) +
  labs(
    title = "ROH frequency landscape (50 kb windows)",
    subtitle = subtitle_txt,
    x = "Chromosome",
    y = "% individuals with ROH in window"
  ) +
  theme_bw(base_size=12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=8),
    plot.title = element_text(face="bold", size=14)
  )

if (nrow(topN) > 0) {
  p <- p +
    geom_point(data=topN, aes(x=pos, y=pct_indiv),
               inherit.aes=FALSE, size=2.0, color="black") +
    geom_label_repel(
      data=topN,
      aes(x=pos, y=pct_indiv, label=label),
      inherit.aes=FALSE,
      size=2.8,
      label.size=0.2,
      min.segment.length=0,
      box.padding=0.30,
      point.padding=0.20,
      max.overlaps=Inf
    )
}

ggsave(out_png, plot=p, width=18, height=6.8, dpi=320)
message("DONE: ", out_png)

# Script: 25_plot_manhattan_fst.R
# Purpose: Generate weighted FST Manhattan-style plot across pairwise breed comparisons.
# Input:
#   - windowed weighted FST tables
# Output:
#   - final FST Manhattan figure
# Main packages:
#   - ggplot2
#   - data.table
# Notes:
#   - Final weighted FST plotting script used in manuscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript v7p3p2_LOCKED.R <raw_dir> <summary.BOTH.tsv> <karyo.tsv> <out.png> <weighted|mean>")
}
rawdir <- args[1]
sumtab <- args[2]
karyo  <- args[3]
outpng <- args[4]
prefer <- tolower(args[5])

if(!(prefer %in% c("weighted","mean"))) stop("prefer must be 'weighted' or 'mean'")

pairs <- c("PE_vs_Nubian","PE_vs_Bengal","PE_vs_TG","PE_vs_GBG","PE_vs_HBG")

pair_cols <- c(
  "PE_vs_Nubian" = "#0B1F5B",
  "PE_vs_Bengal" = "#123A7C",
  "PE_vs_TG"     = "#1C5D86",
  "PE_vs_GBG"    = "#1F7A6B",
  "PE_vs_HBG"    = "#2E8B57"
)

thr_cols <- c(
  "PE_vs_Nubian" = "#FFB000",
  "PE_vs_Bengal" = "#FF3B30",
  "PE_vs_TG"     = "#00C2FF",
  "PE_vs_GBG"    = "#B6FF00",
  "PE_vs_HBG"    = "#B66DFF"
)

pretty_pair <- function(x){
  x <- gsub("^PE_vs_", "PE vs ", x)
  x <- gsub("Bengal", "Black Bengal", x)
  x
}

map <- data.table(
  CHROM = sprintf("NC_030%03d.1", 808:836),
  chr_num = 1:29
)

ky_raw <- fread(karyo, header=FALSE)
setnames(ky_raw, c("CHROM","LEN"))
ky_raw[, CHROM := gsub("\r","", trimws(CHROM))]
ky_raw[, LEN := as.numeric(LEN)]

ky <- merge(map, ky_raw, by="CHROM", all.x=TRUE)
setorder(ky, chr_num)
if (any(is.na(ky$LEN))) stop("Missing LEN in karyo after merge")

ky[, offset := c(0, head(cumsum(LEN), -1))]
ky[, chr_end := offset + LEN]
genome_end <- max(ky$chr_end)
axis_df <- ky[, .(chr_num, center = offset + LEN/2, label = paste0("chr ", chr_num))]

# ---- thresholds from BOTH summary ----
st <- fread(sumtab)
setnames(st, tolower(names(st)))

thr_col <- if(prefer=="weighted") "p99_weighted" else "p99_mean"
if(!(thr_col %in% names(st))) stop(paste("Missing", thr_col, "in", sumtab))

thr_map <- setNames(st[[thr_col]], st$pair)
thr_dt <- data.table(Pair=pairs, thr=as.numeric(thr_map[pairs]))

# ---- read FST (LOCKED column) ----
read_one <- function(pair){
  fn <- file.path(rawdir, paste0(pair, ".50kb.windowed.weir.fst"))
  if (!file.exists(fn)) stop(paste("Missing:", fn))
  dt <- fread(fn)
  setnames(dt, tolower(names(dt)))

  chrom_col <- if ("chrom" %in% names(dt)) "chrom" else names(dt)[grepl("chrom", names(dt))][1]
  start_col <- if ("bin_start" %in% names(dt)) "bin_start" else names(dt)[grepl("start", names(dt))][1]
  end_col   <- if ("bin_end"   %in% names(dt)) "bin_end"   else names(dt)[grepl("end", names(dt))][1]

  fst_col <- if(prefer=="weighted") "weighted_fst" else "mean_fst"
  if(!(fst_col %in% names(dt))){
    stop(paste("ERROR:", fst_col, "not found in", fn, "\nCols:", paste(names(dt), collapse=", ")))
  }

  out <- data.table(
    CHROM = gsub("\r","", trimws(as.character(dt[[chrom_col]]))),
    BIN_START = as.integer(dt[[start_col]]),
    BIN_END   = as.integer(dt[[end_col]]),
    FST = suppressWarnings(as.numeric(dt[[fst_col]])),
    Pair = pair
  )

  out <- merge(out, map, by="CHROM", all.x=FALSE)
  out <- merge(out, ky[,.(chr_num, offset)], by="chr_num", all.x=FALSE)
  out[, pos := ((BIN_START + BIN_END)/2) + offset]
  out
}

all <- rbindlist(lapply(pairs, read_one), use.names=TRUE, fill=TRUE)
all_plot <- all[!is.na(FST) & is.finite(FST) & FST >= 0]

# fill NA thresholds (if any) from empirical quantile
for(i in 1:nrow(thr_dt)){
  if (is.na(thr_dt$thr[i])) {
    thr_dt$thr[i] <- quantile(all_plot[Pair==thr_dt$Pair[i], FST], 0.99, na.rm=TRUE)
  }
}

title_metric <- if(prefer=="weighted") "weighted" else "mean"

main <- ggplot(all_plot, aes(x=pos, y=FST, color=Pair)) +
  geom_point(size=0.55, alpha=0.26) +
  geom_hline(data=thr_dt, aes(yintercept=thr),
             linetype="dashed", linewidth=1.25, alpha=0.95,
             color=thr_cols[thr_dt$Pair]) +
  scale_color_manual(values=pair_cols,
                     labels=c("PE vs Nubian","PE vs Black Bengal","PE vs TG","PE vs GBG","PE vs HBG")) +
  scale_x_continuous(breaks=axis_df$center, labels=axis_df$label, expand=c(0.002,0)) +
  coord_cartesian(xlim=c(0, genome_end)) +
  labs(
    title=paste0("Genome-wide windowed ", title_metric, " FST (50 kb) — PE vs reference breeds (overlay)"),
    x="Chromosome",
    y=paste0("Windowed ", title_metric, " FST (50 kb)")
  ) +
  theme_bw(base_size=12) +
  theme(
    legend.position="bottom",
    legend.title=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    plot.title=element_text(face="bold"),
    axis.title=element_text(face="bold"),
    axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=9),
    plot.margin = margin(t=10, r=10, b=5, l=10)
  )

leg <- data.table(Pair=pairs, y=rev(seq_along(pairs)))
leg[, label := paste0("p99 threshold (dashed) — ", pretty_pair(Pair), " [", title_metric, "]")]
leg[, col := thr_cols[Pair]]

legend_panel <- ggplot(leg) +
  geom_segment(aes(x=0.02, xend=0.22, y=y, yend=y),
               linewidth=1.6, linetype="dashed", color=leg$col) +
  geom_text(aes(x=0.24, y=y, label=label),
            hjust=0, size=3.8) +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0.5, length(pairs)+0.5)) +
  theme_void(base_size=12) +
  theme(plot.margin = margin(t=0, r=10, b=10, l=10))

p <- main / legend_panel + plot_layout(heights=c(8.2, 1.4))
ggsave(outpng, p, width=22, height=11.7, dpi=600, bg="white")

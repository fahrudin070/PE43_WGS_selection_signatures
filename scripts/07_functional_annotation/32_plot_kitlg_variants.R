suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Paths
fst_bengal_file <- "tmp/KITLG_FST_Bengal_2Mb.tsv"
fst_nubian_file <- "tmp/KITLG_FST_Nubian_2Mb.tsv"
thr_bengal_file <- "tmp/THR_FST_Bengal_top1pct.txt"
thr_nubian_file <- "tmp/THR_FST_Nubian_top1pct.txt"
out_png <- "plots/PanelD_KITLG_FST_zoom_clean.png"

# Read data
bengal <- fread(fst_bengal_file)
nubian <- fread(fst_nubian_file)

# If no header, set colnames (vcftools windowed format)
if (!("BIN_START" %in% names(bengal))) {
  setnames(bengal, c("CHROM","BIN_START","BIN_END","N_VARIANTS","WEIGHTED_FST","MEAN_FST"))
}
if (!("BIN_START" %in% names(nubian))) {
  setnames(nubian, c("CHROM","BIN_START","BIN_END","N_VARIANTS","WEIGHTED_FST","MEAN_FST"))
}

# Ensure numeric
bengal[, WEIGHTED_FST := as.numeric(WEIGHTED_FST)]
nubian[, WEIGHTED_FST := as.numeric(WEIGHTED_FST)]

# Add group labels for legend
bengal[, Contrast := "PE vs Black Bengal"]
nubian[, Contrast := "PE vs Nubian"]

df <- rbindlist(list(bengal, nubian), use.names=TRUE, fill=TRUE)

# Read thresholds
thr_bengal <- suppressWarnings(as.numeric(gsub("[^0-9eE+\\.-]", "", readLines(thr_bengal_file)[1])))
thr_nubian <- suppressWarnings(as.numeric(gsub("[^0-9eE+\\.-]", "", readLines(thr_nubian_file)[1])))

if (is.na(thr_bengal) || is.na(thr_nubian)) {
  stop("Threshold file(s) could not be parsed. Check tmp/THR_FST_* files.")
}

# Threshold lines as a data.frame so they can appear in legend
thr_df <- data.table(
  y = c(thr_bengal, thr_nubian),
  ThrLabel = c("Top 1% cutoff (PE vs Bengal)", "Top 1% cutoff (PE vs Nubian)")
)

# Plot
p <- ggplot(df, aes(x=BIN_START/1e6, y=WEIGHTED_FST, color=Contrast)) +
  geom_point(size=2.2, alpha=0.75) +
  geom_hline(data=thr_df, aes(yintercept=y, linetype=ThrLabel),
             linewidth=0.9, show.legend=TRUE) +
  scale_color_manual(values=c("PE vs Black Bengal"="#d73027",
                              "PE vs Nubian"="#4575b4")) +
  scale_linetype_manual(values=c("Top 1% cutoff (PE vs Bengal)"="dashed",
                                 "Top 1% cutoff (PE vs Nubian)"="dashed")) +
  theme_classic(base_size=16) +
  labs(
    x="Position (Mb)",
    y="Weighted FST (50 kb windows)",
    title="KITLG region (2 Mb): PE vs Bengal vs Nubian",
    subtitle=paste0("Dashed lines = genome-wide top 1% FST cutoffs (Bengal=", signif(thr_bengal,3),
                    "; Nubian=", signif(thr_nubian,3), ")"),
    color="FST contrast",
    linetype="Genome-wide threshold"
  ) +
  theme(
    plot.title=element_text(face="bold"),
    legend.position="right"
  )

ggsave(out_png, p, width=10.5, height=6, dpi=400)
message("Saved: ", out_png)

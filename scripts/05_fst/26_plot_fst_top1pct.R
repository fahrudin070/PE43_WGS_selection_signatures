library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

plot_fst <- function(infile, thr_file, outpng, title){
  thr <- as.numeric(readLines(thr_file)[1])

  df <- read_tsv(infile, show_col_types=FALSE) %>%
    filter(!is.na(WEIGHTED_FST)) %>%
    mutate(CHR_NUM = as.integer(str_extract(CHROM, "(?<=NC_030)\\d+"))) %>%
    arrange(CHR_NUM, BIN_START) %>%
    mutate(CHR_GROUP = ifelse(CHR_NUM %% 2 == 0, "even", "odd"))

  p <- ggplot(df, aes(x=BIN_START/1e6, y=WEIGHTED_FST)) +
    geom_point(aes(color=CHR_GROUP), size=0.4, alpha=0.8) +
    facet_grid(. ~ CHROM, scales="free_x", space="free_x") +
    geom_hline(yintercept=thr, linetype="dashed") +
    theme_bw(base_size=10) +
    theme(
      legend.position="none",
      strip.text.x = element_text(size=6, angle=90),
      panel.spacing.x = unit(0.06, "lines")
    ) +
    labs(x="Position (Mb)", y="Windowed FST (50 kb)", title=title,
         subtitle=paste0("Dashed line = top 1% threshold (", signif(thr,4), ")"))

  ggsave(outpng, p, width=14, height=4.5, dpi=300)
}

plot_fst(
  "tables/FST_PE43_vs_Nubian17_50kb.windowed.weir.fst",
  "tmp/THR_FST_Nubian_top1pct.txt",
  "plots/Fig_STEP10A_FST_manhattan_PE43_vs_Nubian17_top1pct.png",
  "Genome-wide windowed FST (PE43 vs Nubian; 50 kb)"
)

plot_fst(
  "tables/FST_PE43_vs_Bengal6_50kb.windowed.weir.fst",
  "tmp/THR_FST_Bengal_top1pct.txt",
  "plots/Fig_STEP10B_FST_manhattan_PE43_vs_Bengal6_top1pct.png",
  "Genome-wide windowed FST (PE43 vs Black Bengal; 50 kb)"
)

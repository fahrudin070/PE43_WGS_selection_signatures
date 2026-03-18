# Script: 35_circos_plot_main.R
# Purpose: Generate integrative Circos plot combining ROH, FST, overlap regions, and candidate/core genes.
# Input:
#   - chromosome size table
#   - ROH summary tables
#   - FST top regions
#   - overlap/candidate gene tables
# Output:
#   - final Circos figure for manuscript
# Main packages:
#   - circlize
#   - data.table
#   - grid
# Notes:
#   - Main integrative visualization script used in manuscript
suppressPackageStartupMessages({
  library(circlize)
  library(data.table)
})

dir.create("plots", showWarnings = FALSE)

# ----- load sizes ordered Chr1..Chr29
chr <- fread("10_circos/chr_sizes_autosome.chr1to29.tsv", header=FALSE)
setnames(chr, c("chr","len"))

# label: Chr1..Chr29
chr[, label := {
  num <- chr
  num <- gsub("^NC_0308","", num)
  num <- gsub("\\.1$","", num)
  paste0("Chr", (as.integer(num) - 7))
}]

# ----- regions
fst <- fread("10_circos/PE_FST_top1pct_regions.supermerged.bed", header=FALSE)
setnames(fst, c("chr","start","end"))

roh <- fread("10_circos/PE_ROH_thr30.bed", header=FALSE)
setnames(roh, c("chr","start","end"))

ovl <- fread("10_circos/PE_FSTxROH_overlap_regions.bed", header=FALSE)
setnames(ovl, c("chr","start","end"))

# ----- genes
core <- fread("10_circos/core6_genes.bed", header=FALSE)
setnames(core, c("chr","start","end","gene"))
core[, mid := (start+end)/2]

cand50 <- fread("10_circos/cand50_genes.bed", header=FALSE)
setnames(cand50, c("chr","start","end","gene"))
cand50[, mid := (start+end)/2]

cons10 <- fread("10_circos/consensus5of5_10genes.bed", header=FALSE)
setnames(cons10, c("chr","start","end","gene"))
cons10[, mid := (start+end)/2]

# keep only chromosomes in chr sizes (safety)
keep_chr <- chr$chr
fst  <- fst [chr %in% keep_chr]
roh  <- roh [chr %in% keep_chr]
ovl  <- ovl [chr %in% keep_chr]
core <- core[chr %in% keep_chr]
cand50 <- cand50[chr %in% keep_chr]
cons10 <- cons10[chr %in% keep_chr]

# ----- brighter colors (high contrast)
COL_BG   <- "#08101F"   # dark navy
COL_FST  <- "#FF1FBF"   # neon magenta
COL_ROH  <- "#00D9FF"   # neon cyan
COL_OVL  <- "#FFE600"   # neon yellow
COL_GENE <- "#B8FF6A"   # neon green (candidate genes)
COL_TXT  <- "#F2F6FF"   # near-white

draw_circos <- function(){
  circos.clear()
  par(bg=COL_BG)

  # start at Chr1 on top
  circos.par(
    start.degree = 90,
    gap.degree = 2,
    track.margin = c(0.004, 0.004),
    points.overflow.warning = FALSE,
    cell.padding = c(0, 0, 0, 0)
  )

  circos.initialize(factors = chr$chr,
                    xlim = cbind(rep(0, nrow(chr)), chr$len))

  # Track 1: chromosome labels
  circos.trackPlotRegion(ylim=c(0,1), bg.border=NA, track.height=0.06,
    panel.fun=function(x,y){
      i <- which(chr$chr == CELL_META$sector.index)
      circos.text(CELL_META$xcenter, 0.5, chr$label[i],
                  col=COL_TXT, cex=0.60,
                  facing="clockwise", niceFacing=TRUE)
    })

  # Track 2: FST blocks
  circos.trackPlotRegion(ylim=c(0,1), bg.border=NA, track.height=0.10,
    panel.fun=function(x,y){
      this <- fst[chr == CELL_META$sector.index]
      if(nrow(this)>0) circos.rect(this$start, 0, this$end, 1, col=COL_FST, border=NA)
    })

  # Track 3: ROH blocks
  circos.trackPlotRegion(ylim=c(0,1), bg.border=NA, track.height=0.10,
    panel.fun=function(x,y){
      this <- roh[chr == CELL_META$sector.index]
      if(nrow(this)>0) circos.rect(this$start, 0, this$end, 1, col=COL_ROH, border=NA)
    })

  # Track 4: Overlap highlight blocks (FST×ROH)
  circos.trackPlotRegion(ylim=c(0,1), bg.border=NA, track.height=0.08,
    panel.fun=function(x,y){
      this <- ovl[chr == CELL_META$sector.index]
      if(nrow(this)>0) circos.rect(this$start, 0, this$end, 1, col=COL_OVL, border=NA)
    })

  # Track 5: all candidate genes (cand50) + consensus10 as marks (NOT labels, supaya tidak chaos)
  circos.trackPlotRegion(ylim=c(0,1), bg.border=NA, track.height=0.10,
    panel.fun=function(x,y){
      c50 <- cand50[chr == CELL_META$sector.index]
      if(nrow(c50)>0) circos.points(c50$mid, rep(0.45, nrow(c50)),
                                    pch=16, cex=0.35, col=COL_GENE)

      c10 <- cons10[chr == CELL_META$sector.index]
      if(nrow(c10)>0) circos.points(c10$mid, rep(0.70, nrow(c10)),
                                    pch=16, cex=0.45, col=COL_TXT)
    })

  # Track 6: core6 labels (only these get text labels)
  circos.trackPlotRegion(ylim=c(0,1), bg.border=NA, track.height=0.14,
    panel.fun=function(x,y){
      this <- core[chr == CELL_META$sector.index]
      if(nrow(this)>0){
        circos.points(this$mid, rep(0.25, nrow(this)), pch=16, cex=0.80, col=COL_OVL)
        circos.text(this$mid, rep(0.78, nrow(this)), this$gene,
                    col=COL_TXT, cex=0.60,
                    facing="clockwise", niceFacing=TRUE)
      }
    })
}

png("plots/Circos_PE_integrative_FINAL.png", width=2600, height=2600, res=300, bg=COL_BG)
draw_circos()
dev.off()

pdf("plots/Circos_PE_integrative_FINAL.pdf", width=9, height=9, bg=COL_BG)
draw_circos()
dev.off()

message("DONE: plots/Circos_PE_integrative_FINAL.png/.pdf")

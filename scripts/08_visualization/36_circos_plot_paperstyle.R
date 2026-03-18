suppressPackageStartupMessages({
  library(circlize)
  library(data.table)
  library(grDevices)
})

dir.create("plots", showWarnings = FALSE)

message("[1] Loading inputs...")

# ---------- chromosome sizes (29 autosomes) ----------
chr_sizes <- fread("10_circos/chr_sizes_autosome.tsv", header = FALSE)
setnames(chr_sizes, c("chr","len"))
stopifnot(nrow(chr_sizes) == 29)

# order: NC_030808.1..NC_030836.1  => Chr1..Chr29
chr_sizes[, chr_num := as.integer(sub("\\.1$", "", sub("^NC_0308", "", chr)))]
setorder(chr_sizes, chr_num)
chr_sizes[, label := paste0("Chr-", chr_num - 7)]   # 8->Chr-1 ... 36->Chr-29

valid_chr <- chr_sizes$chr
chr_label <- setNames(chr_sizes$label, chr_sizes$chr)

# ---------- ROH islands ----------
roh <- fread("10_circos/PE_ROH_thr30.bed", header = FALSE)
setnames(roh, c("chr","start","end"))
roh <- roh[chr %in% valid_chr]

# ---------- FST regions (per comparison) ----------
fst_bengal <- fread("10_circos/FST_PE_vs_Bengal.supermerged.bed", header = FALSE)
setnames(fst_bengal, c("chr","start","end"))
fst_bengal <- fst_bengal[chr %in% valid_chr]

fst_nubian <- fread("10_circos/FST_PE_vs_Nubian.supermerged.bed", header = FALSE)
setnames(fst_nubian, c("chr","start","end"))
fst_nubian <- fst_nubian[chr %in% valid_chr]

# ---------- all genes ticks (midpoint tick) ----------
genes <- fread("10_circos/all_genes_ticks.bed", header = FALSE)
setnames(genes, c("chr","start","end","gene"))
genes <- genes[chr %in% valid_chr]
genes[, mid := start]  # already midpoint (your file is mid..mid+1)

# ---------- core6 genes ----------
core <- fread("10_circos/core6_genes.bed", header = FALSE)
setnames(core, c("chr","start","end","gene"))
core <- core[chr %in% valid_chr]
core[, mid := (start + end)/2]

# Priority gene name (KITLG as purple in example)
priority_gene <- "KITLG"

message("[2] Drawing circos...")

# ---------- style (paper-like) ----------
COL_ROH      <- "black"
COL_BENGAL   <- "red"
COL_NUBIAN   <- "blue"
COL_CAND     <- "#00B894"   # green-ish
COL_PRIORITY <- "#8E44AD"   # purple
COL_TICK     <- adjustcolor("grey35", alpha.f = 0.40)

# output
png("plots/Circos_PE_paperstyle.png", width = 2400, height = 2400, res = 300)
par(bg = "white", mar = c(1,1,2,1), xpd = NA)

circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = 2,
  track.margin = c(0.004, 0.004),
  points.overflow.warning = FALSE,
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(
  factors = chr_sizes$chr,
  xlim = cbind(rep(0, nrow(chr_sizes)), chr_sizes$len)
)

# ===== Track 1: chromosome labels (outer) =====
circos.trackPlotRegion(
  ylim = c(0,1),
  bg.border = NA,
  track.height = 0.07,
  panel.fun = function(x, y){
    ch  <- CELL_META$sector.index
    lab <- chr_label[[ch]]
    circos.text(
      CELL_META$xcenter, 0.5, lab,
      facing = "clockwise", niceFacing = TRUE,
      cex = 0.65, font = 2
    )
  }
)

# ===== Track 2: ROH islands (black blocks) =====
circos.trackPlotRegion(
  ylim = c(0,1),
  bg.border = NA,
  track.height = 0.09,
  panel.fun = function(x,y){
    ch <- CELL_META$sector.index
    d  <- roh[chr == ch]
    if(nrow(d) > 0){
      circos.rect(d$start, 0, d$end, 1, col = COL_ROH, border = NA)
    }
  }
)

# ===== Track 3: FST outlier windows PE vs Bengal (red) =====
circos.trackPlotRegion(
  ylim = c(0,1),
  bg.border = NA,
  track.height = 0.08,
  panel.fun = function(x,y){
    ch <- CELL_META$sector.index
    d  <- fst_bengal[chr == ch]
    if(nrow(d) > 0){
      circos.rect(d$start, 0, d$end, 1, col = adjustcolor(COL_BENGAL, 0.55), border = NA)
    }
  }
)

# ===== Track 4: FST outlier windows PE vs Nubian (blue) =====
circos.trackPlotRegion(
  ylim = c(0,1),
  bg.border = NA,
  track.height = 0.08,
  panel.fun = function(x,y){
    ch <- CELL_META$sector.index
    d  <- fst_nubian[chr == ch]
    if(nrow(d) > 0){
      circos.rect(d$start, 0, d$end, 1, col = adjustcolor(COL_NUBIAN, 0.55), border = NA)
    }
  }
)

# ===== Track 5: ALL genes as ticks + candidate dots =====
# NOTE: This is “all genes on every chromosome”.
# If it becomes too dense, later we can downsample, but this matches your request.
circos.trackPlotRegion(
  ylim = c(0,1),
  bg.border = NA,
  track.height = 0.11,
  panel.fun = function(x,y){
    ch <- CELL_META$sector.index
    d  <- genes[chr == ch]
    if(nrow(d) > 0){
      # tiny tick marks (vertical)
      circos.segments(d$mid, 0.05, d$mid, 0.35, col = COL_TICK, lwd = 0.4)
    }

    # show core genes as green dots (candidate genes)
    cc <- core[chr == ch]
    if(nrow(cc) > 0){
      circos.points(cc$mid, rep(0.60, nrow(cc)), pch = 16, cex = 0.9, col = COL_CAND)
    }
  }
)

# ===== Track 6: core gene labels (KITLG as priority purple) =====
circos.trackPlotRegion(
  ylim = c(0,1),
  bg.border = NA,
  track.height = 0.16,
  panel.fun = function(x,y){
    ch <- CELL_META$sector.index
    cc <- core[chr == ch]
    if(nrow(cc) > 0){
      for(i in seq_len(nrow(cc))){
        g  <- cc$gene[i]
        mx <- cc$mid[i]

        col_txt <- ifelse(g == priority_gene, COL_PRIORITY, "black")
        cex_txt <- ifelse(g == priority_gene, 0.85, 0.65)
        font_tx <- ifelse(g == priority_gene, 2, 1)

        circos.text(mx, 0.5, paste0(g, " (", chr_label[[ch]], ")"),
                    facing = "inside", niceFacing = TRUE,
                    col = col_txt, cex = cex_txt, font = font_tx)
        # connector
        circos.segments(mx, 0.02, mx, 0.40, col = "grey30", lwd = 0.6)
      }
    }
  }
)

# ----- Legend (paper style, bottom) -----
par(xpd = NA)
legend("bottom", inset = -0.03, horiz = TRUE, bty = "n", cex = 0.85,
       legend = c("ROH islands (PE)", "FST outlier windows: PE vs Bengal",
                  "FST outlier windows: PE vs Nubian", "Candidate genes (core overlap)",
                  "Priority candidate gene: KITLG"),
       fill   = c(COL_ROH, COL_BENGAL, COL_NUBIAN, COL_CAND, COL_PRIORITY))

dev.off()
circos.clear()

message("DONE: plots/Circos_PE_paperstyle.png")

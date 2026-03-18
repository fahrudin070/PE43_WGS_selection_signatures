args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 4){
  stop("Usage: Rscript plot_roh_occurrence_pretty.R <POP> <NPOP> <FRAC> <GAPTAG>\nExample: Rscript ... PE43 43 0.2 20pct")
}
pop  <- args[1]
npop <- as.numeric(args[2])
frac <- as.numeric(args[3])          # 0.2 or 0.3
tag  <- args[4]                      # "20pct" or "30pct"

# input files
occ_file <- paste0("roh_occ_", pop, ".tsv")
isl_file <- paste0("ROH_islands_", pop, "_", tag, ".bed")

df <- read.table(occ_file, header=TRUE, sep="\t")
thr <- frac * npop

# ===== Fig2: Boxplot occurrence per SNP =====
png(paste0("Fig2_Boxplot_ROHocc_", pop, "_", tag, ".png"), width=2600, height=1500, res=250)
par(mar=c(6,5,4,2))
boxplot(df$ROH_OCC,
        main=paste0("ROH occurrence per SNP (", pop, ")  |  threshold=", tag, " (", thr, " individuals)"),
        ylab="Number of individuals with ROH covering SNP",
        xlab="",
        col="grey90",
        border="grey40",
        outline=FALSE)
stripchart(df$ROH_OCC, method="jitter", pch=16, cex=0.25, col=rgb(0,0,0,0.15), add=TRUE)
abline(h=thr, lwd=2, col="red")
legend("topright",
       legend=c("Threshold line"),
       lwd=2, col="red", bty="n")
dev.off()

# ===== Prepare Manhattan coordinates =====
chr <- df$CHR
pos <- df$POS
chrs <- sort(unique(chr))
offset <- 0
x <- numeric(nrow(df))
ticks <- numeric(length(chrs))
cols <- rep(c("grey30","grey65"), length.out=length(chrs))

for(i in seq_along(chrs)){
  c <- chrs[i]
  idx <- which(chr==c)
  x[idx] <- pos[idx] + offset
  ticks[i] <- mean(range(x[idx]))
  offset <- max(x[idx])
}

# ===== read islands (if any) =====
has_islands <- file.exists(isl_file) && (file.info(isl_file)$size > 20)
islands <- NULL
if(has_islands){
  islands <- read.table(isl_file, header=TRUE, sep="\t")
  # convert islands to x coords using same offsets
  chr_offsets <- setNames(rep(0, length(chrs)), chrs)
  off <- 0
  for(i in seq_along(chrs)){
    c <- chrs[i]
    idx <- which(chr==c)
    chr_offsets[as.character(c)] <- off
    off <- max(x[idx])
  }
  islands$XSTART <- islands$START + chr_offsets[as.character(islands$CHR)]
  islands$XEND   <- islands$END   + chr_offsets[as.character(islands$CHR)]
}

# ===== Fig3: Manhattan ROH occurrence + islands highlight =====
png(paste0("Fig3_Manhattan_ROHocc_", pop, "_", tag, ".png"), width=4200, height=1600, res=250)
par(mar=c(6,5,4,2))

# plot points by chr with alternating colors
plot(NA, xlim=c(min(x), max(x)), ylim=c(0, max(df$ROH_OCC, thr)*1.05),
     xaxt="n",
     xlab="Chromosome",
     ylab="ROH occurrence (individual count)",
     main=paste0("ROH occurrence per SNP (", pop, ") | islands=", tag))

for(i in seq_along(chrs)){
  c <- chrs[i]
  idx <- which(chr==c)
  points(x[idx], df$ROH_OCC[idx], pch=16, cex=0.20, col=cols[i])
}

axis(1, at=ticks, labels=chrs, cex.axis=0.9)
abline(h=thr, col="red", lwd=2)

# highlight islands
if(!is.null(islands) && nrow(islands) > 0){
  for(i in 1:nrow(islands)){
    rect(islands$XSTART[i], -1, islands$XEND[i], max(df$ROH_OCC, thr)*1.05,
         col=rgb(1,0,0,0.12), border=NA)
  }
  legend("topright",
         legend=c("Threshold", "ROH islands"),
         lwd=c(2, NA),
         pch=c(NA, 15),
         col=c("red", rgb(1,0,0,0.35)),
         pt.cex=2,
         bty="n")
} else {
  legend("topright", legend=c("Threshold"), lwd=2, col="red", bty="n")
}
dev.off()

cat("Saved pretty plots for", pop, tag, "\n")

args <- commandArgs(trailingOnly=TRUE)
K <- as.integer(args[1])

qfile <- sprintf("ALL_BREEDS_LDpruned.numeric.%d.Q", K)
fam <- "../plink/ALL_BREEDS_LDpruned.numeric.fam"
grp <- "../meta/sample_groups.tsv"

Q <- as.matrix(read.table(qfile))
famdat <- read.table(fam, stringsAsFactors=FALSE)
iid <- famdat$V2

g <- read.table(grp, sep="\t", stringsAsFactors=FALSE)
colnames(g) <- c("IID","BREED")

m <- match(iid, g$IID)
breed <- g$BREED[m]
if(any(is.na(breed))) stop("Ada IID di FAM yang tidak ketemu di sample_groups.tsv")

ord <- order(breed, iid)
Q <- Q[ord,]
breed <- breed[ord]

# warna default R yang cukup jelas
cols <- rainbow(K)

png(sprintf("ADMIXTURE_K%d_barplot.png",K), width=2600, height=900, res=200)
par(mar=c(9,4,2,1))

bp <- barplot(t(Q), col=cols, border=NA, space=0, axes=FALSE)
axis(2)

# garis batas antar breed
r <- rle(breed)
cuts <- cumsum(r$lengths)
abline(v=cuts + 0.5, col="grey70", lwd=1)

# label breed di bawah
mid <- cuts - r$lengths/2
text(x=mid, y=-0.03, labels=r$values, srt=45, adj=1, xpd=NA, cex=0.9)

dev.off()


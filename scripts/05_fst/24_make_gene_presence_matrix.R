suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2){
  stop("Usage: Rscript make_gene_presence_matrix_weighted.R <genes_dir> <out_prefix>")
}
genes_dir <- args[1]
out_prefix <- args[2]

pairs <- c("PE_vs_Nubian","PE_vs_Bengal","PE_vs_TG","PE_vs_GBG","PE_vs_HBG")

lst <- lapply(pairs, function(p){
  f <- file.path(genes_dir, paste0(p, ".WEIGHTED.genes.tsv"))
  if(!file.exists(f)) stop(paste("Missing:", f))
  dt <- fread(f, header=FALSE)
  setnames(dt, "gene")
  dt[, pair := p]
  unique(dt)
})

all <- rbindlist(lst)
mat <- dcast(all, gene ~ pair, fun.aggregate=length, value.var="pair")
for(p in pairs) mat[[p]] <- as.integer(mat[[p]] > 0)
mat[, presence_count := rowSums(.SD), .SDcols=pairs]
setorder(mat, -presence_count, gene)

f_matrix <- paste0(out_prefix, ".WEIGHTED.gene_presence_matrix.tsv")
fwrite(mat, f_matrix, sep="\t")

fwrite(mat[presence_count==5, .(gene)],
       paste0(out_prefix, ".WEIGHTED.core5of5.genes.txt"),
       col.names=FALSE)

fwrite(mat[presence_count>=4, .(gene)],
       paste0(out_prefix, ".WEIGHTED.core4of5.genes.txt"),
       col.names=FALSE)

cat("Wrote:\n", f_matrix, "\n")

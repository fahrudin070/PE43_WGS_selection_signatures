suppressPackageStartupMessages({
  library(ggvenn)
  library(readr)
})

# ---- input files ----
f_roh    <- "tables/VENN_A_ROH_islands.txt"
f_nubian <- "tables/VENN_B_ROH_FST_Nubian.txt"
f_bengal <- "tables/VENN_C_ROH_FST_Bengal.txt"

# ---- read as unique gene sets ----
read_set <- function(f){
  x <- read_lines(f)
  x <- trimws(x)
  x <- x[x != ""]
  unique(x)
}

sets <- list(
  "ROH islands (thr30)"         = read_set(f_roh),
  "ROH ∩ FST (PE vs Nubian)"    = read_set(f_nubian),
  "ROH ∩ FST (PE vs Bengal)"    = read_set(f_bengal)
)

# ---- counts (print to console) ----
cat("Set sizes:\n")
print(sapply(sets, length))

core <- Reduce(intersect, sets)
cat("\nCore intersection (3-way):", length(core), "\n")

# ---- plot ----
p <- ggvenn(
  sets,
  show_percentage = FALSE,
  show_elements = FALSE,
  stroke_size = 0.8,
  set_name_size = 4
)

ggsave("figures/Venn_ROH_FST_Nubian_Bengal.png", plot = p, width = 9, height = 6, dpi = 300)
ggsave("figures/Venn_ROH_FST_Nubian_Bengal.pdf", plot = p, width = 9, height = 6)

# ---- export core list ----
writeLines(core, "tables/CORE_genes_ROH_FST_Nubian_Bengal.txt")

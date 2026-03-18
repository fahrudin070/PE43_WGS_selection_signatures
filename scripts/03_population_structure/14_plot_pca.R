# Script: 14_plot_pca.R
# Purpose: Generate PCA plot for PE goats and reference breeds.
# Input:
#   - PCA eigenvec/eigenval files
#   - sample/population metadata
# Output:
#   - final PCA figure
# Main packages:
#   - ggplot2
#   - data.table
# Notes:
#   - Final PCA plotting version used in manuscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(maps)
  library(mapdata)
  library(cowplot)
  library(ggrepel)
  library(viridis)
})

# -------------------------
# Input PCA
# -------------------------
pca <- read.table("ALL_BREEDS_PCA.eigenvec", header = FALSE)
colnames(pca)[1:6] <- c("FID","IID","PC1","PC2","PC3","PC4")

eig <- scan("ALL_BREEDS_PCA.eigenval")
pc1p <- round(100*eig[1]/sum(eig), 2)
pc2p <- round(100*eig[2]/sum(eig), 2)

# -------------------------
# SRR mapping (Nubian / Black Bengal)
# -------------------------
map_srr <- readr::read_tsv("srr_to_breed.tsv", show_col_types = FALSE) %>%
  mutate(
    Breed = case_when(
      Breed %in% c("BlackBengal","Black Bengal") ~ "Black Bengal",
      Breed %in% c("Nubian") ~ "Nubian",
      TRUE ~ Breed
    )
  )

# -------------------------
# PCA breed labels (FULL NAMES for legend)
# -------------------------
pca2 <- pca %>%
  left_join(map_srr, by = c("IID"="IID")) %>%
  mutate(
    BreedFull = case_when(
      str_detect(IID, "^PEB") ~ "Peranakan Etawa",
      str_detect(IID, "^GBG") ~ "Guizhou black goats",
      str_detect(IID, "^HBG") ~ "Hezhang black goats",
      str_detect(IID, "^TG")  ~ "Tashi goats",
      str_detect(IID, "^SRR") & !is.na(Breed) ~ Breed,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(BreedFull))

# Legend label (country/context)
legend_label <- function(b) {
  dplyr::case_when(
    b == "Peranakan Etawa"       ~ "Indonesia:\nPeranakan Etawa (PE)",
    b == "Nubian" ~ "India (proxy ancestry):\nNubian (NB)",
    b == "Black Bengal"          ~ "India:\nBlack Bengal (BB)",
    b == "Guizhou black goats"   ~ "China:\nGuizhou black goats (GBG)",
    b == "Hezhang black goats"   ~ "China:\nHezhang black goats (HBG)",
    b == "Tashi goats"           ~ "China:\nTashi goats (TG)",
    TRUE ~ b
  )
}

breed_order <- c(
  "Peranakan Etawa",
  "Nubian",
  "Black Bengal",
  "Guizhou black goats",
  "Hezhang black goats",
  "Tashi goats"
)

pca2$BreedFull <- factor(pca2$BreedFull, levels = breed_order)
pca2 <- pca2 %>%
  mutate(BreedLegend = factor(map_chr(as.character(BreedFull), legend_label),
                              levels = map_chr(breed_order, legend_label)))

# palette
pal <- viridis(length(levels(pca2$BreedLegend)), option = "D", end = 0.95)
names(pal) <- levels(pca2$BreedLegend)

# -------------------------
# PCA panel
# -------------------------
p_pca <- ggplot(pca2, aes(PC1, PC2, color = BreedLegend)) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.4, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.4, alpha = 0.6) +
  geom_point(alpha = 0.9, size = 1.8) +
  stat_ellipse(aes(group = BreedLegend),
               level = 0.95, type = "norm",
               linewidth = 0.65, show.legend = FALSE, na.rm = TRUE) +
  scale_color_manual(values = pal) +
  guides(color = guide_legend(title = "Goat Breed")) +
  labs(
    x = paste0("PC1 (", pc1p, "%)"),
    y = paste0("PC2 (", pc2p, "%)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face="bold"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold")
  )

# -------------------------
# MAP panel (FINAL)
#   - Hanya 6 breed (PE, NB, BB, GBG, HBG, TG)
#   - Tidak ada ancestry silang
# -------------------------
world <- map_data("world")

geo_all <- readr::read_tsv("breed_geo_final.tsv", show_col_types = FALSE)

geo_breed <- geo_all %>%
  mutate(
    # samakan nama supaya match dengan legend/palette
    BreedFull = case_when(
      BreedFull %in% c("Guizhou Black") ~ "Guizhou black goats",
      BreedFull %in% c("Hezhang Black") ~ "Hezhang black goats",
      BreedFull %in% c("Tashi Goat")    ~ "Tashi goats",
      TRUE ~ BreedFull
    ),
    BreedFull = factor(BreedFull, levels = breed_order),
    BreedLegend = factor(map_chr(as.character(BreedFull), legend_label),
                         levels = levels(pca2$BreedLegend))
  )

# zoom bbox
pad_x <- 18
pad_y <- 12
xlim <- range(geo_breed$Longitude, na.rm = TRUE) + c(-pad_x, pad_x)
ylim <- range(geo_breed$Latitude,  na.rm = TRUE) + c(-pad_y, pad_y)

p_map <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey92", color = "white", linewidth = 0.2) +

  # breed points (colored)
  geom_point(data = geo_breed, aes(Longitude, Latitude, color = BreedLegend),
             size = 3.6, alpha = 0.95) +
  ggrepel::geom_text_repel(
    data = geo_breed, aes(Longitude, Latitude, label = Abbrev),
    size = 4.2, fontface = "bold",
    box.padding = 0.55, point.padding = 0.35,
    segment.size = 0.25, segment.alpha = 0.7,
    show.legend = FALSE
  ) +

  coord_quickmap(xlim = xlim, ylim = ylim) +
  scale_color_manual(values = pal) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold")
  ) +
  labs(x = "Longitude", y = "Latitude")

final <- cowplot::plot_grid(p_pca, p_map, nrow = 1, rel_widths = c(1.2, 1))

ggsave("PCA_modern_with_map_v3.png", final, width = 15, height = 6.4, dpi = 300)
ggsave("PCA_modern_with_map_v3.pdf", final, width = 15, height = 6.4)

#!/usr/bin/env Rscript

# Script: h1h2_dotplot.R
# Purpose: Parse minimap2 PAF file and generate a synteny dotplot colored by identity percentage.

suppressPackageStartupMessages({
library(tidyverse)
})

# Set up file paths and load top 11 largest contigs for hap1, hap2
paf_file    <- "03_assembly/results/h1h2/hap1_vs_hap2.paf"
h1_ids_file <- "03_assembly/results/hifiasm/post_hifi_stats/top11_hap1_ids.txt"
h2_ids_file <- "03_assembly/results/hifiasm/post_hifi_stats/top11_hap2_ids.txt"
out_dir     <- "03_assembly/results/h1h2"

top11_h1 <- read_lines(h1_ids_file)
top11_h2 <- read_lines(h2_ids_file)

# Load .paf file from alignment and filter to top 11 contigs for plotting
# also restrict to large enough alignment blocks, large enough aggregate query alignment
min_align_len <- 2000     # -m filter
min_query_len <- 500000   # -q filter

paf_cols <- c("qname", "qlen", "qstart", "qend", "strand",
              "tname", "tlen", "tstart", "tend", 
              "nmatch", "alen", "mapq")

paf <- read_tsv(
  paf_file, 
  col_names = paf_cols, 
  col_select = 1:12, 
  show_col_types = FALSE
)

# Compute alignment identity %
paf_top11 <- paf %>%
  filter(qname %in% top11_h1 & tname %in% top11_h2) %>%
  filter(alen >= min_align_len) %>%                  # Apply -m filter
  group_by(qname) %>%
  filter(sum(alen) >= min_query_len) %>%              # Apply -q filter
  ungroup() %>%
  mutate(
    pct_identity = (nmatch / alen) * 100,
    qstart_plot = if_else(strand == "-", qend, qstart),
    qend_plot   = if_else(strand == "-", qstart, qend)
  )

# Generate Dotplot
# Sort contigs by length/ID order in factors for clean plot axes
paf_top11 <- paf_top11 %>%
  mutate(
    qname = factor(qname, levels = rev(top11_h1)), # Hap1 on Y-axis
    tname = factor(tname, levels = top11_h2)       # Hap2 on X-axis
  )

plot_title <- paste0(
  "Hap1 vs Hap2 Alignment Dotplot (Top 11 Contigs)\n",
  "Alignments: ", nrow(paf_top11)
)

dotplot <- ggplot(paf_top11) +
  geom_segment(
    aes(
      x = tstart / 1e6, 
      xend = tend / 1e6, 
      y = qstart_plot / 1e6, 
      yend = qend_plot / 1e6, 
      color = pct_identity
    ),
    linewidth = 0.8
  ) +
  facet_grid(qname ~ tname, scales = "free", space = "free") +
  scale_color_viridis_c(
    option = "plasma",
    name = "Percent Identity (%)"
  ) +
  theme_bw() +
  labs(
    title = plot_title,
    x = "Hap2 Reference Contigs (Mb)",
    y = "Hap1 Query Contigs (Mb)"
  ) +
  theme(
    strip.text.x = element_text(size = 7, angle = 90),
    strip.text.y = element_text(size = 7, angle = 0),
    axis.text = element_text(size = 6),
    panel.spacing = unit(0.1, "lines"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
    legend.position = "right"
  )

# Save plot output as png
png_out <- file.path(out_dir, "hap1_vs_hap2_top11_dotplot.png")
ggsave(png_out, plot = dotplot, width = 12, height = 12, dpi = 300)
cat("Saved dotplot image to:", png_out, "\n")

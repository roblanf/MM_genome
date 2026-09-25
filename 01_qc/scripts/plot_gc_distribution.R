#!/usr/bin/env Rscript

# Libraries
# installs the packages you need in case you don't have them
# Auto-install missing packages to local user library (~/R)
required_packages <- c("ggplot2", "scales")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) {
  install.packages(new_packages, repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

#input/outputs
input_file <- "01_qc/results/gc/gc_pcts.tsv"
output_plot <- "01_qc/results/gc/gc_distribution.png"

# Read the first 500,000 rows (columns: read_id, length, gc) not all of them to save power
# seqkit fx2tab outputs GC as a percentage (0-100)
data <- read.delim(
  input_file, 
  header = FALSE, 
  nrows = 500000,
  col.names = c("read_id", "length", "gc")
)

# Calculate mean GC for vertical reference line
mean_gc <- mean(data$gc, na.rm = TRUE)

# Generate GC distribution histogram
p <- ggplot(data, aes(x = gc)) +
  geom_histogram(
    binwidth = 1, 
    fill = "#2b5c8f", 
    color = "white", 
    linewidth = 0.1
  ) +
  geom_vline(
    xintercept = mean_gc, 
    color = "#d95f02", 
    linetype = "dashed", 
    linewidth = 0.8
  ) +
  annotate(
    "text", 
    x = mean_gc + 2, 
    y = Inf, 
    label = paste0("Mean GC: ", round(mean_gc, 1), "%"), 
    vjust = 2, 
    hjust = 0, 
    color = "#d95f02", 
    fontface = "bold"
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 5), 
    limits = c(0, 100)
  ) +
  scale_y_continuous(labels = comma) +
  labs(
    title = expression(italic("Eucalyptus x phylacis") ~ "- Per-read GC Content Distribution"),
    subtitle = "Subsampled to first 500,000 reads",
    x = "GC Content (%)",
    y = "Number of Reads"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
    axis.title = element_text(face = "bold")
  )

# Save high-resolution plot
ggsave(
  filename = output_plot, 
  plot = p, 
  width = 8, 
  height = 5, 
  dpi = 300
)

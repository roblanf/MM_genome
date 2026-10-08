#!/usr/bin/env Rscript
# Purpose: New approach similar to initial_trials smoothing, rolling mean for AAACCCT telomere repeat profiles

library(tidyverse)
#for rolling mean:
library(zoo)

tidk_dir <- "03_assembly/results/tidk"
out_dir  <- "03_assembly/results/tidk/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

assemblies <- c("primary.p_ctg", "hap1.p_ctg", "hap2.p_ctg")

motif <- "AAACCCT"
motif_color <- "#228B22"

for (asm in assemblies) {
  file_path <- file.path(tidk_dir, sprintf("MM_assembly.%s_%s_telomeric_repeat_windows.tsv", asm, motif))
  
  if (file.exists(file_path)) {
    df <- read_tsv(file_path, show_col_types = FALSE)
    
    # Prefix mapping for hifiasm contig names to p, h1, h2
    prefix <- if (grepl("hap1", asm)) "h1" else if (grepl("hap2", asm)) "h2" else "p"
    target_contigs <- sprintf("%stg0000%02dl", prefix, 1:11)
    
    # Filter data for top 11 contigs and calculate total count & smoothed rolling mean for trendline
    plot_data <- df %>%
      filter(id %in% target_contigs) %>%
      mutate(
        total_count = forward_repeat_number + reverse_repeat_number,
        position_mb = window / 1e6,
        contig_label = factor(id, levels = target_contigs, labels = paste0("C", 1:11))
      ) %>%
      group_by(contig_label) %>%
      mutate(
        smoothed = rollmean(total_count, k = 5, fill = NA, align = "center")
      ) %>%
      # Handling NAs at boundaries for smoothing
      mutate(smoothed = ifelse(is.na(smoothed), total_count, smoothed)) %>%
      ungroup()
    
    p <- ggplot(plot_data, aes(x = position_mb)) +
      # 1. Raw count points
      geom_point(aes(y = total_count), color = "black", size = 0.4, alpha = 0.5) +
      # 2. 5-window smoothed rolling mean trend line
      geom_line(aes(y = smoothed), color = motif_color, linewidth = 0.8) +
      facet_wrap(~ contig_label, scales = "free_x", ncol = 1, strip.position = "left") +
      theme_classic(base_size = 11) +
      theme(
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 8),
        axis.text.y = element_text(size = 7),
        panel.spacing = unit(0.3, "lines"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, face = "italic", hjust = 0.5)
      ) +
      labs(
        title = sprintf("Repeat distribution for %s (%s)", asm, motif),
        subtitle = "Method: Motif count mapped across 10kb windows using 'tidk search' with smoothed rolling mean trendline",
        x = "Position (Mb)",
        y = "Motif count"
      )
    
    out_png <- file.path(out_dir, sprintf("%s_%s_fingerprint.png", asm, motif))
    ggsave(out_png, plot = p, width = 8, height = 10, dpi = 300)
    cat(sprintf("Saved updated fingerprint plot: %s\n", out_png))
  } else {
    cat(sprintf("Warning: File not found: %s\n", file_path))
  }
}

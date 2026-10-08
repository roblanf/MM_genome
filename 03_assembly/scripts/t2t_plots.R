#!/usr/bin/env Rscript
# Purpose: Uses smoothed trendline with rolling mean to get AAACCCT telomere repeat profiles for top11 longest contigs in each assembly

library(tidyverse)
#for rolling mean:
library(zoo)

tidk_dir <- "03_assembly/results/tidk"
stats_dir <- "03_assembly/results/hifiasm/post_hifi_stats"
out_dir  <- "03_assembly/results/tidk/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

assemblies <- c("primary.p_ctg", "hap1.p_ctg", "hap2.p_ctg")

motif <- "AAACCCT"

assemblies_config <- tibble(
  asm_id    = c("primary.p_ctg", "hap1.p_ctg", "hap2.p_ctg"),
  id_prefix = c("primary", "hap1", "hap2"),
  asm_label = c("Primary Assembly", "Haplotype 1 (Hap1)", "Haplotype 2 (Hap2)"),
  line_col  = c("#D97706", "#228B22", "#2563EB")
)

for (i in 1:nrow(assemblies_config)) {
  asm      <- assemblies_config$asm_id[i]
  prefix   <- assemblies_config$id_prefix[i]
  label    <- assemblies_config$asm_label[i]
  line_col <- assemblies_config$line_col[i]

  tsv_path <- file.path(tidk_dir, sprintf("MM_assembly.%s_%s_telomeric_repeat_windows.tsv", asm, motif)) 
  id_file  <- file.path(stats_dir, sprintf("top11_%s_ids.txt", prefix))
    
  if (file.exists(tsv_path) && file.exists(id_file)) {
    # Read 11 longest contig IDs in exact order from get_top11.sh
    top11_ids <- readLines(id_file) %>% trimws()
    top11_ids <- top11_ids[top11_ids != ""]

    # Map IDs to ranked labels C1...C11
    rank_map <- setNames(paste0("C", 1:length(top11_ids)), top11_ids)
    
    df <- read_tsv(tsv_path, show_col_types = FALSE)

    # Filter, calculate rolling mean, and enforce C1...C11 ordering
    plot_df <- df %>%
      filter(id %in% top11_ids) %>%
      mutate(total_count = forward_repeat_number + reverse_repeat_number) %>%
      group_by(id) %>%
      mutate(
        smoothed = rollmean(total_count, k = 5, fill = NA, align = "center"),
        smoothed = ifelse(is.na(smoothed), total_count, smoothed)
      ) %>%
      ungroup() %>%
      mutate(
        contig_label = factor(rank_map[id], levels = paste0("C", 1:11)),
        position_mb  = window / 1e6
      )
    
    # Plot for different assemblies
    p <- ggplot(plot_df, aes(x = position_mb)) +
      # 5-window smoothed rolling mean trend line
      geom_line(aes(y = smoothed), color = line_col, linewidth = 0.8) +
      # Raw count points on top
      geom_point(aes(y = smoothed), color = line_col, size = 0.4, alpha = 0.5) +
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
    cat(sprintf("Warning: Input files not found for %s. Checked:\n  TSV: %s\n  IDs: %s\n", asm, tsv_path, id_file))
  }
}

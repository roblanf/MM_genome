#!/usr/bin/env Rscript
# Purpose: Generate side-by-side comparative T2T telomere fingerprint plot for presentation (Primary vs Hap1 vs Hap2)

library(tidyverse)
library(zoo)

tidk_dir <- "03_assembly/results/tidk"
out_dir  <- "03_assembly/results/tidk/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

motif <- "AAACCCT"

# Set up for labels, colours, prefixes used in side-by-side plot
assemblies_config <- tibble(
  asm_id     = c("primary.p_ctg", "hap1.p_ctg", "hap2.p_ctg"),
  asm_label  = factor(c("Primary", "Haplotype 1 (Hap1)", "Haplotype 2 (Hap2)"),
                      levels = c("Primary", "Haplotype 1 (Hap1)", "Haplotype 2 (Hap2)")),
  prefix     = c("p", "h1", "h2"),
  line_color = c("#D97706", "#228B22", "#2563EB") # Orange (Primary), Green (Hap1), Blue (Hap2)
)

combined_data_list <- list()

for (i in 1:nrow(assemblies_config)) {
  asm        <- assemblies_config$asm_id[i]
  asm_disp   <- assemblies_config$asm_label[i]
  pfx        <- assemblies_config$prefix[i]

# Now most of the code from the t2t_plots.R script can be re-used  
  file_path <- file.path(tidk_dir, sprintf("MM_assembly.%s_%s_telomeric_repeat_windows.tsv", asm, motif))
  
  if (file.exists(file_path)) {
    df <- read_tsv(file_path, show_col_types = FALSE)
    
    target_contigs <- sprintf("%stg0000%02dl", pfx, 1:11)
    
    asm_data <- df %>%
      filter(id %in% target_contigs) %>%
      mutate(
        total_count   = forward_repeat_number + reverse_repeat_number,
        position_mb   = window / 1e6,
        contig_label  = factor(id, levels = target_contigs, labels = paste0("C", 1:11)),
        assembly_name = asm_disp
      ) %>%
      group_by(contig_label) %>%
      mutate(
        smoothed = rollmean(total_count, k = 5, fill = NA, align = "center"),
        smoothed = ifelse(is.na(smoothed), total_count, smoothed)
      ) %>%
      ungroup()
    
    combined_data_list[[asm]] <- asm_data
  } else {
    cat(sprintf("Warning: File not found: %s\n", file_path))
  }
}

plot_df <- bind_rows(combined_data_list)

p <- ggplot(plot_df, aes(x = position_mb)) +
  # 1. Smoothed line in background (colored by Assembly)
  geom_line(aes(y = smoothed, color = assembly_name), linewidth = 0.7, alpha = 0.85) +
  # 2. Raw black points on top (was hidden in t2t_plots because of wrong ordering)
  geom_point(aes(y = total_count), color = "black", size = 0.35, alpha = 0.5) +
  # Facet grid: Contigs on y-axis, assembly on x-axis
  facet_grid(contig_label ~ assembly_name, scales = "free_x") +
  scale_color_manual(values = setNames(assemblies_config$line_color, assemblies_config$asm_label)) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#F3F4F6", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    strip.text.y = element_text(angle = 0, face = "bold", size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 6),
    panel.spacing = unit(0.2, "lines"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "italic", hjust = 0.5)
  ) +
  labs(
    title = sprintf("Repeat distributions for primary, haplotype assemblies of telomere (%s)", motif),
    subtitle = "Motif count for top 11 contigs mapped accross 10Kb windows with smoothed rolling mean trendline",
    x = "Position (Mb)",
    y = "Motif Count"
  )

out_png <- file.path(out_dir, "t2t_telomere_side_by_side_comparison.png")
ggsave(out_png, plot = p, width = 14, height = 11, dpi = 300)
cat(sprintf("Saved presentation plot: %s\n", out_png))

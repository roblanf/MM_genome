#!/usr/bin/env Rscript
# Purpose: Single-motif plot for AAACCCT to visually confirm T2T contigs

library(tidyverse)

tidk_dir <- "03_assembly/results/tidk"
out_dir  <- "03_assembly/results/tidk/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

assemblies <- c("primary.p_ctg", "hap1.p_ctg", "hap2.p_ctg")

for (asm in assemblies) {
  file_path <- file.path(tidk_dir, sprintf("MM_assembly.%s_AAACCCT_telomeric_repeat_windows.tsv", asm))
  # Look for the key metrics needed for plotting in tsv files
  if (file.exists(file_path)) {
    dat <- read_tsv(file_path, show_col_types = FALSE)
    colnames(dat)[1:4] <- c("contig", "start", "end", "count")
    
    # Select top 11 largest contigs (n = 11 for Eucalyptus) and arrange for plotting
    top_contigs <- dat %>%
      group_by(contig) %>%
      summarise(max_len = max(end)) %>%
      arrange(desc(max_len)) %>%
      slice_head(n = 11) %>%
      pull(contig)
        
    plot_data <- dat %>%
      filter(contig %in% top_contigs) %>%
      mutate(
        position_mb = start / 1e6,
        contig_label = factor(contig, levels = top_contigs, labels = paste0("C", 1:length(top_contigs)))
      )
    # Use facet wrap to make plots for 11 top contigs
    p <- ggplot(plot_data, aes(x = position_mb, y = count)) +
      geom_line(color = "#D97706", linewidth = 0.5) +
      facet_wrap(~ contig_label, scales = "free_x", ncol = 1, strip.position = "left") +
      theme_classic(base_size = 11) +
      theme(
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0.3, "lines")
      ) +
      labs(
        title = sprintf("Telomere Repeat Distribution Profile: %s", asm),
        subtitle = "AAACCCT motif identified using `tidk explore`, mapped using `tidk search` with 10Kb windows (bins)",
        x = "Position (Mb)",
        y = "Motif occurrences"
      )
    
    out_png <- file.path(out_dir, sprintf("%s_t2t_aaaccct.png", asm))
    ggsave(out_png, plot = p, width = 8, height = 10, dpi = 300)
    cat(sprintf("Saved T2T plot: %s\n", out_png))
  }
}

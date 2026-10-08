#!/usr/bin/env Rscript
# Purpose: Single-motif plot for AAACCCT to visually confirm T2T contigs

library(tidyverse)

tidk_dir <- "03_assembly/results/tidk"
out_dir  <- "03_assembly/results/tidk/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

assemblies <- c("primary.p_ctg", "hap1.p_ctg", "hap2.p_ctg")

#Look at both forward and reverse telomere so we get peaks at both ends, not just one
for (asm in assemblies) {
  fwd_file <- file.path(tidk_dir, sprintf("MM_assembly.%s_AAACCCT_telomeric_repeat_windows.tsv", asm))
  rev_file <- file.path(tidk_dir, sprintf("MM_assembly.%s_TTTAGGG_telomeric_repeat_windows.tsv", asm))
  
# Look for the key metrics needed for plotting in tsv files
# Updated with numeric parsing from tsv files as error thought to be from issues in count_fwd + count_rev before
  if (file.exists(fwd_file) && file.exists(rev_file)) {
    # Skip header/comment lines in case that was an issue
    fwd <- read_tsv(fwd_file, comment = "#", show_col_types = FALSE)
    rev <- read_tsv(rev_file, comment = "#", show_col_types = FALSE)
    
    # Standardize first 4 columns: contig, start, end, count
    colnames(fwd)[1:4] <- c("contig", "start", "end", "count_fwd")
    colnames(rev)[1:4] <- c("contig", "start", "end", "count_rev")

    # Filter out non-numeric header rows if present and parse counts
    fwd_clean <- fwd %>%
      filter(!is.na(suppressWarnings(as.numeric(start)))) %>%
      mutate(start = as.numeric(start), end = as.numeric(end), count_fwd = as.numeric(count_fwd))
      
    rev_clean <- rev %>%
      filter(!is.na(suppressWarnings(as.numeric(start)))) %>%
      mutate(start = as.numeric(start), end = as.numeric(end), count_rev = as.numeric(count_rev))

    # Merge forward & reverse counts
    combined <- fwd_clean %>%
      inner_join(rev_clean, by = c("contig", "start", "end")) %>%
      mutate(count = coalesce(count_fwd, 0) + coalesce(count_rev, 0))
    
    # Focus on top 11 contigs, sort by length, descending to arrange for plotting
    top_contigs <- combined %>%
      group_by(contig) %>%
      summarise(max_len = max(end, na.rm = TRUE)) %>%
      arrange(desc(max_len)) %>%
      slice_head(n = 11) %>%
      pull(contig)
    
    plot_data <- combined %>%
      filter(contig %in% top_contigs) %>%
      mutate(
        position_mb = start / 1e6,
        contig_label = factor(contig, levels = top_contigs, labels = paste0("C", 1:length(top_contigs)))
      )

    # Use facet wrap to make plots for 11 top contigs
    # Fixed ylim to avoid cutting out the real end peaks
    p <- ggplot(plot_data, aes(x = position_mb, y = count)) +
      geom_line(color = "#D97706", linewidth = 0.4) +
      ylim(0, max(plot_data$count)) +
      facet_wrap(~ contig_label, scales = "free_x", ncol = 1, strip.position = "left") +
      theme_classic(base_size = 11) +
      theme(
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0.3, "lines"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9)
      ) +
      labs(
        title = sprintf("Telomere Repeat Distribution Profile: %s", asm),
        subtitle = "AAACCCT telomere motif identified using `tidk explore`, mapped using `tidk search` with 10Kb windows (bins)",
        x = "Position (Mb)",
        y = "Motif count per window"
      )
    
    out_png <- file.path(out_dir, sprintf("%s_t2t_aaaccct.png", asm))
    ggsave(out_png, plot = p, width = 8, height = 10, dpi = 300)
    cat(sprintf("Saved T2T plot: %s\n", out_png))
  }
}

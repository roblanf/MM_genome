#!/usr/bin/env Rscript

library(tidyverse)
library(viridis)

# --- 1. Arguments & Setup ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript generate_parental_plots.R <input_dir> <output_dir>; expects 4 specifically named input files in input dir, check script for names")
}

input_dir  <- args[1]
output_dir <- args[2]

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message(paste("Processing data from:", input_dir))
message(paste("Saving results to:", output_dir))

# --- 2. The Neutral Lookup Table ---
lookup <- tribble(
  ~chr_num, ~hap1_contig, ~hap2_contig,
  "Chr01", "h1tg000007l", "h2tg000010l",
  "Chr02", "h1tg000010l", "h2tg000008l",
  "Chr03", "h2tg000009l", "h1tg000006l",
  "Chr04", "h2tg000011l", "h1tg000004l",
  "Chr05", "h1tg000003l", "h2tg000307l",
  "Chr06", "h1tg000011l", "h2tg000001l",
  "Chr07", "h2tg000002l", "h1tg000001l",
  "Chr08", "h1tg000008l", "h2tg000003l",
  "Chr09", "h2tg000006l", "h1tg000009l",
  "Chr10", "h1tg000002l", "h2tg000004l",
  "Chr11", "h1tg000005l", "h2tg000007l"
)

# --- 3. Loading & Processing Functions ---
load_and_tidy <- function(v_file, d_file, assembly_label, window_label) {
  v_path <- file.path(input_dir, v_file)
  d_path <- file.path(input_dir, d_file)
  
  cols <- c("contig", "start", "end", "kmer_count", "bases_covered", "window_size", "fraction")
  
  v <- read_tsv(v_path, col_names = cols, show_col_types = FALSE) %>% 
    mutate(parent = "E. virginea", assembly = assembly_label, window = window_label)
  
  d <- read_tsv(d_path, col_names = cols, show_col_types = FALSE) %>% 
    mutate(parent = "E. decipiens", assembly = assembly_label, window = window_label)
  
  bind_rows(v, d) %>%
    left_join(
      bind_rows(
        lookup %>% select(chr_num, contig = hap1_contig),
        lookup %>% select(chr_num, contig = hap2_contig)
      ) %>% distinct(),
      by = "contig"
    ) %>%
    filter(!is.na(chr_num))
}

calculate_proportions <- function(data_set) {
  data_set %>%
    select(chr_num, assembly, start, end, parent, kmer_count) %>%
    pivot_wider(names_from = parent, values_from = kmer_count) %>%
    rename(v_hits = `E. virginea`, d_hits = `E. decipiens`) %>%
    mutate(
      total = v_hits + d_hits,
      `E. virginea` = ifelse(total > 0, v_hits / total, 0.5),
      `E. decipiens` = ifelse(total > 0, d_hits / total, 0.5)
    ) %>%
    pivot_longer(cols = c(`E. virginea`, `E. decipiens`), 
                 names_to = "parent", 
                 values_to = "proportion")
}

# --- 4. Load and Prepare Data ---
message("Reading raw BED coverage files...")
data_10k <- bind_rows(
  load_and_tidy("h1_v_10k.txt", "h1_d_10k.txt", "Hap1", "10K"),
  load_and_tidy("h2_v_10k.txt", "h2_d_10k.txt", "Hap2", "10K")
)

data_100k <- bind_rows(
  load_and_tidy("h1_v_100k.txt", "h1_d_100k.txt", "Hap1", "100K"),
  load_and_tidy("h2_v_100k.txt", "h2_d_100k.txt", "Hap2", "100K")
)

message("Calculating proportions...")
data_prop_10k  <- calculate_proportions(data_10k)
data_prop_100k <- calculate_proportions(data_100k)

# --- 5. Plotting Functions ---

# Function for Raw Counts (Free Y-axis)
plot_raw <- function(df, window_lab) {
  ggplot(df, aes(x = start / 1e6, y = kmer_count, color = parent)) +
    geom_line(alpha = 1.0, linewidth = ifelse(window_lab == "10K", 0.2, 1.0)) +
    geom_point(alpha = 0.5, size = ifelse(window_lab == "10K", 0.2, 0.8)) +
    facet_grid(chr_num ~ assembly, scales = "free") +
    scale_color_viridis_d(option = "viridis", end = 0.8) +
    theme_minimal(base_size = 10) +
    labs(title = paste("Raw Parental K-mer Density (", window_lab, " Windows)", sep=""),
         x = "Position (Mb)", y = "K-mer Hits per Window", color = "Parent") +
    theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))
}

# Function for Proportions (Fixed 0-1 Y-axis)
plot_prop <- function(df, window_lab) {
  ggplot(df, aes(x = start / 1e6, y = proportion, color = parent)) +
    geom_line(alpha = 1.0, linewidth = ifelse(window_lab == "10K", 0.2, 1.0)) +
    geom_point(alpha = 0.3, size = ifelse(window_lab == "10K", 0.2, 0.8)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.5) +
    facet_grid(chr_num ~ assembly, scales = "free_x") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_color_viridis_d(option = "viridis", end = 0.8) +
    theme_minimal(base_size = 10) +
    labs(title = paste("Parental K-mer Proportion (", window_lab, " Windows)", sep=""),
         subtitle = "0.5 = neutral/shared. Values > 0.5 indicate parental dominance.",
         x = "Position (Mb)", y = "Proportion of Total Hits", color = "Parent") +
    theme(legend.position = "bottom", strip.text.y = element_text(angle = 0), panel.grid.minor = element_blank())
}

# --- 6. Generate and Save All Plots ---
message("Generating and saving plots...")

# Common dimensions for 11-row layout
w <- 24
h <- 16

# 1. Raw 10K
ggsave(file.path(output_dir, "raw_density_10k.png"), plot_raw(data_10k, "10K"), width = w, height = h, dpi = 300)

# 2. Raw 100K
ggsave(file.path(output_dir, "raw_density_100k.png"), plot_raw(data_100k, "100K"), width = w, height = h, dpi = 300)

# 3. Proportion 10K
ggsave(file.path(output_dir, "proportion_10k.png"), plot_prop(data_prop_10k, "10K"), width = w, height = h, dpi=300)

# 4. Proportion 100K
ggsave(file.path(output_dir, "proportion_100k.png"), plot_prop(data_prop_100k, "100K"), width = w, height = h, dpi = 300)

message("Success! All plots are in the output directory.")

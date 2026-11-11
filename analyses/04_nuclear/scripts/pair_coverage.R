#!/usr/bin/env Rscript

# ============================================================
#  Pairwise coverage between the contigs of two haplotypes
#  from a minimap2 PAF
#
#  Input : PAF of hap1 (query) aligned to hap2 (target), e.g.
#          minimap2 -cx asm20 --secondary=no -t 128 hap2.purged.fa hap1.purged.fa > h1_vs_h2.paf
#  Output: <prefix>.pairs.tsv    one row per (hap1 contig, hap2 contig) pair
#          <prefix>.contigs.tsv  one row per contig of either haplotype
#
#  Usage: Rscript pair_coverage.R -i h1_vs_h2.paf -p h1_vs_h2
#         Rscript pair_coverage.R -i h1_vs_h2.paf -p h1_vs_h2 --min-aln 10000 --min-id 0.9 --min-cov 0.01
# ============================================================

# --- Load libraries and functions ---------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

# compare and sort contig names byte-wise (same as Python)
invisible(Sys.setlocale("LC_COLLATE", "C"))

# disable scientific notation
options(scipen=999)

# function: extract an optional PAF tag (e.g. "de:f:") from every line
f_extract_tag <- function(ls_lines, tag) {
  match <- regexpr(paste0("\t", tag, "[^\t]*"), ls_lines, perl=TRUE)

  ls_value <- rep(NA_character_, length(ls_lines))
  ls_value[match != -1] <- substring(regmatches(ls_lines, match), nchar(tag) + 2)

  return(ls_value)
}

# function: read a PAF file into a data.table (12 standard columns + de and tp tags)
f_read_paf <- function(fn_paf) {
  ls_lines <- readLines(fn_paf)

  # keep lines with at least 12 columns
  ls_lines <- ls_lines[grepl("^([^\t]*\t){11}", ls_lines, perl=TRUE)]

  # standard columns
  df_paf <- data.table::fread(text=sub("^((?:[^\t]*\t){11}[^\t]*).*$", "\\1", ls_lines, perl=TRUE),
                              sep="\t", header=FALSE, quote="",
                              col.names=c("qname", "qlen", "qstart", "qend", "strand",
                                          "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq"))

  # optional tags
  df_paf$de <- as.numeric(f_extract_tag(ls_lines, "de:f:"))
  df_paf$tp <- f_extract_tag(ls_lines, "tp:A:")

  return(df_paf)
}

# function: total length covered by a set of (possibly overlapping) intervals
f_merged_len <- function(ls_start, ls_end) {
  idx <- order(ls_start, ls_end)
  ls_start <- ls_start[idx]
  ls_end <- ls_end[idx]

  # an interval starts a new block if it begins after every previous end
  ls_prev_end <- c(-Inf, cummax(ls_end)[-length(ls_end)])
  ls_block <- cumsum(ls_start > ls_prev_end)

  return(sum(tapply(ls_end, ls_block, max) - tapply(ls_start, ls_block, min)))
}

# function: "partner:cov" list of every partner covering >= min_cov, best first
f_partner_list <- function(ls_partner, ls_cov, min_cov) {
  idx <- order(-ls_cov)
  idx <- idx[ls_cov[idx] >= min_cov]

  return(paste(sprintf("%s:%.2f", ls_partner[idx], ls_cov[idx]), collapse=","))
}

# function: convert TRUE/FALSE into yes/no
f_yes_no <- function(x) {
  return(ifelse(x, "yes", "no"))
}

# --- Argument parsing ----------------------------------------
option_list <- list(
  make_option(c("-i", "--paf"), type="character", default=NULL,
              help="PAF of hap1 (query) aligned to hap2 (target) [required]", metavar="FILE"),
  make_option(c("-p", "--prefix"), type="character", default="pairs",
              help="output prefix [default %default]", metavar="PREFIX"),
  make_option(c("--min-aln"), dest="min_aln", type="integer", default=10000,
              help="min alignment block length in bp [default %default]"),
  make_option(c("--min-id"), dest="min_id", type="double", default=0.90,
              help="min identity (1 - de, or matches/alen) [default %default]"),
  make_option(c("--min-cov"), dest="min_cov", type="double", default=0.01,
              help="report pairs covering >= this fraction of either contig [default %default]")
)

# parse the arguments
opt <- parse_args(OptionParser(option_list=option_list))

# stop the code if input file is invalid
if (is.null(opt$paf)) {
  stop("-i/--paf is required.\n
        Usage: Rscript pair_coverage.R -i h1_vs_h2.paf -p h1_vs_h2")
}

if (!file.exists(opt$paf)) {
  stop(paste("PAF file not found:", opt$paf))
}

# --- Output files --------------------------------------------
fn_pairs   <- paste0(opt$prefix, ".pairs.tsv")
fn_contigs <- paste0(opt$prefix, ".contigs.tsv")

# --- Read and filter alignments ------------------------------
df_paf <- f_read_paf(opt$paf)

# remove short alignments
df_paf <- df_paf[alen >= opt$min_aln]

# identity: 1 - de tag if present, else matches / alignment length
df_paf$identity <- ifelse(is.na(df_paf$de), df_paf$nmatch / df_paf$alen, 1 - df_paf$de)
df_paf <- df_paf[identity >= opt$min_id]

# --- Summarise per pair of contigs ---------------------------
df_pairs <- df_paf[, .(h1_len     = qlen[1],
                       h2_len     = tlen[1],
                       n_aln      = .N,
                       aligned_bp = sum(alen),
                       identity   = sum(identity * alen) / sum(alen),
                       h1_cov     = f_merged_len(qstart, qend) / qlen[1],
                       h2_cov     = f_merged_len(tstart, tend) / tlen[1]),
                   by=.(h1=qname, h2=tname)]

# best partner of every contig (highest coverage of the contig itself)
df_pairs[, best_for_h1 := h2 == h2[which.max(h1_cov)], by=h1]
df_pairs[, best_for_h2 := h1 == h1[which.max(h2_cov)], by=h2]
df_pairs[, reciprocal := best_for_h1 & best_for_h2]

# --- Write pairs table ---------------------------------------
df_out <- df_pairs[order(-aligned_bp)]
df_out <- df_out[h1_cov >= opt$min_cov | h2_cov >= opt$min_cov]
df_out <- df_out[, .(h1, h1_len, h2, h2_len, n_aln, aligned_bp, identity, h1_cov, h2_cov,
                     best_for_h1, best_for_h2, reciprocal)]

df_out$identity    <- sprintf("%.4f", df_out$identity)
df_out$h1_cov      <- sprintf("%.3f", df_out$h1_cov)
df_out$h2_cov      <- sprintf("%.3f", df_out$h2_cov)
df_out$best_for_h1 <- f_yes_no(df_out$best_for_h1)
df_out$best_for_h2 <- f_yes_no(df_out$best_for_h2)
df_out$reciprocal  <- f_yes_no(df_out$reciprocal)

data.table::fwrite(df_out, file=fn_pairs, sep="\t", quote=FALSE)

# --- Summarise per contig ------------------------------------
# fraction of each contig covered by the other haplotype (all partners)
df_total_h1 <- df_paf[, .(total_cov=f_merged_len(qstart, qend) / qlen[1]), by=.(contig=qname)]
df_total_h2 <- df_paf[, .(total_cov=f_merged_len(tstart, tend) / tlen[1]), by=.(contig=tname)]

# hap1 contigs
df_h1 <- df_pairs[, .(hap              = "hap1",
                      length           = h1_len[1],
                      best_partner     = h2[which.max(h1_cov)],
                      best_partner_cov = max(h1_cov),
                      reciprocal       = reciprocal[which.max(h1_cov)],
                      partners         = f_partner_list(h2, h1_cov, opt$min_cov)),
                  by=.(contig=h1)]
df_h1$total_cov <- df_total_h1$total_cov[match(df_h1$contig, df_total_h1$contig)]

# hap2 contigs
df_h2 <- df_pairs[, .(hap              = "hap2",
                      length           = h2_len[1],
                      best_partner     = h1[which.max(h2_cov)],
                      best_partner_cov = max(h2_cov),
                      reciprocal       = reciprocal[which.max(h2_cov)],
                      partners         = f_partner_list(h1, h2_cov, opt$min_cov)),
                  by=.(contig=h2)]
df_h2$total_cov <- df_total_h2$total_cov[match(df_h2$contig, df_total_h2$contig)]

# combine and sort by haplotype and length
df_contigs <- rbind(df_h1, df_h2)
df_contigs <- df_contigs[order(hap, -length)]
df_contigs <- df_contigs[, .(contig, hap, length, total_cov, best_partner, best_partner_cov, reciprocal, partners)]

# --- Write contigs table -------------------------------------
df_contigs$total_cov        <- sprintf("%.3f", df_contigs$total_cov)
df_contigs$best_partner_cov <- sprintf("%.3f", df_contigs$best_partner_cov)
df_contigs$reciprocal       <- f_yes_no(df_contigs$reciprocal)

data.table::fwrite(df_contigs, file=fn_contigs, sep="\t", quote=FALSE)

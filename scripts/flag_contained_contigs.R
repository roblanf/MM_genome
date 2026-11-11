#!/usr/bin/env Rscript

# ============================================================
#  Flag contigs that are (mostly) contained in a LONGER contig
#  of the same haplotype
#
#  Input : PAF from a self-alignment, e.g.
#          minimap2 -cx asm20 -DP -m200 -t16 hap1.fa hap1.fa > h1_self.paf
#  Output: TSV on stdout, one row per flagged contig (most-covered first)
#          contig  length  frac_covered  max_identity  contained_in
#
#  Usage: Rscript flag_contained_contigs.R -i h1_self.paf > h1_dups.tsv
#         Rscript flag_contained_contigs.R -i h1_self.paf --min-aln 10000 --min-id 0.9 --min-frac 0.8
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

# function: the (up to) n targets receiving the most aligned bases
#           (alignments are first sorted by position, so ties are broken as in the Python script)
f_top_targets <- function(ls_start, ls_end, ls_identity, ls_target, n=5) {
  idx <- order(ls_start, ls_end, ls_identity, ls_target)

  df_target <- data.table::data.table(target=ls_target[idx], len=(ls_end - ls_start)[idx])
  df_target <- df_target[, .(len=sum(len)), by=target]
  df_target <- df_target[order(-len)]

  return(paste(head(df_target$target, n), collapse=","))
}

# --- Argument parsing ----------------------------------------
option_list <- list(
  make_option(c("-i", "--paf"), type="character", default=NULL,
              help="PAF from a minimap2 self-alignment [required]", metavar="FILE"),
  make_option(c("--min-aln"), dest="min_aln", type="integer", default=10000,
              help="min alignment block length in bp [default %default]"),
  make_option(c("--min-id"), dest="min_id", type="double", default=0.90,
              help="min identity, 1 - de tag or matches/alen [default %default]"),
  make_option(c("--min-frac"), dest="min_frac", type="double", default=0.80,
              help="flag contig if >= this fraction of its length is covered [default %default]")
)

# parse the arguments
opt <- parse_args(OptionParser(option_list=option_list))

# stop the code if input file is invalid
if (is.null(opt$paf)) {
  stop("-i/--paf is required.\n
        Usage: Rscript flag_contained_contigs.R -i h1_self.paf > h1_dups.tsv")
}

if (!file.exists(opt$paf)) {
  stop(paste("PAF file not found:", opt$paf))
}

# --- Read and filter alignments ------------------------------
df_paf <- f_read_paf(opt$paf)

# remove self-hits and short alignments
df_paf <- df_paf[qname != tname & alen >= opt$min_aln]

# only count alignments to a longer contig, so the shorter copy of a
# duplicated pair is flagged, not both (ties broken by name)
df_paf <- df_paf[tlen > qlen | (tlen == qlen & tname > qname)]

# identity: 1 - de tag if present, else matches / alignment length
df_paf$identity <- ifelse(is.na(df_paf$de), df_paf$nmatch / df_paf$alen, 1 - df_paf$de)
df_paf <- df_paf[identity >= opt$min_id]

# --- Summarise per contig ------------------------------------
df_contig <- df_paf[, .(length       = qlen[1],
                        frac_covered = f_merged_len(qstart, qend) / qlen[1],
                        max_identity = max(identity),
                        contained_in = f_top_targets(qstart, qend, identity, tname)),
                    by=.(contig=qname)]

# keep contigs that are mostly contained in another contig
df_contig <- df_contig[frac_covered >= opt$min_frac]
df_contig <- df_contig[order(-frac_covered)]

# --- Output --------------------------------------------------
df_contig$frac_covered <- sprintf("%.3f", df_contig$frac_covered)
df_contig$max_identity <- sprintf("%.4f", df_contig$max_identity)

data.table::fwrite(df_contig, file="", sep="\t", quote=FALSE)

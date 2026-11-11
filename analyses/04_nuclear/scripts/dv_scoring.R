#!/usr/bin/env Rscript

# ============================================================
#  Parental assignment table from alignments of hap1/hap2 to
#  two candidate parent references
#
#  Every contig is scored on ALL its primary alignments to each
#  reference, not on the single alignment with the most matches.
#
#  Per contig and per reference:
#    aligned %  : fraction of the contig covered by primary alignments (merged intervals)
#    identity % : length-weighted identity over those alignments (1 - de tag if present, else matches/block)
#    chrom      : the reference chromosome receiving the most aligned bases
#  Derived: delta identity (decipiens - virginea), aligned ratio (decipiens / virginea), a call,
#  the partner contig in the other haplotype (largest summed alignment length in the hap1-vs-hap2 PAF),
#  the partner's call, and whether the two calls are opposite (reciprocal).
#
#  Usage: Rscript dv_scoring.R --h1h2 hap1_vs_hap2.paf \
#             --h1-dec hap1_vs_decipiens.paf --h1-vir hap1_vs_virginea.paf \
#             --h2-dec hap2_vs_decipiens.paf --h2-vir hap2_vs_virginea.paf \
#             [--blob1 hap1_stats.tsv --blob2 hap2_stats.tsv] \
#             [--min-aln 1000] [--min-id 0.80] [--margin 0.5] [--min-aligned 0.2] [--all-alignments] \
#             --md combined_parental_comparison.md --tsv combined_parental_comparison.tsv
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

# function: the target receiving the most aligned bases
f_top_target <- function(ls_target, ls_len) {
  df_target <- data.table::data.table(target=ls_target, len=ls_len)[, .(len=sum(len)), by=target]

  return(df_target$target[which.max(df_target$len)])
}

# function: score every contig of a PAF against one reference
f_score_paf <- function(fn_paf, min_aln, min_id, all_alignments) {
  df_paf <- f_read_paf(fn_paf)

  # keep primary alignments only (lines without a tp tag are kept)
  if (!all_alignments) {
    df_paf <- df_paf[is.na(tp) | tp == "P"]
  }

  # remove short and divergent alignments
  df_paf <- df_paf[alen >= min_aln]
  df_paf$identity <- ifelse(is.na(df_paf$de), df_paf$nmatch / df_paf$alen, 1 - df_paf$de)
  df_paf <- df_paf[identity >= min_id]

  # summarise per contig
  df_score <- df_paf[, .(len      = qlen[1],
                         aligned  = f_merged_len(qstart, qend) / qlen[1],
                         identity = 100 * sum(identity * alen) / sum(alen),
                         chrom    = f_top_target(tname, qend - qstart),
                         n        = .N),
                     by=.(contig=qname)]

  return(df_score)
}

# function: best partner of each hap1 contig (and of each hap2 contig) by summed alignment block length
f_partners <- function(fn_paf) {
  df_paf <- f_read_paf(fn_paf)

  df_sum <- df_paf[, .(alen=sum(alen)), by=.(qname, tname)]
  df_best_h1 <- df_sum[, .(partner=tname[which.max(alen)]), by=.(contig=qname)]
  df_best_h2 <- df_sum[, .(partner=qname[which.max(alen)]), by=.(contig=tname)]

  df_len_h1 <- unique(df_paf[, .(contig=qname, len=qlen)])
  df_len_h2 <- unique(df_paf[, .(contig=tname, len=tlen)])

  return(list(best_h1=df_best_h1, best_h2=df_best_h2, len_h1=df_len_h1, len_h2=df_len_h2))
}

# function: contig -> gc%, cov, phylum from a blobtools stats table; tolerant to column names/positions
f_read_blob <- function(fn_blob) {
  df_blob <- data.table::data.table(contig=character(), gc=numeric(), cov=numeric(), phylum=character())
  if (is.null(fn_blob)) {
    return(df_blob)
  }

  ls_rows <- strsplit(readLines(fn_blob), "\t", fixed=TRUE)
  ls_header <- tolower(trimws(ls_rows[[1]]))

  # function: index of the first column matching (or starting with) any of the names
  f_col <- function(ls_names, default) {
    for (i in seq_along(ls_header)) {
      if (any(ls_header[i] == ls_names | startsWith(ls_header[i], ls_names))) {
        return(i)
      }
    }
    return(default)
  }

  col_id     <- f_col(c("name", "id", "contig"), 2)
  col_gc     <- f_col(c("gc"), 3)
  col_cov    <- f_col(c("cov"), 5)
  col_phylum <- f_col(c("phylum", "bestsum_phylum"), 6)

  # keep complete rows (and skip a repeated header line)
  ls_rows <- ls_rows[-1]
  ls_rows <- ls_rows[sapply(ls_rows, function(row) {
    length(row) >= max(col_id, col_gc, col_cov, col_phylum) && tolower(trimws(row[1])) != "index"
  })]
  if (length(ls_rows) == 0) {
    return(df_blob)
  }

  df_blob <- data.table::data.table(contig = trimws(sapply(ls_rows, `[`, col_id)),
                                    gc     = suppressWarnings(as.numeric(sapply(ls_rows, `[`, col_gc))),
                                    cov    = suppressWarnings(as.numeric(sapply(ls_rows, `[`, col_cov))),
                                    phylum = trimws(sapply(ls_rows, `[`, col_phylum)))

  # skip rows with non-numeric values; gc as a fraction -> percentage
  df_blob <- df_blob[!is.na(gc) & !is.na(cov)]
  df_blob$gc <- ifelse(df_blob$gc <= 1, df_blob$gc * 100, df_blob$gc)

  return(df_blob)
}

# function: assign a contig to a parent
f_call <- function(dec_aligned, vir_aligned, dec_id, vir_id, min_aligned, margin) {
  has_dec <- !is.na(dec_aligned)
  has_vir <- !is.na(vir_aligned)

  if (!has_dec && !has_vir) {
    return("no alignment")
  }

  if (max(ifelse(has_dec, dec_aligned, 0), ifelse(has_vir, vir_aligned, 0)) < min_aligned) {
    return("unaligned")
  }

  if (!has_dec) {
    return("virginea")
  }
  if (!has_vir) {
    return("decipiens")
  }

  delta <- dec_id - vir_id
  if (delta >= margin) {
    return("decipiens")
  }
  if (delta <= -margin) {
    return("virginea")
  }

  return("ambiguous")
}

# function: format a number, or "-" if missing
f_fmt <- function(x, digits=2, sign=FALSE) {
  fmt <- paste0("%", ifelse(sign, "+", ""), ".", digits, "f")

  return(ifelse(is.na(x), "-", sprintf(fmt, x)))
}

# function: pad a string to a given width (character-aware)
f_pad <- function(x, width, right=FALSE) {
  padding <- strrep(" ", pmax(0, width - nchar(x)))

  return(ifelse(right, paste0(padding, x), paste0(x, padding)))
}

# --- Argument parsing ----------------------------------------
option_list <- list(
  make_option(c("--h1h2"), type="character", default=NULL,
              help="PAF of hap1 (query) vs hap2 (target) [required]", metavar="FILE"),
  make_option(c("--h1-dec"), dest="h1_dec", type="character", default=NULL,
              help="PAF of hap1 vs decipiens reference [required]", metavar="FILE"),
  make_option(c("--h1-vir"), dest="h1_vir", type="character", default=NULL,
              help="PAF of hap1 vs virginea reference [required]", metavar="FILE"),
  make_option(c("--h2-dec"), dest="h2_dec", type="character", default=NULL,
              help="PAF of hap2 vs decipiens reference [required]", metavar="FILE"),
  make_option(c("--h2-vir"), dest="h2_vir", type="character", default=NULL,
              help="PAF of hap2 vs virginea reference [required]", metavar="FILE"),
  make_option(c("--blob1"), type="character", default=NULL,
              help="blobtools stats table of hap1", metavar="FILE"),
  make_option(c("--blob2"), type="character", default=NULL,
              help="blobtools stats table of hap2", metavar="FILE"),
  make_option(c("--min-aln"), dest="min_aln", type="integer", default=1000,
              help="min alignment block length (bp) [default %default]"),
  make_option(c("--min-id"), dest="min_id", type="double", default=0.80,
              help="min identity of an alignment to be counted [default %default]"),
  make_option(c("--margin"), type="double", default=0.5,
              help="|delta identity| (percentage points) needed for a call [default %default]"),
  make_option(c("--min-aligned"), dest="min_aligned", type="double", default=0.2,
              help="min aligned fraction to either reference for a call [default %default]"),
  make_option(c("--all-alignments"), dest="all_alignments", action="store_true", default=FALSE,
              help="count secondary alignments too (default: primary only)"),
  make_option(c("--md"), type="character", default="combined_parental_comparison.md",
              help="output Markdown table [default %default]", metavar="FILE"),
  make_option(c("--tsv"), type="character", default="combined_parental_comparison.tsv",
              help="output TSV table [default %default]", metavar="FILE")
)

# parse the arguments
opt <- parse_args(OptionParser(option_list=option_list))

# stop the code if any required input is missing
ls_required <- c("h1h2", "h1_dec", "h1_vir", "h2_dec", "h2_vir")
ls_missing <- ls_required[sapply(ls_required, function(x) is.null(opt[[x]]))]
if (length(ls_missing) > 0) {
  stop(paste("Missing required arguments:", paste0("--", gsub("_", "-", ls_missing), collapse=", ")))
}

for (fn in c(opt$h1h2, opt$h1_dec, opt$h1_vir, opt$h2_dec, opt$h2_vir, opt$blob1, opt$blob2)) {
  if (!file.exists(fn)) {
    stop(paste("File not found:", fn))
  }
}

# --- Score the alignments ------------------------------------
df_h1_dec <- f_score_paf(opt$h1_dec, opt$min_aln, opt$min_id, opt$all_alignments)
df_h1_vir <- f_score_paf(opt$h1_vir, opt$min_aln, opt$min_id, opt$all_alignments)
df_h2_dec <- f_score_paf(opt$h2_dec, opt$min_aln, opt$min_id, opt$all_alignments)
df_h2_vir <- f_score_paf(opt$h2_vir, opt$min_aln, opt$min_id, opt$all_alignments)

ls_partners <- f_partners(opt$h1h2)

df_blob1 <- f_read_blob(opt$blob1)
df_blob2 <- f_read_blob(opt$blob2)

# --- Combine per haplotype -----------------------------------
df_rows <- data.table::data.table()

ls_haps <- list(hap1=list(dec=df_h1_dec, vir=df_h1_vir, partner=ls_partners$best_h1, len=ls_partners$len_h1, blob=df_blob1),
                hap2=list(dec=df_h2_dec, vir=df_h2_vir, partner=ls_partners$best_h2, len=ls_partners$len_h2, blob=df_blob2))

for (hap in names(ls_haps)) {
  df_dec     <- ls_haps[[hap]]$dec
  df_vir     <- ls_haps[[hap]]$vir
  df_partner <- ls_haps[[hap]]$partner
  df_len     <- ls_haps[[hap]]$len
  df_blob    <- ls_haps[[hap]]$blob

  # every contig seen in any of the inputs
  ls_contigs <- sort(unique(c(df_dec$contig, df_vir$contig, df_len$contig)))

  idx_dec  <- match(ls_contigs, df_dec$contig)
  idx_vir  <- match(ls_contigs, df_vir$contig)
  idx_len  <- match(ls_contigs, df_len$contig)
  idx_part <- match(ls_contigs, df_partner$contig)
  idx_blob <- match(ls_contigs, df_blob$contig)

  # contig length: from the decipiens scores, else virginea, else the hap1-vs-hap2 PAF, else 0
  length <- df_dec$len[idx_dec]
  length <- ifelse(is.na(length), df_vir$len[idx_vir], length)
  length <- ifelse(is.na(length), df_len$len[idx_len], length)
  length <- ifelse(is.na(length), 0, length)

  df_hap <- data.table::data.table(contig      = ls_contigs,
                                   hap         = hap,
                                   length      = length,
                                   gc          = df_blob$gc[idx_blob],
                                   cov         = df_blob$cov[idx_blob],
                                   phylum      = ifelse(is.na(idx_blob), "NA", df_blob$phylum[idx_blob]),
                                   dec_chrom   = ifelse(is.na(idx_dec), "-", df_dec$chrom[idx_dec]),
                                   vir_chrom   = ifelse(is.na(idx_vir), "-", df_vir$chrom[idx_vir]),
                                   dec_aligned = ifelse(is.na(idx_dec), 0, 100 * df_dec$aligned[idx_dec]),
                                   vir_aligned = ifelse(is.na(idx_vir), 0, 100 * df_vir$aligned[idx_vir]),
                                   dec_id      = df_dec$identity[idx_dec],
                                   vir_id      = df_vir$identity[idx_vir],
                                   partner     = ifelse(is.na(idx_part), "-", df_partner$partner[idx_part]))

  df_hap$delta <- df_hap$dec_id - df_hap$vir_id
  df_hap$ratio <- ifelse(!is.na(idx_dec) & !is.na(idx_vir) & df_hap$vir_aligned > 0,
                         df_hap$dec_aligned / df_hap$vir_aligned, NA)

  # call the parent of every contig
  df_hap$call <- sapply(seq_along(ls_contigs), function(i) {
    f_call(df_dec$aligned[idx_dec[i]], df_vir$aligned[idx_vir[i]],
           df_hap$dec_id[i], df_hap$vir_id[i], opt$min_aligned, opt$margin)
  })

  df_rows <- rbind(df_rows, df_hap)
}

# partner call and reciprocity
idx_partner <- match(df_rows$partner, df_rows$contig)
df_rows$partner_call  <- ifelse(is.na(idx_partner), "-", df_rows$call[idx_partner])
df_rows$partner_delta <- df_rows$delta[idx_partner]
df_rows$reciprocal    <- ifelse((df_rows$call == "decipiens" & df_rows$partner_call == "virginea") |
                                (df_rows$call == "virginea" & df_rows$partner_call == "decipiens"),
                                "yes", "no")

# order of the calls for the summary (first appearance)
ls_call_levels <- unique(df_rows$call)

# sort by haplotype and length
df_rows <- df_rows[order(hap, -length)]

# --- Format the table ----------------------------------------
ls_header <- c("Contig", "Hap", "Length", "GC%", "Cov", "Phylum", "decipiens chrom", "virginea chrom",
               "decipiens aligned %", "virginea aligned %", "decipiens ID %", "virginea ID %", "\u0394 ID",
               "Aligned ratio", "Call", "Partner", "Partner \u0394 ID", "Partner call", "Reciprocal")

df_cells <- data.table::data.table(df_rows$contig,
                                   df_rows$hap,
                                   formatC(df_rows$length, format="d", big.mark=","),
                                   f_fmt(df_rows$gc, 1),
                                   f_fmt(df_rows$cov, 1),
                                   df_rows$phylum,
                                   df_rows$dec_chrom,
                                   df_rows$vir_chrom,
                                   f_fmt(df_rows$dec_aligned, 1),
                                   f_fmt(df_rows$vir_aligned, 1),
                                   f_fmt(df_rows$dec_id),
                                   f_fmt(df_rows$vir_id),
                                   f_fmt(df_rows$delta, sign=TRUE),
                                   f_fmt(df_rows$ratio),
                                   df_rows$call,
                                   df_rows$partner,
                                   f_fmt(df_rows$partner_delta, sign=TRUE),
                                   df_rows$partner_call,
                                   df_rows$reciprocal)
colnames(df_cells) <- ls_header

# --- Write TSV -----------------------------------------------
data.table::fwrite(df_cells, file=opt$tsv, sep="\t", quote=FALSE)

# --- Write Markdown ------------------------------------------
# column widths and alignment (numbers are right-aligned)
ls_width <- sapply(seq_along(ls_header), function(i) max(nchar(ls_header[i]), nchar(df_cells[[i]])))
ls_right <- seq_along(ls_header) %in% c(3, 4, 5, 9, 10, 11, 12, 13, 14, 17)

# function: format one row of the table
f_md_row <- function(ls_cell) {
  return(paste0("| ", paste(f_pad(ls_cell, ls_width, ls_right), collapse=" | "), " |"))
}

ls_separator <- ifelse(ls_right, paste0(strrep("-", ls_width - 1), ":"), paste0(":", strrep("-", ls_width - 1)))

ls_md <- c(f_md_row(ls_header),
           paste0("| ", paste(ls_separator, collapse=" | "), " |"))
if (nrow(df_cells) > 0) {
  ls_md <- c(ls_md, apply(df_cells, 1, f_md_row))
}

writeLines(ls_md, opt$md, useBytes=TRUE)

# --- Summary -------------------------------------------------
ls_calls <- table(factor(df_rows$call, levels=ls_call_levels))
message(nrow(df_rows), " contigs scored; calls: {", paste0("'", names(ls_calls), "': ", ls_calls, collapse=", "), "}; ",
        "reciprocal partner pairs: ", sum(df_rows$reciprocal == "yes") %/% 2)

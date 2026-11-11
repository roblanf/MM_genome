#!/usr/bin/env Rscript

# ============================================================
#  Assign every contig of hap1 and hap2 to one of two parental genomes
#
#  Rules, in priority order:
#    chromosome : contig listed in --chrom (chromosome-scale contigs) -> parent given there
#    group      : fragment whose best partner in the other haplotype is a listed chromosome contig
#                 covering >= --min-cov of it -> the OPPOSITE parent of that partner (the fragment
#                 belongs to its own haplotype's copy of that chromosome, which has a gap there)
#    intrinsic  : ungrouped fragment -> its own evidence: delta ID call from dv_scoring.tsv and/or
#                 hapmer majority >= --hapmer-majority from hapmers.count; if both exist they must agree
#    unassigned : nothing confident
#
#  Also reports, for group fragments, whether their own evidence is consistent:
#    consistent      : intrinsic call == group parent
#    identical-filler: >= 99% identical to the partner chromosome (both parents carry it); resembling
#                      the partner's parent is expected and not a conflict
#    CHECK           : < 99% identical to the partner and intrinsic call == partner's parent -> possible
#                      misphased fragment; verify with read depth (half depth = other parent's allele)
#
#  Usage: Rscript split_parents.R --pairs both.pairs.tsv --chrom parent_of_chromosome.tsv \
#             --fai hap1.purged.fa.fai,hap2.purged.fa.fai [--dv dv_scoring.tsv] [--hapmers hapmers.count] \
#             [--min-cov 0.5] [--hapmer-majority 0.9] --out parent_of_contig.tsv
#  Writes <out>, plus virginea.txt / decipiens.txt / unassigned.txt next to it.
# ============================================================

# --- Load libraries and functions ---------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

# compare and sort contig names byte-wise (same as Python)
invisible(Sys.setlocale("LC_COLLATE", "C"))

# disable scientific notation
options(scipen=999)

# the opposite parent
ls_opposite <- c(virginea="decipiens", decipiens="virginea")

# function: read the parent of every chromosome-scale contig (TSV: contig <tab> virginea|decipiens)
f_read_chrom <- function(fn_chrom) {
  ls_lines <- readLines(fn_chrom)
  ls_lines <- ls_lines[trimws(ls_lines) != "" & !startsWith(ls_lines, "#")]

  ls_fields <- strsplit(ls_lines, "[[:space:]]+")
  ls_contig <- sapply(ls_fields, `[`, 1)
  ls_parent <- sapply(ls_fields, `[`, 2)

  # stop the code if a parent is invalid
  ls_invalid <- which(!ls_parent %in% names(ls_opposite))
  if (length(ls_invalid) > 0) {
    stop(paste("parent must be virginea or decipiens:", ls_lines[ls_invalid[1]]))
  }

  return(setNames(ls_parent, ls_contig))
}

# function: contig lengths from a samtools faidx .fai (or a FASTA file)
f_read_lengths <- function(fn_input) {
  if (endsWith(fn_input, ".fai")) {
    df_fai <- data.table::fread(fn_input, sep="\t", header=FALSE, select=1:2, col.names=c("contig", "length"))
    return(setNames(df_fai$length, df_fai$contig))
  }

  # FASTA: count the bases of every record
  ls_lines <- readLines(fn_input)
  ls_is_header <- startsWith(ls_lines, ">")

  ls_contig <- sapply(strsplit(sub("^>", "", ls_lines[ls_is_header]), "[[:space:]]+"), `[`, 1)
  ls_record <- cumsum(ls_is_header)[!ls_is_header]
  ls_length <- tapply(nchar(trimws(ls_lines[!ls_is_header])), ls_record, sum)

  ls_out <- setNames(rep(0, length(ls_contig)), ls_contig)
  ls_out[as.integer(names(ls_length))] <- ls_length

  return(ls_out)
}

# --- Argument parsing ----------------------------------------
option_list <- list(
  make_option(c("--pairs"), type="character", default=NULL,
              help="pairs table from pair_coverage.R [required]", metavar="FILE"),
  make_option(c("--chrom"), type="character", default=NULL,
              help="TSV: contig <tab> virginea|decipiens [required]", metavar="FILE"),
  make_option(c("--fai"), type="character", default=NULL,
              help="comma-separated samtools faidx .fai of hap1 and hap2 (or the FASTA files) [required]", metavar="FILES"),
  make_option(c("--dv"), type="character", default=NULL,
              help="TSV from dv_scoring.R", metavar="FILE"),
  make_option(c("--hapmers"), type="character", default=NULL,
              help="Merqury hapmers.count table", metavar="FILE"),
  make_option(c("--min-cov"), dest="min_cov", type="double", default=0.5,
              help="min coverage of a fragment by its chromosome partner [default %default]"),
  make_option(c("--hapmer-majority"), dest="hapmer_majority", type="double", default=0.9,
              help="min fraction of hapmers from one parent for a hapmer call [default %default]"),
  make_option(c("--min-hapmers"), dest="min_hapmers", type="integer", default=1000,
              help="min number of hapmers for a hapmer call [default %default]"),
  make_option(c("--out"), type="character", default="parent_of_contig.tsv",
              help="output table [default %default]", metavar="FILE")
)

# parse the arguments
opt <- parse_args(OptionParser(option_list=option_list))

# stop the code if any required input is missing
if (is.null(opt$pairs) || is.null(opt$chrom) || is.null(opt$fai)) {
  stop("--pairs, --chrom and --fai are required.\n
        Usage: Rscript split_parents.R --pairs both.pairs.tsv --chrom parent_of_chromosome.tsv --fai hap1.fa.fai,hap2.fa.fai --out parent_of_contig.tsv")
}

ls_fai <- unlist(strsplit(opt$fai, split=","))
for (fn in c(opt$pairs, opt$chrom, ls_fai, opt$dv, opt$hapmers)) {
  if (!file.exists(fn)) {
    stop(paste("File not found:", fn))
  }
}

# --- Read inputs ---------------------------------------------
# parent of every chromosome-scale contig
ls_parent_chrom <- f_read_chrom(opt$chrom)

# contig lengths
ls_length <- c()
for (fn in ls_fai) {
  ls_length <- c(ls_length, f_read_lengths(fn))
}
ls_length <- ls_length[!duplicated(names(ls_length), fromLast=TRUE)]

# best partner of every contig from the pairs table (both directions)
df_pairs <- data.table::fread(opt$pairs, sep="\t", colClasses=list(character=c("h1", "h2")))
df_best <- rbind(data.table::data.table(contig=df_pairs$h1, partner=df_pairs$h2, cov=df_pairs$h1_cov, identity=df_pairs$identity),
                 data.table::data.table(contig=df_pairs$h2, partner=df_pairs$h1, cov=df_pairs$h2_cov, identity=df_pairs$identity))
df_best <- df_best[, .SD[which.max(cov)], by=contig]

# intrinsic evidence: delta ID call
ls_dv_call <- c()
if (!is.null(opt$dv)) {
  df_dv <- data.table::fread(opt$dv, sep="\t", colClasses="character")
  df_dv <- df_dv[Call %in% names(ls_opposite)]
  ls_dv_call <- setNames(df_dv$Call, df_dv$Contig)
}

# intrinsic evidence: hapmer majority
ls_hap_call <- c()
if (!is.null(opt$hapmers)) {
  df_hapmers <- data.table::fread(opt$hapmers, sep="\t")
  col_vir <- grep("virginea", colnames(df_hapmers), value=TRUE)[1]
  col_dec <- grep("decipiens", colnames(df_hapmers), value=TRUE)[1]

  n_vir <- as.numeric(df_hapmers[[col_vir]])
  n_dec <- as.numeric(df_hapmers[[col_dec]])
  n_total <- n_vir + n_dec

  ls_call <- ifelse(n_total >= opt$min_hapmers & n_vir / n_total >= opt$hapmer_majority, "virginea",
             ifelse(n_total >= opt$min_hapmers & n_dec / n_total >= opt$hapmer_majority, "decipiens", NA))
  ls_hap_call <- setNames(ls_call, df_hapmers$Contig)
  ls_hap_call <- ls_hap_call[!is.na(ls_hap_call)]
}

# --- Assign every contig -------------------------------------
# sort contigs by haplotype (first two characters of the name) and length
df_contig <- data.table::data.table(contig=names(ls_length), length=as.numeric(ls_length))
df_contig$hap <- ifelse(startsWith(df_contig$contig, "h1"), "hap1", "hap2")
df_contig <- df_contig[order(substr(contig, 1, 2), -length)]

# best partner, coverage and identity
idx_best <- match(df_contig$contig, df_best$contig)
df_contig$partner          <- ifelse(is.na(idx_best), "-", df_best$partner[idx_best])
df_contig$partner_cov      <- ifelse(is.na(idx_best), 0, df_best$cov[idx_best])
df_contig$partner_identity <- ifelse(is.na(idx_best), 0, df_best$identity[idx_best])

# intrinsic evidence
df_contig$dv_call     <- ifelse(df_contig$contig %in% names(ls_dv_call), ls_dv_call[df_contig$contig], NA)
df_contig$hapmer_call <- ifelse(df_contig$contig %in% names(ls_hap_call), ls_hap_call[df_contig$contig], NA)

# intrinsic call: the available evidence must agree
df_contig$intrinsic <- ifelse(is.na(df_contig$dv_call), df_contig$hapmer_call,
                       ifelse(is.na(df_contig$hapmer_call) | df_contig$hapmer_call == df_contig$dv_call, df_contig$dv_call, NA))

# apply the rules in priority order
df_contig$rule <- "unassigned"
df_contig$parent <- "unassigned"
df_contig$check <- "-"

for (i in seq_len(nrow(df_contig))) {
  contig  <- df_contig$contig[i]
  partner <- df_contig$partner[i]
  intr    <- df_contig$intrinsic[i]

  if (contig %in% names(ls_parent_chrom)) {
    df_contig$rule[i]   <- "chromosome"
    df_contig$parent[i] <- ls_parent_chrom[contig]

  } else if (partner %in% names(ls_parent_chrom) && df_contig$partner_cov[i] >= opt$min_cov) {
    df_contig$rule[i]   <- "group"
    df_contig$parent[i] <- ls_opposite[ls_parent_chrom[partner]]

    # is the fragment's own evidence consistent with its group?
    if (is.na(intr)) {
      df_contig$check[i] <- "no intrinsic call"
    } else if (intr == df_contig$parent[i]) {
      df_contig$check[i] <- "consistent"
    } else if (df_contig$partner_identity[i] >= 0.99) {
      df_contig$check[i] <- "identical-filler"
    } else {
      df_contig$check[i] <- "CHECK"
    }

  } else if (!is.na(intr)) {
    df_contig$rule[i]   <- "intrinsic"
    df_contig$parent[i] <- intr
  }
}

# --- Output --------------------------------------------------
df_out <- data.table::data.table(contig           = df_contig$contig,
                                 hap              = df_contig$hap,
                                 length           = df_contig$length,
                                 rule             = df_contig$rule,
                                 parent           = df_contig$parent,
                                 partner          = df_contig$partner,
                                 partner_cov      = sprintf("%.2f", df_contig$partner_cov),
                                 partner_identity = sprintf("%.3f", df_contig$partner_identity),
                                 dv_call          = ifelse(is.na(df_contig$dv_call), "-", df_contig$dv_call),
                                 hapmer_call      = ifelse(is.na(df_contig$hapmer_call), "-", df_contig$hapmer_call),
                                 check            = df_contig$check)

data.table::fwrite(df_out, file=opt$out, sep="\t", quote=FALSE)

# list of contigs per parent, next to the output table
outdir <- dirname(normalizePath(opt$out))
for (parent in c("virginea", "decipiens", "unassigned")) {
  writeLines(df_contig$contig[df_contig$parent == parent], file.path(outdir, paste0(parent, ".txt")))
}

# --- Summary -------------------------------------------------
df_summary <- df_contig[, .(n=.N, length=sum(length)), by=.(parent, rule)]
df_summary <- df_summary[order(parent, rule)]

cat(sprintf("%-11s %-11s n=%4d  %7.2f Mb\n", df_summary$parent, df_summary$rule, df_summary$n, df_summary$length / 1e6), sep="")
cat(sprintf("fragments flagged CHECK: %d\n", sum(df_contig$check == "CHECK")))

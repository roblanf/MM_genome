#!/usr/bin/env Rscript

# ============================================================
#  Summarise BUSCO (compleasm) and k-mer (Merqury) results
#  across the hifiasm -s sweep
#
#  Expected layout (run from 02_assembly/homcov55/):
#    compleasm/<label>/summary.txt
#    merqury/<label>/<label>_merq.qv
#    merqury/<label>/<label>_merq.completeness.stats
#    merqury/<label>/<label>_merq.spectra-cn.hist
#  where <label> looks like s035_hap1, s035_hap2, s035_union, ...
#
#  Prints three Markdown tables: BUSCO completeness, Merqury QV / k-mer
#  completeness, and a copy-number breakdown of the het-peak and hom-peak
#  k-mers from spectra-cn.hist.
#
#  Usage: Rscript scripts/summarise_s_runs.R [--compleasm compleasm] [--merqury merqury] [--het 28] [--hom 55]
# ============================================================

# --- Load libraries and functions ---------------------------
suppressPackageStartupMessages(library(optparse))

# disable scientific notation
options(scipen=999)

# order of the assemblies within one -s value
ls_asm_order <- c(hap1=0, hap2=1, p_ctg=2, union=3)

# function: sort labels as s030_hap1, s030_hap2, s030_union, s035_hap1, ...; unknown labels go last
f_sort_labels <- function(ls_label) {
  ls_match <- regmatches(ls_label, regexec("^s([0-9]+)_(\\w+)$", ls_label, perl=TRUE))

  ls_s_val <- sapply(ls_match, function(m) if (length(m) == 3) as.numeric(m[2]) / 100 else 9)
  ls_asm   <- sapply(ls_match, function(m) {
    if (length(m) == 3 && m[3] %in% names(ls_asm_order)) ls_asm_order[m[3]] else if (length(m) == 3) 50 else 99
  })

  return(ls_label[order(ls_s_val, ls_asm, ls_label)])
}

# function: pad a string to a given width (character-aware)
f_pad <- function(x, width, right=FALSE) {
  padding <- strrep(" ", pmax(0, width - nchar(x)))

  return(ifelse(right, paste0(padding, x), paste0(x, padding)))
}

# function: print a Markdown table; align is a vector of "<" (left) or ">" (right) per column
f_print_table <- function(ls_header, ls_rows, ls_align) {
  ls_width <- sapply(seq_along(ls_header), function(i) max(nchar(ls_header[i]), nchar(sapply(ls_rows, `[`, i))))
  ls_right <- ls_align == ">"

  # function: format one row of the table
  f_row <- function(ls_cell) {
    return(paste0("| ", paste(f_pad(ls_cell, ls_width, ls_right), collapse=" | "), " |"))
  }

  ls_separator <- ifelse(ls_right, paste0(strrep("-", ls_width - 1), ":"), strrep("-", ls_width))

  cat(f_row(ls_header), "\n", sep="")
  cat("| ", paste(ls_separator, collapse=" | "), " |\n", sep="")
  for (row in ls_rows) {
    cat(f_row(row), "\n", sep="")
  }
  cat("\n")
}

# function: thousands separator
f_comma <- function(x) {
  return(formatC(as.numeric(x), format="d", big.mark=","))
}

# --- compleasm -----------------------------------------------
# function: S/D/F/I/M as "xx.xx% (N)" and the total N from a compleasm summary
f_parse_compleasm <- function(fn_summary) {
  ls_stats <- c()
  total <- "0"

  for (line in readLines(fn_summary)) {
    m <- regmatches(line, regexec("([SDFIM]):([0-9.]+)%,\\s+([0-9]+)", line, perl=TRUE))[[1]]
    if (length(m) == 4) {
      ls_stats[m[2]] <- paste0(m[3], "% (", m[4], ")")
    }

    n <- regmatches(trimws(line), regexec("^[Nn]:([0-9]+)", trimws(line), perl=TRUE))[[1]]
    if (length(n) == 2) {
      total <- n[2]
    }
  }

  return(list(stats=ls_stats, total=total))
}

# function: print the BUSCO table
f_busco_table <- function(dir_compleasm) {
  ls_files <- Sys.glob(file.path(dir_compleasm, "*", "summary.txt"))
  if (length(ls_files) == 0) {
    cat("No compleasm summaries found under ", dir_compleasm, "/*/summary.txt\n\n", sep="")
    return(invisible(NULL))
  }

  # sort by label
  ls_labels <- basename(dirname(ls_files))
  ls_files <- ls_files[match(f_sort_labels(ls_labels), ls_labels)]

  ls_rows <- list()
  for (fn in ls_files) {
    label <- basename(dirname(fn))
    summary <- f_parse_compleasm(fn)

    ls_values <- sapply(c("S", "D", "F", "I", "M"), function(k) {
      if (k %in% names(summary$stats)) summary$stats[k] else "0.00% (0)"
    })
    ls_rows[[length(ls_rows) + 1]] <- c(label, ls_values, summary$total)
  }

  cat("### BUSCO completeness (compleasm, eudicotyledons_odb12)\n\n")
  f_print_table(c("Dataset", "S (%, (N))", "D (%, (N))", "F (%, (N))", "I (%, (N))", "M (%, (N))", "Total"),
                ls_rows, c("<", ">", ">", ">", ">", ">", ">"))
}

# --- merqury -------------------------------------------------
# function: read a whitespace-delimited Merqury table (first column = assembly name) with >= 5 fields
f_read_merqury <- function(fn_input) {
  ls_rows <- strsplit(trimws(readLines(fn_input)), "[[:space:]]+")
  ls_rows <- ls_rows[sapply(ls_rows, length) >= 5]

  return(setNames(ls_rows, sapply(ls_rows, `[`, 1)))
}

# function: print the QV and k-mer completeness table; handles multi-line files (two-assembly mode)
f_merqury_table <- function(dir_merqury) {
  ls_files <- Sys.glob(file.path(dir_merqury, "*", "*_merq.qv"))
  if (length(ls_files) == 0) {
    cat("No Merqury .qv files found under ", dir_merqury, "/*/*_merq.qv\n\n", sep="")
    return(invisible(NULL))
  }

  # sort by label
  ls_labels <- basename(dirname(ls_files))
  ls_files <- ls_files[match(f_sort_labels(ls_labels), ls_labels)]

  ls_rows <- list()
  for (fn_qv in ls_files) {
    label <- basename(dirname(fn_qv))
    fn_comp <- sub("\\.qv$", ".completeness.stats", fn_qv)

    # <asm>  asm_only_kmers  total_kmers  QV  error_rate
    ls_qv <- f_read_merqury(fn_qv)

    # <asm>  all  found  total  completeness%
    ls_comp <- list()
    if (file.exists(fn_comp)) {
      ls_comp <- f_read_merqury(fn_comp)
    }

    for (asm in names(ls_qv)) {
      p <- ls_qv[[asm]]
      c <- ls_comp[[asm]]
      if (is.null(c)) {
        c <- ls_comp[[tolower(asm)]]
      }

      ls_rows[[length(ls_rows) + 1]] <- c(
        if (length(ls_qv) == 1) label else paste0(label, " (", asm, ")"),
        p[4],                                                             # QV
        sprintf("%.2e", as.numeric(p[5])),                                # error rate
        if (is.null(c)) "NA" else c[5],                                   # completeness %
        if (is.null(c)) "NA" else paste(f_comma(c[3]), "/", f_comma(c[4])),
        f_comma(p[2])                                                     # assembly-only (erroneous) k-mers
      )
    }
  }

  cat("### k-mer completeness and QV (Merqury, k=21)\n\n")
  f_print_table(c("Dataset", "QV", "Error rate", "Completeness (%)", "Solid k-mers found / total", "Asm-only k-mers"),
                ls_rows, c("<", ">", ">", ">", ">", ">"))
}

# function: where do the het-peak and hom-peak read k-mers land in the assembly?
#           from <label>_merq.spectra-cn.hist (columns: Copies, kmer_multiplicity, Count)
#           windows: het = [0.5*het, 1.5*het), hom = [0.8*hom, 1.25*hom)
f_copy_number_table <- function(dir_merqury, het_peak, hom_peak) {
  # two-assembly runs write <label>_merq.spectra-cn.hist (the union); single-assembly
  # runs only write <label>_merq.<asm>.spectra-cn.hist, so fall back to that
  ls_dirs <- Sys.glob(file.path(dir_merqury, "*"))
  ls_dirs <- ls_dirs[match(f_sort_labels(basename(ls_dirs)), basename(ls_dirs))]

  ls_files <- c()
  for (dir in ls_dirs) {
    label <- basename(dir)
    fn_union <- file.path(dir, paste0(label, "_merq.spectra-cn.hist"))

    if (file.exists(fn_union)) {
      ls_files <- c(ls_files, fn_union)
    } else {
      ls_files <- c(ls_files, sort(Sys.glob(file.path(dir, paste0(label, "_merq.*.spectra-cn.hist")))))
    }
  }

  if (length(ls_files) == 0) {
    cat("No Merqury spectra-cn.hist files found under ", dir_merqury, "\n\n", sep="")
    return(invisible(NULL))
  }

  het_lo <- 0.5 * het_peak
  het_hi <- 1.5 * het_peak
  hom_lo <- 0.8 * hom_peak
  hom_hi <- 1.25 * hom_peak

  # function: copy-number bucket of a k-mer
  f_bucket <- function(ls_copies) {
    return(ifelse(ls_copies == "read-only", "0", ifelse(ls_copies %in% c("1", "2"), ls_copies, "3+")))
  }

  # function: percentage of a bucket within a window
  f_pct <- function(ls_count, ls_bucket, bucket, total) {
    return(sprintf("%5.1f", 100 * sum(ls_count[ls_bucket %in% bucket]) / total))
  }

  ls_rows <- list()
  for (fn in ls_files) {
    label <- basename(dirname(fn))

    ls_fields <- strsplit(readLines(fn), "\t", fixed=TRUE)
    ls_fields <- ls_fields[sapply(ls_fields, length) >= 3]
    ls_fields <- ls_fields[sapply(ls_fields, `[`, 1) != "Copies"]

    ls_bucket <- f_bucket(sapply(ls_fields, `[`, 1))
    ls_mult   <- as.numeric(sapply(ls_fields, `[`, 2))
    ls_count  <- floor(as.numeric(sapply(ls_fields, `[`, 3)))

    is_het <- ls_mult >= het_lo & ls_mult < het_hi
    is_hom <- !is_het & ls_mult >= hom_lo & ls_mult < hom_hi

    het_total <- max(sum(ls_count[is_het]), 1)
    hom_total <- max(sum(ls_count[is_hom]), 1)

    ls_rows[[length(ls_rows) + 1]] <- c(
      label,
      f_pct(ls_count[is_het], ls_bucket[is_het], "0", het_total),
      f_pct(ls_count[is_het], ls_bucket[is_het], "1", het_total),
      f_pct(ls_count[is_het], ls_bucket[is_het], c("2", "3+"), het_total),
      f_pct(ls_count[is_hom], ls_bucket[is_hom], "0", hom_total),
      f_pct(ls_count[is_hom], ls_bucket[is_hom], "1", hom_total),
      f_pct(ls_count[is_hom], ls_bucket[is_hom], "2", hom_total),
      f_pct(ls_count[is_hom], ls_bucket[is_hom], "3+", hom_total)
    )
  }

  cat(sprintf("### Copy number of het-peak (%.0f-%.0fx) and hom-peak (%.0f-%.0fx) k-mers in the assembly (%%)\n\n",
              het_lo, het_hi, hom_lo, hom_hi))
  f_print_table(c("Dataset", "het: 0x", "het: 1x", "het: >=2x", "hom: 0x", "hom: 1x", "hom: 2x", "hom: >=3x"),
                ls_rows, c("<", ">", ">", ">", ">", ">", ">", ">"))
  cat("Expected for one haplotype: het ~50% 0x / ~50% 1x, hom ~100% 1x. ",
      "Expected for hap1+hap2 together: het ~100% 1x, hom ~100% 2x; ",
      "het k-mers at >=2x or hom k-mers at >=3x are duplicated regions.\n\n", sep="")
}

# --- Argument parsing ----------------------------------------
option_list <- list(
  make_option(c("--compleasm"), type="character", default="compleasm",
              help="directory holding <label>/summary.txt [default %default]", metavar="DIR"),
  make_option(c("--merqury"), type="character", default="merqury",
              help="directory holding <label>/<label>_merq.* [default %default]", metavar="DIR"),
  make_option(c("--het"), type="double", default=28,
              help="heterozygous k-mer peak (GenomeScope kcov) [default %default]"),
  make_option(c("--hom"), type="double", default=55,
              help="homozygous k-mer peak (2 x kcov) [default %default]")
)

# parse the arguments
opt <- parse_args(OptionParser(option_list=option_list))

# --- Print the tables ----------------------------------------
f_busco_table(opt$compleasm)
f_merqury_table(opt$merqury)
f_copy_number_table(opt$merqury, opt$het, opt$hom)

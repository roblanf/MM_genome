#!/usr/bin/env Rscript

# ============================================================
#  Draw a circular gene map of an oatk organelle assembly
#  (mitochondrion or plastid)
#
#  Inputs are the two files written by oatk / pathfinder:
#    prefix.mito.ctg.fasta + prefix.mito.ctg.bed     (mitochondrion)
#    prefix.pltd.ctg.fasta + prefix.pltd.ctg.bed     (plastid)
#  FASTA headers carry length=, circular=, path=; the BED is BED6
#  (seq, start, end, gene, score, strand) with a '#' header line.
#
#  Usage: Rscript oatk_circular_map.R -f x.mito.ctg.fasta -b x.mito.ctg.bed -o mito_map.png
#         Rscript oatk_circular_map.R -f x.pltd.ctg.fasta -b x.pltd.ctg.bed -o pltd_map.svg
#         Rscript oatk_circular_map.R -f x.mito.ctg.fasta -b x.mito.ctg.bed --contig ctg000001c --genbank mito.gb
#
#  Requires: optparse, data.table, circlize
# ============================================================

# --- Load libraries and functions ---------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(circlize))

# disable scientific notation
options(scipen=999)

# --- Gene categories (OatkDB embryophyta gene names) and colours
df_mito_categories <- data.table::data.table(
  label   = c("NADH dehydrogenase (complex I)", "Succinate dehydrogenase (II)", "Ubiquinol cytochrome c reductase (III)",
              "Cytochrome c oxidase (IV)", "ATP synthase (V)", "Cytochrome c biogenesis", "Ribosomal proteins",
              "Maturase / other proteins", "rRNA", "tRNA"),
  colour  = c("#ffd166", "#ef8354", "#d62828", "#c77dff", "#06a77d", "#a0522d", "#118ab2", "#6c757d", "#ff006e", "#2b2d42"),
  pattern = c("^nad", "^sdh", "^cob", "^cox", "^atp", "^ccm", "^rp[sl]", "^(matR|mttB|rtl|orf)", "^rrn", "^trn")
)

df_pltd_categories <- data.table::data.table(
  label   = c("Photosystem I", "Photosystem II", "Cytochrome b/f complex", "ATP synthase", "NADH dehydrogenase",
              "RubisCO large subunit", "RNA polymerase", "Ribosomal proteins", "Chlorophyll biosynthesis", "Other genes",
              "Hypothetical ORFs (ycf)", "rRNA", "tRNA"),
  colour  = c("#1b5e20", "#4caf50", "#a5d6a7", "#00acc1", "#ffd166", "#5d4037", "#d62828", "#118ab2", "#8e24aa",
              "#6c757d", "#adb5bd", "#ff006e", "#2b2d42"),
  pattern = c("^psa", "^psb", "^pet", "^atp", "^ndh", "^rbcL", "^rpo", "^rp[sl]", "^chl",
              "^(matK|clpP|accD|cemA|ccsA|infA|tufA|cys)", "^ycf", "^rrn", "^trn")
)

ls_other <- c(label="Other", colour="#cccccc")
pattern_pltd_markers <- "^(psa|psb|pet|ndh|rbcL|rpo|ycf|matK|clpP|accD|cemA|ccsA|rrn16|rrn23|rrn4\\.5)"

# function: category (label and colour) of a gene; the first matching pattern wins
f_category <- function(gene, df_categories) {
  for (i in seq_len(nrow(df_categories))) {
    if (grepl(df_categories$pattern[i], gene, ignore.case=TRUE, perl=TRUE)) {
      return(c(label=df_categories$label[i], colour=df_categories$colour[i]))
    }
  }

  return(ls_other)
}

# --- Parsers ---
# function: read a FASTA file into a named list of records (seq, length, circular, path, header)
f_read_fasta <- function(fn_fasta) {
  ls_lines <- readLines(fn_fasta)
  ls_is_header <- startsWith(ls_lines, ">")

  ls_header <- sub("^>", "", ls_lines[ls_is_header])
  ls_name <- sapply(strsplit(ls_header, "[[:space:]]+"), `[`, 1)

  # concatenate the sequence lines of every record
  ls_record <- cumsum(ls_is_header)[!ls_is_header]
  ls_chunk <- trimws(ls_lines[!ls_is_header])
  ls_seq <- rep("", length(ls_name))
  if (length(ls_chunk) > 0) {
    ls_concat <- tapply(ls_chunk, ls_record, paste, collapse="")
    ls_seq[as.integer(names(ls_concat))] <- ls_concat
  }
  ls_seq <- toupper(ls_seq)

  ls_contigs <- list()
  for (i in seq_along(ls_name)) {
    # header tags (key=value)
    ls_tags <- regmatches(ls_header[i], gregexpr("(\\w+)=(\\S+)", ls_header[i], perl=TRUE))[[1]]
    ls_tags <- setNames(sub("^\\w+=", "", ls_tags), sub("=.*$", "", ls_tags))

    length <- if ("length" %in% names(ls_tags)) as.integer(ls_tags["length"]) else nchar(ls_seq[i])

    # oatk naming: ctg000001c / ctg000001l
    circular <- if ("circular" %in% names(ls_tags)) tolower(ls_tags["circular"]) == "true" else endsWith(ls_name[i], "c")

    ls_contigs[[ls_name[i]]] <- list(seq=ls_seq[i], length=length, circular=unname(circular),
                                     path=if ("path" %in% names(ls_tags)) unname(ls_tags["path"]) else "",
                                     header=ls_header[i])
  }

  return(ls_contigs)
}

# function: read a BED6 file into a data.table (seq, start, end, gene, score, strand)
f_read_bed <- function(fn_bed, min_score=0) {
  ls_lines <- readLines(fn_bed)
  ls_lines <- ls_lines[trimws(ls_lines) != "" & !startsWith(ls_lines, "#")]

  # split by tab, or by whitespace when there are fewer than 6 columns
  ls_fields <- strsplit(ls_lines, "\t", fixed=TRUE)
  idx <- vapply(ls_fields, length, integer(1)) < 6
  ls_fields[idx] <- strsplit(ls_lines[idx], "[[:space:]]+")
  ls_fields <- ls_fields[vapply(ls_fields, length, integer(1)) >= 6]

  df_bed <- data.table::data.table(seq    = vapply(ls_fields, `[`, character(1), 1),
                                   start  = as.integer(vapply(ls_fields, `[`, character(1), 2)),
                                   end    = as.integer(vapply(ls_fields, `[`, character(1), 3)),
                                   gene   = vapply(ls_fields, `[`, character(1), 4),
                                   score  = suppressWarnings(as.numeric(vapply(ls_fields, `[`, character(1), 5))),
                                   strand = vapply(ls_fields, `[`, character(1), 6))
  df_bed$score <- ifelse(is.na(df_bed$score), 0, df_bed$score)
  df_bed <- df_bed[score >= min_score]

  # make sure start <= end
  idx <- df_bed$start > df_bed$end
  df_bed[idx, c("start", "end") := .(end, start)]

  return(df_bed[order(seq, start, end, gene, score, strand)])
}

# --- Helpers ---
# function: pick a 1/2/5 x 10^k tick interval giving roughly target_ticks ticks
f_nice_interval <- function(length, target_ticks=20) {
  raw <- length / target_ticks
  exp <- 10^floor(log10(raw))

  for (m in c(1, 2, 5, 10)) {
    if (m * exp >= raw) {
      return(m * exp)
    }
  }

  return(10 * exp)
}

# function: sliding-window GC fraction; returns window centres (x) and GC values (y)
f_gc_content <- function(seq, window) {
  n <- nchar(seq)
  if (n == 0) {
    return(list(x=numeric(0), y=numeric(0)))
  }

  ls_start <- seq(0, n - 1, by=window)
  ls_end <- pmin(ls_start + window, n)
  ls_chunk <- substring(seq, ls_start + 1, ls_end)

  ls_gc <- nchar(gsub("[^GC]", "", ls_chunk))
  ls_base <- nchar(gsub("[^ACGT]", "", ls_chunk))

  idx <- ls_base > 0
  return(list(x=((ls_start + ls_end) / 2)[idx], y=(ls_gc / ls_base)[idx]))
}

# function: reverse complement
f_revcomp <- function(seq) {
  return(intToUtf8(rev(utf8ToInt(chartr("ACGTNacgtn", "TGCANtgcan", seq)))))
}

# function: locate the largest pair of inverted repeats (plastid IRA/IRB)
#           returns list(a=c(start, end), b=c(start, end)) as 0-based half-open intervals on the
#           original sequence (an interval with start > end wraps around the origin),
#           or NULL if no inverted repeat >= min_len is found
f_find_inverted_repeat <- function(seq, circular=TRUE, k=31, min_len=1000, max_gap=100) {
  L <- nchar(seq)
  if (L < 2 * min_len) {
    return(NULL)
  }

  # doubling lets a repeat cross the origin
  s <- if (circular) paste0(seq, seq) else seq
  n <- nchar(s)
  rc <- f_revcomp(s)

  # every k-mer of the forward and reverse-complement sequence (0-based positions)
  df_fwd <- data.table::data.table(kmer=substring(s, 1:(n - k + 1), k:n), p=0:(n - k))
  df_rc  <- data.table::data.table(kmer=substring(rc, 1:(n - k + 1), k:n), i=0:(n - k))
  df_hit <- merge(df_fwd, df_rc, by="kmer", allow.cartesian=TRUE)

  # drop palindromic self-matches
  df_hit <- df_hit[abs(p - (n - i - k)) >= k]
  if (nrow(df_hit) == 0) {
    return(NULL)
  }

  # diagonals, ordered by first appearance (same as the Python script)
  df_hit$d <- df_hit$p - df_hit$i
  df_hit[, first_i := min(i), by=d]
  df_hit <- df_hit[order(first_i, d, p)]

  # runs of k-mers on the same diagonal, broken by gaps > max_gap
  df_hit[, new_run := c(TRUE, diff(p) > max_gap), by=d]
  df_hit$run <- cumsum(df_hit$new_run)

  df_run <- df_hit[, .(start=min(p), end=max(p) + k), by=.(run, d)]
  df_run$length <- df_run$end - df_run$start
  df_run <- df_run[length >= min_len]
  if (nrow(df_run) == 0) {
    return(NULL)
  }

  best <- df_run[which.max(length)]
  if (best$length > L / 2 || (circular && best$length >= L - 1)) {
    return(NULL)
  }

  # partner copy: i = p - d, forward position q = n - i - k
  a0 <- best$start
  a1 <- best$end
  b0 <- n - (a1 - k - best$d) - k
  b1 <- n - (a0 - best$d)

  if (circular) {
    a0 <- a0 %% L
    a1 <- a1 %% L
    b0 <- b0 %% L
    b1 <- b1 %% L
  }

  return(list(a=c(a0, a1), b=c(b0, b1)))
}

# function: length of an arc on a circular sequence
f_arc_len <- function(start, end, L) {
  if (start == end) {
    return(0)
  }

  return((end - start) %% L)
}

# function: format a length in bp / kb / Mb
f_fmt_bp <- function(x) {
  if (x >= 1000000) {
    return(sprintf("%g Mb", x / 1e6))
  }
  if (x >= 1000) {
    return(sprintf("%g kb", x / 1000))
  }

  return(sprintf("%g bp", x))
}

# function: thousands separator
f_comma <- function(x) {
  return(formatC(x, format="d", big.mark=","))
}

# function: guess the organelle from the file name, then from the gene names
f_detect_organelle <- function(fn_fasta, df_bed) {
  name <- tolower(fn_fasta)
  if (grepl("pltd", name) || grepl("plastid", name) || grepl("chloro", name)) {
    return("pltd")
  }
  if (grepl("mito", name)) {
    return("mito")
  }

  if (any(grepl(pattern_pltd_markers, df_bed$gene, ignore.case=TRUE, perl=TRUE))) {
    return("pltd")
  }

  return("mito")
}

# --- Optional GenBank export (e.g. for OGDRAW, Geneious, SnapGene) ---
# function: write the contigs and their annotation as GenBank
f_write_genbank <- function(ls_contigs, df_bed, fn_out) {
  today <- Sys.Date()
  date <- sprintf("%02d-%s-%s", as.integer(format(today, "%d")), toupper(month.abb[as.integer(format(today, "%m"))]), format(today, "%Y"))

  ls_out <- c()
  for (name in names(ls_contigs)) {
    contig <- ls_contigs[[name]]
    L <- contig$length
    topo <- if (contig$circular) "circular" else "linear"

    ls_out <- c(ls_out,
                sprintf("LOCUS       %-23s %7d bp    DNA     %-8s PLN %s", name, L, topo, date),
                sprintf("DEFINITION  %s oatk organelle assembly, %s.", name, topo),
                sprintf("ACCESSION   %s", name),
                sprintf("VERSION     %s", name),
                "KEYWORDS    .", "SOURCE      .", "  ORGANISM  .", "            .",
                "FEATURES             Location/Qualifiers",
                sprintf("     source          1..%d", L))
    if (contig$path != "") {
      ls_out <- c(ls_out, sprintf("                     /note=\"oatk path %s\"", contig$path))
    }

    df_feat <- df_bed[seq == name]
    for (i in seq_len(nrow(df_feat))) {
      # oatk writes nhmmscan's 1-based alignment coordinates into the BED
      # (not 0-based BED convention), so they are used as-is; ends can be +/-1 bp
      start <- max(1, df_feat$start[i])
      end <- min(L, df_feat$end[i])
      if (end < start) {
        next
      }

      gene <- df_feat$gene[i]
      loc <- if (df_feat$strand[i] == "+") sprintf("%d..%d", start, end) else sprintf("complement(%d..%d)", start, end)
      key <- if (startsWith(tolower(gene), "trn")) "tRNA" else if (startsWith(tolower(gene), "rrn")) "rRNA" else "CDS"

      ls_out <- c(ls_out,
                  sprintf("     gene            %s", loc),
                  sprintf("                     /gene=\"%s\"", gene),
                  sprintf("     %-16s%s", key, loc),
                  sprintf("                     /gene=\"%s\"", gene),
                  sprintf("                     /note=\"oatk hmmannot score %g\"", df_feat$score[i]))
    }

    # sequence in blocks of 60 bases (6 x 10)
    ls_out <- c(ls_out, "ORIGIN")
    seq <- tolower(contig$seq)
    for (i in seq_len(ceiling(nchar(seq) / 60)) * 60 - 59) {
      chunk <- substring(seq, i, min(i + 59, nchar(seq)))
      ls_block <- substring(chunk, seq(1, nchar(chunk), by=10), pmin(seq(10, nchar(chunk) + 9, by=10), nchar(chunk)))
      ls_out <- c(ls_out, sprintf("%9d %s", i, paste(ls_block, collapse=" ")))
    }
    ls_out <- c(ls_out, "//")
  }

  writeLines(ls_out, fn_out)
  message("wrote ", fn_out)
}

# --- Plot ---
# function: open the output graphics device based on the file extension
f_open_device <- function(fn_out, width, height, dpi) {
  ext <- tolower(tools::file_ext(fn_out))

  if (ext == "png") {
    png(fn_out, width=width, height=height, units="in", res=dpi)
  } else if (ext == "svg") {
    svg(fn_out, width=width, height=height)
  } else if (ext == "pdf") {
    pdf(fn_out, width=width, height=height)
  } else {
    stop(paste("unsupported output format (use .png, .svg or .pdf):", fn_out))
  }
}

# function: draw the circular map
f_draw <- function(ls_contigs, df_bed, fn_out, organelle="mito", title=NULL, label_size=7, no_trn=FALSE,
                   no_labels=FALSE, gc_window=1000, dpi=300, figsize=9, show_ir=TRUE) {
  df_categories <- if (organelle == "pltd") df_pltd_categories else df_mito_categories

  ls_names <- names(ls_contigs)
  single <- length(ls_names) == 1
  ls_length <- sapply(ls_contigs, `[[`, "length")
  total <- sum(ls_length)

  # --- genes, labels and categories of every contig ----------
  df_genes <- df_bed[seq %in% ls_names]
  df_genes$start <- pmax(0, pmin(df_genes$start, ls_length[df_genes$seq]))
  df_genes$end <- pmax(0, pmin(df_genes$end, ls_length[df_genes$seq]))
  df_genes <- df_genes[end > start]

  ls_category <- lapply(df_genes$gene, f_category, df_categories=df_categories)
  df_genes$label <- vapply(ls_category, `[[`, character(1), "label")
  df_genes$colour <- vapply(ls_category, `[[`, character(1), "colour")
  df_genes$deg <- 360 * (df_genes$end - df_genes$start) / total

  # categories used, in the order of the legend
  ls_used <- unique(df_genes[, .(label, colour)])
  ls_used <- ls_used[match(intersect(c(df_categories$label, ls_other["label"]), ls_used$label), ls_used$label)]

  # gene labels (tRNA labels can be hidden)
  df_labels <- df_genes
  if (no_labels) {
    df_labels <- df_labels[0]
  } else if (no_trn) {
    df_labels <- df_labels[!startsWith(tolower(gene), "trn")]
  }

  # --- GC content and inverted repeats of every contig -------
  ls_gc <- list()
  ls_ir <- list()
  for (name in ls_names) {
    contig <- ls_contigs[[name]]
    L <- contig$length

    if (nchar(contig$seq) > 0 && gc_window > 0 && nchar(contig$seq) > gc_window) {
      ls_gc[[name]] <- f_gc_content(contig$seq, gc_window)
    }

    # plastid quadripartite structure (LSC / IRB / SSC / IRA)
    if (organelle == "pltd" && show_ir && nchar(contig$seq) > 0 && L >= 20000) {
      ir <- f_find_inverted_repeat(contig$seq, circular=contig$circular)

      if (is.null(ir)) {
        message(name, "\tno inverted repeat >= 1 kb found")
        next
      }

      a <- ir$a
      b <- ir$b
      if (contig$circular) {
        # clockwise from the end of one copy you reach the start of the other:
        # the two single-copy regions are (a1 -> b0) and (b1 -> a0)
        sc1 <- c(a[2], b[1])
        sc2 <- c(b[2], a[1])
        if (f_arc_len(sc1[1], sc1[2], L) >= f_arc_len(sc2[1], sc2[2], L)) {
          lsc <- sc1; ssc <- sc2; irb <- b; ira <- a
        } else {
          lsc <- sc2; ssc <- sc1; irb <- a; ira <- b
        }
        df_regions <- data.table::data.table(name=c("LSC", "IRB", "SSC", "IRA"),
                                             start=c(lsc[1], irb[1], ssc[1], ira[1]),
                                             end=c(lsc[2], irb[2], ssc[2], ira[2]),
                                             colour=c("#e0e0e0", "#9e9e9e", "#e0e0e0", "#9e9e9e"))
      } else {
        if (a[1] > b[1]) {
          tmp <- a; a <- b; b <- tmp
        }
        lsc <- c(0, a[1]); ssc <- c(a[2], b[1]); irb <- a; ira <- b
        df_regions <- data.table::data.table(name=c("IR", "IR"), start=c(irb[1], ira[1]), end=c(irb[2], ira[2]),
                                             colour=c("#9e9e9e", "#9e9e9e"))
      }

      ls_ir[[name]] <- df_regions
      message(name, "\tinverted repeat: IRA ", ira[1] + 1, "-", ifelse(ira[2] == 0, L, ira[2]),
              " / IRB ", irb[1] + 1, "-", ifelse(irb[2] == 0, L, irb[2]),
              " (", f_comma(f_arc_len(irb[1], irb[2], L)), " bp), LSC ", f_comma(f_arc_len(lsc[1], lsc[2], L)),
              " bp, SSC ", f_comma(f_arc_len(ssc[1], ssc[2], L)), " bp")
    }
  }
  ir_found <- length(ls_ir) > 0

  # --- circular layout ---------------------------------------
  f_open_device(fn_out, width=figsize + 3, height=figsize, dpi=dpi)
  layout(matrix(1:2, nrow=1), widths=c(figsize, 3))
  par(mar=c(0.5, 0.5, 0.5, 0.5))

  circos.clear()
  circos.par(start.degree=90, gap.degree=ifelse(single, 0, 4), track.margin=c(0.005, 0.005),
             cell.padding=c(0, 0, 0, 0), points.overflow.warning=FALSE)
  circos.initialize(sectors=ls_names, xlim=cbind(0, ls_length))

  # --- gene labels: outer ring, italic for the - strand ------
  if (nrow(df_labels) > 0) {
    circos.labels(sectors=df_labels$seq, x=(df_labels$start + df_labels$end) / 2, labels=df_labels$gene,
                  cex=label_size / 12, font=ifelse(df_labels$strand == "-", 3, 1), side="outside",
                  line_col="#555555", line_lwd=0.4, connection_height=0.03, labels_height=0.07)
  }

  # --- outer scale -------------------------------------------
  circos.track(ylim=c(0, 1), track.height=0.02, bg.border=NA, panel.fun=function(x, y) {
    name <- CELL_META$sector.index
    L <- ls_length[name]
    circular <- ls_contigs[[name]]$circular

    major <- f_nice_interval(L, target_ticks=if (single) 20 else max(4, as.integer(20 * L / total)))
    ls_at <- seq(0, L, by=major)
    ls_labels <- ifelse(ls_at > 0, paste0(ls_at / 1000, " kb"), "0")

    # hide the end label of a single circular contig
    if (circular && single) {
      ls_labels[length(ls_labels)] <- ""
    }

    circos.axis(h="top", major.at=ls_at, labels=FALSE, minor.ticks=4, major.tick.length=1, lwd=0.8)

    # horizontal tick labels, anchored away from the circle
    for (i in which(ls_labels != "")) {
      theta <- circlize(ls_at[i], 2.5, sector.index=name)[1, "theta"] * pi / 180
      adj_x <- if (cos(theta) > 0.3) 0 else if (cos(theta) < -0.3) 1 else 0.5
      adj_y <- if (sin(theta) > 0.3) 0 else if (sin(theta) < -0.3) 1 else 0.5

      circos.text(ls_at[i], 2.5, ls_labels[i], facing="downward", adj=c(adj_x, adj_y), cex=6 / 12)
    }
  })

  # --- gene tracks: + strand outer, - strand inner -----------
  # function: draw the genes of one strand as a track
  f_gene_track <- function(which_strand) {
    circos.track(ylim=c(0, 1), track.height=0.09, bg.col="#f7f7f7", bg.border="#cccccc", bg.lwd=0.5,
                 panel.fun=function(x, y) {
      name <- CELL_META$sector.index
      df_sector <- df_genes[seq == name & strand == which_strand]
      deg_per_bp <- (CELL_META$cell.start.degree - CELL_META$cell.end.degree) / ls_length[name]

      for (i in seq_len(nrow(df_sector))) {
        start <- df_sector$start[i]
        end <- df_sector$end[i]
        colour <- df_sector$colour[i]

        if (df_sector$deg[i] < 1) {
          # too short for an arrow head
          circos.rect(start, 0, end, 1, col=colour, border="black", lwd=0.4)
        } else {
          circos.arrow(start, end, y=0.5, width=0.6, arrow.head.width=1,
                       arrow.head.length=min(1.5, df_sector$deg[i] / 2) / deg_per_bp,
                       arrow.position=if (which_strand == "+") "end" else "start",
                       col=colour, border="black", lwd=0.4)
        }
      }
    })
  }

  f_gene_track("+")
  f_gene_track("-")

  # --- GC content --------------------------------------------
  if (length(ls_gc) > 0) {
    # y-limits per sector
    ls_ylim <- t(sapply(ls_names, function(name) {
      if (is.null(ls_gc[[name]])) {
        return(c(0, 1))
      }

      y <- ls_gc[[name]]$y
      pad <- max(0.02, (max(y) - min(y)) * 0.1)
      return(c(max(0, min(y) - pad), min(1, max(y) + pad)))
    }))

    circos.track(ylim=ls_ylim, track.height=0.13, bg.col=NA, bg.lwd=0.5,
                 bg.border=ifelse(ls_names %in% names(ls_gc), "#999999", NA), panel.fun=function(x, y) {
      name <- CELL_META$sector.index
      if (is.null(ls_gc[[name]])) {
        return(invisible(NULL))
      }

      x <- ls_gc[[name]]$x
      y <- ls_gc[[name]]$y
      mean_gc <- mean(y)

      circos.lines(x, pmax(y, mean_gc), area=TRUE, baseline=mean_gc, col="#2a9d8f", border=NA)
      circos.lines(x, pmin(y, mean_gc), area=TRUE, baseline=mean_gc, col="#9b2226", border=NA)
      circos.lines(x, rep(mean_gc, length(x)), col="#444444", lwd=0.5)
    })
  }

  # --- plastid quadripartite structure -----------------------
  if (ir_found) {
    circos.track(ylim=c(0, 1), track.height=0.05, bg.col=NA, bg.lwd=0.5,
                 bg.border=ifelse(ls_names %in% names(ls_ir), "#999999", NA), panel.fun=function(x, y) {
      name <- CELL_META$sector.index
      df_regions <- ls_ir[[name]]
      if (is.null(df_regions)) {
        return(invisible(NULL))
      }

      L <- ls_length[name]
      for (i in seq_len(nrow(df_regions))) {
        s0 <- df_regions$start[i]
        s1 <- df_regions$end[i]
        if (f_arc_len(s0, s1, L) == 0) {
          next
        }

        # region ends exactly at the origin
        if (s1 == 0) {
          s1 <- L
        }

        if (s0 < s1) {
          circos.rect(s0, 0, s1, 1, col=df_regions$colour[i], border="black", lwd=0.3)
        } else {
          # wraps around the origin
          circos.rect(s0, 0, L, 1, col=df_regions$colour[i], border="black", lwd=0.3)
          circos.rect(0, 0, s1, 1, col=df_regions$colour[i], border="black", lwd=0.3)
        }
      }
    })
  }

  # --- region and sector labels ------------------------------
  if (ir_found || !single) {
    circos.track(ylim=c(0, 1), track.height=0.12, bg.border=NA, panel.fun=function(x, y) {
      name <- CELL_META$sector.index
      L <- ls_length[name]

      # LSC / IRB / SSC / IRA labels
      df_regions <- ls_ir[[name]]
      if (!is.null(df_regions)) {
        for (i in seq_len(nrow(df_regions))) {
          ln <- f_arc_len(df_regions$start[i], df_regions$end[i], L)
          if (ln == 0) {
            next
          }

          mid <- (df_regions$start[i] + ln / 2) %% L
          circos.text(mid, 0.75, sprintf("%s\n%.1f kb", df_regions$name[i], ln / 1000), facing="downward", cex=6 / 12)
        }
      }

      # sector label for multi-contig figures
      if (!single) {
        topo <- if (ls_contigs[[name]]$circular) "circular" else "linear"
        circos.text(CELL_META$xcenter, 0.2, paste0(name, "\n", f_fmt_bp(L), ", ", topo), facing="downward", cex=7 / 12)
      }
    })
  }

  # --- centre text -------------------------------------------
  org_name <- if (organelle == "pltd") "plastid" else "mitochondrion"
  if (single) {
    contig <- ls_contigs[[1]]
    centre <- paste0(org_name, "\n", ls_names[1], "\n", f_comma(contig$length), " bp\n",
                     if (contig$circular) "circular" else "linear")
  } else {
    centre <- paste0(org_name, "\n", length(ls_names), " contigs\n", f_comma(total), " bp total")
  }
  if (!is.null(title)) {
    centre <- paste0(title, "\n\n", centre)
  }
  text(0, 0, centre, cex=10 / 12)

  circos.clear()

  # --- legend ------------------------------------------------
  ls_legend <- ls_used$label
  ls_fill <- ls_used$colour
  ls_border <- rep("black", length(ls_legend))

  ls_legend <- c(ls_legend, paste0("GC above mean (", gc_window, " bp windows)"), "GC below mean")
  ls_fill <- c(ls_fill, "#2a9d8f", "#9b2226")
  ls_border <- c(ls_border, NA, NA)

  if (ir_found) {
    ls_legend <- c(ls_legend, "Inverted repeats (IRA/IRB)", "Single-copy regions (LSC/SSC)")
    ls_fill <- c(ls_fill, "#9e9e9e", "#e0e0e0")
    ls_border <- c(ls_border, "black", "black")
  }

  par(mar=c(0, 0, 0, 0))
  plot.new()
  legend("left", legend=ls_legend, fill=ls_fill, border=ls_border, bty="n", cex=7 / 12, y.intersp=1.2,
         title="outer ring: + strand\ninner ring: - strand (italic labels)\nscale: kb", title.adj=0)

  invisible(dev.off())
  message("wrote ", fn_out)
}

# --- Argument parsing ----------------------------------------
option_list <- list(
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="prefix.{mito,pltd}.ctg.fasta from oatk/pathfinder [required]", metavar="FILE"),
  make_option(c("-b", "--bed"), type="character", default=NULL,
              help="prefix.{mito,pltd}.ctg.bed from oatk/pathfinder [required]", metavar="FILE"),
  make_option(c("-o", "--out"), type="character", default="organelle_map.png",
              help="output image (.png, .svg, .pdf) [default %default]", metavar="FILE"),
  make_option(c("--organelle"), type="character", default="auto",
              help="gene colour scheme: auto, mito or pltd; auto = guess from the file name, then from gene names [default %default]"),
  make_option(c("--contig"), type="character", default=NULL,
              help="comma-separated contig name(s) to draw [default all]", metavar="NAMES"),
  make_option(c("--min-len"), dest="min_len", type="integer", default=0,
              help="skip contigs shorter than this [default %default]"),
  make_option(c("--min-score"), dest="min_score", type="double", default=0,
              help="skip BED features below this HMM score [default %default]"),
  make_option(c("--gc-window"), dest="gc_window", type="integer", default=1000,
              help="GC-content window in bp; 0 disables the track [default %default]"),
  make_option(c("--no-ir"), dest="no_ir", action="store_true", default=FALSE,
              help="plastid: do not detect/draw the LSC-IRB-SSC-IRA ring"),
  make_option(c("--no-trn"), dest="no_trn", action="store_true", default=FALSE,
              help="do not label tRNA genes"),
  make_option(c("--no-labels"), dest="no_labels", action="store_true", default=FALSE,
              help="draw genes without text labels"),
  make_option(c("--label-size"), dest="label_size", type="double", default=7,
              help="gene label font size [default %default]"),
  make_option(c("--title"), type="character", default=NULL,
              help="title printed in the centre of the map"),
  make_option(c("--figsize"), type="double", default=9,
              help="figure width/height in inches [default %default]"),
  make_option(c("--dpi"), type="integer", default=300,
              help="resolution of PNG output [default %default]"),
  make_option(c("--genbank"), type="character", default=NULL,
              help="also write the selected contigs + annotation as GenBank (for OGDRAW etc.)", metavar="FILE")
)

# parse the arguments
opt <- parse_args(OptionParser(option_list=option_list))

# stop the code if any required input is missing
if (is.null(opt$fasta) || is.null(opt$bed)) {
  stop("-f/--fasta and -b/--bed are required.\n
        Usage: Rscript oatk_circular_map.R -f x.mito.ctg.fasta -b x.mito.ctg.bed -o mito_map.png")
}

for (fn in c(opt$fasta, opt$bed)) {
  if (!file.exists(fn)) {
    stop(paste("File not found:", fn))
  }
}

if (!opt$organelle %in% c("auto", "mito", "pltd")) {
  stop("--organelle must be auto, mito or pltd")
}

# --- Read inputs ---------------------------------------------
ls_contigs <- f_read_fasta(opt$fasta)
df_bed <- f_read_bed(opt$bed, opt$min_score)
organelle <- if (opt$organelle != "auto") opt$organelle else f_detect_organelle(opt$fasta, df_bed)

# select contigs
if (!is.null(opt$contig)) {
  ls_select <- unlist(strsplit(opt$contig, split=","))

  ls_missing <- ls_select[!ls_select %in% names(ls_contigs)]
  if (length(ls_missing) > 0) {
    stop(paste0("contig(s) not found in FASTA: ", paste(ls_missing, collapse=", "), "\n",
                "available: ", paste(names(ls_contigs), collapse=", ")))
  }

  ls_contigs <- ls_contigs[ls_select]
}

if (opt$min_len > 0) {
  ls_contigs <- ls_contigs[sapply(ls_contigs, `[[`, "length") >= opt$min_len]
}

if (length(ls_contigs) == 0) {
  stop("no contigs left to draw")
}

message("organelle: ", organelle)
for (name in names(ls_contigs)) {
  contig <- ls_contigs[[name]]
  message(name, "\t", f_comma(contig$length), " bp\t", if (contig$circular) "circular" else "linear",
          "\t", sum(df_bed$seq == name), " features")
}

ls_unknown <- setdiff(unique(df_bed$seq), names(ls_contigs))
if (length(ls_unknown) > 0 && is.null(opt$contig)) {
  message("note: BED features on contigs absent from the FASTA were ignored: ", paste(sort(ls_unknown), collapse=", "))
}

# --- GenBank export ------------------------------------------
if (!is.null(opt$genbank)) {
  f_write_genbank(ls_contigs, df_bed, opt$genbank)
}

# --- Draw the map --------------------------------------------
f_draw(ls_contigs, df_bed, opt$out, organelle=organelle, title=opt$title, label_size=opt$label_size,
       no_trn=opt$no_trn, no_labels=opt$no_labels, gc_window=opt$gc_window,
       dpi=opt$dpi, figsize=opt$figsize, show_ir=!opt$no_ir)

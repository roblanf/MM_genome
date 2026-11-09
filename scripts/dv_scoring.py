#!/usr/bin/env python3
"""Parental assignment table from alignments of hap1/hap2 to two candidate parent references.

Rewrite of dv_scoring.sh. The difference: every contig is scored on ALL its primary alignments
to each reference, not on the single alignment with the most matches.

Per contig and per reference:
  aligned %  : fraction of the contig covered by primary alignments (merged intervals)
  identity % : length-weighted identity over those alignments (1 - de tag if present, else matches/block)
  chrom      : the reference chromosome receiving the most aligned bases
Derived: delta identity (decipiens - virginea), aligned ratio (decipiens / virginea), a call, the
partner contig in the other haplotype (largest summed alignment length in the hap1-vs-hap2 PAF),
the partner's call, and whether the two calls are opposite (reciprocal).

Usage:
  python3 dv_scoring.py --h1h2 hap1_vs_hap2.paf \
      --h1-dec hap1_vs_decipiens.paf --h1-vir hap1_vs_virginea.paf \
      --h2-dec hap2_vs_decipiens.paf --h2-vir hap2_vs_virginea.paf \
      [--blob1 hap1_stats.tsv --blob2 hap2_stats.tsv] \
      [--min-aln 1000] [--min-id 0.80] [--margin 0.5] [--min-aligned 0.2] [--all-alignments] \
      --md combined_parental_comparison.md --tsv combined_parental_comparison.tsv
"""
import argparse, csv, sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--h1h2", required=True, help="PAF of hap1 (query) vs hap2 (target)")
ap.add_argument("--h1-dec", required=True); ap.add_argument("--h1-vir", required=True)
ap.add_argument("--h2-dec", required=True); ap.add_argument("--h2-vir", required=True)
ap.add_argument("--blob1"); ap.add_argument("--blob2")
ap.add_argument("--min-aln", type=int, default=1000, help="min alignment block length (bp)")
ap.add_argument("--min-id", type=float, default=0.80, help="min identity of an alignment to be counted")
ap.add_argument("--margin", type=float, default=0.5, help="|delta identity| (percentage points) needed for a call")
ap.add_argument("--min-aligned", type=float, default=0.2, help="min aligned fraction to either reference for a call")
ap.add_argument("--all-alignments", action="store_true", help="count secondary alignments too (default: primary only)")
ap.add_argument("--md", default="combined_parental_comparison.md")
ap.add_argument("--tsv", default="combined_parental_comparison.tsv")
a = ap.parse_args()


def merged_len(ivs):
    ivs.sort()
    tot, cs, ce = 0, ivs[0][0], ivs[0][1]
    for s, e in ivs[1:]:
        if s > ce:
            tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs)


def score_paf(path):
    """contig -> dict(len, aligned_frac, identity, chrom, n)"""
    ivs = defaultdict(list); alen = defaultdict(float); match = defaultdict(float)
    per_t = defaultdict(lambda: defaultdict(int)); qlen = {}; n = defaultdict(int)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            tags = f[12:]
            if not a.all_alignments and any(t.startswith("tp:A:") for t in tags) and "tp:A:P" not in tags:
                continue
            q, ql, qs, qe, t, nm, bl = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[9]), int(f[10])
            if bl < a.min_aln:
                continue
            de = next((float(x[5:]) for x in tags if x.startswith("de:f:")), None)
            ident = 1.0 - de if de is not None else nm / bl
            if ident < a.min_id:
                continue
            qlen[q] = ql
            ivs[q].append((qs, qe)); alen[q] += bl; match[q] += ident * bl
            per_t[q][t] += qe - qs; n[q] += 1
    out = {}
    for q in qlen:
        out[q] = dict(len=qlen[q], aligned=merged_len(ivs[q]) / qlen[q], identity=100 * match[q] / alen[q],
                      chrom=max(per_t[q], key=per_t[q].get), n=n[q])
    return out


def partners(path):
    """best partner of each hap1 contig (and of each hap2 contig) by summed alignment block length"""
    s = defaultdict(int); l1, l2 = {}, {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            s[(f[0], f[5])] += int(f[10]); l1[f[0]] = int(f[1]); l2[f[5]] = int(f[6])
    b1, b2 = {}, {}
    for (q, t), v in s.items():
        if q not in b1 or v > s[(q, b1[q])]: b1[q] = t
        if t not in b2 or v > s[(b2[t], t)]: b2[t] = q
    return b1, b2, l1, l2


def blob(path):
    """contig -> (gc%, cov, phylum) from a blobtools stats table; tolerant to column names/positions"""
    if not path:
        return {}
    out = {}
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    hdr = [h.strip().lower() for h in rows[0]]
    def col(*names, default):
        for i, h in enumerate(hdr):
            if any(h == nm or h.startswith(nm) for nm in names):
                return i
        return default
    ci, cg, cc, cp = col("name", "id", "contig", default=1), col("gc", default=2), col("cov", default=4), col("phylum", "bestsum_phylum", default=5)
    for r in rows[1:]:
        if len(r) <= max(ci, cg, cc, cp) or r[0].strip().lower() == "index":
            continue
        try:
            gc = float(r[cg]); gc = gc * 100 if gc <= 1 else gc
            out[r[ci].strip()] = (gc, float(r[cc]), r[cp].strip())
        except ValueError:
            continue
    return out


def call(dec, vir):
    if dec is None and vir is None:
        return "no alignment"
    fa_d = dec["aligned"] if dec else 0.0; fa_v = vir["aligned"] if vir else 0.0
    if max(fa_d, fa_v) < a.min_aligned:
        return "unaligned"
    if dec is None: return "virginea"
    if vir is None: return "decipiens"
    d = dec["identity"] - vir["identity"]
    if d >= a.margin: return "decipiens"
    if d <= -a.margin: return "virginea"
    return "ambiguous"


h1d, h1v, h2d, h2v = score_paf(a.h1_dec), score_paf(a.h1_vir), score_paf(a.h2_dec), score_paf(a.h2_vir)
p1, p2, len1, len2 = partners(a.h1h2)
b1, b2 = blob(a.blob1), blob(a.blob2)

rows = []
for hap, dec_s, vir_s, part, bl, lens in (("hap1", h1d, h1v, p1, b1, len1), ("hap2", h2d, h2v, p2, b2, len2)):
    contigs = sorted(set(dec_s) | set(vir_s) | set(lens))
    for c in contigs:
        d, v = dec_s.get(c), vir_s.get(c)
        L = (d or v or {}).get("len") or lens.get(c, 0)
        gc, cov, phyl = bl.get(c, (None, None, "NA"))
        rows.append(dict(contig=c, hap=hap, length=L, gc=gc, cov=cov, phylum=phyl,
                         dec_chrom=d["chrom"] if d else "-", vir_chrom=v["chrom"] if v else "-",
                         dec_aligned=100 * d["aligned"] if d else 0.0, vir_aligned=100 * v["aligned"] if v else 0.0,
                         dec_id=d["identity"] if d else None, vir_id=v["identity"] if v else None,
                         delta=(d["identity"] - v["identity"]) if (d and v) else None,
                         ratio=(d["aligned"] / v["aligned"]) if (d and v and v["aligned"] > 0) else None,
                         call=call(d, v), partner=part.get(c, "-")))

call_of = {r["contig"]: r["call"] for r in rows}
delta_of = {r["contig"]: r["delta"] for r in rows}
for r in rows:
    p = r["partner"]
    r["partner_call"] = call_of.get(p, "-")
    r["partner_delta"] = delta_of.get(p)
    r["reciprocal"] = "yes" if {r["call"], r["partner_call"]} == {"decipiens", "virginea"} else "no"

fmt = lambda x, p=2: "-" if x is None else f"{x:.{p}f}"
hdr = ["Contig", "Hap", "Length", "GC%", "Cov", "Phylum", "decipiens chrom", "virginea chrom", "decipiens aligned %",
       "virginea aligned %", "decipiens ID %", "virginea ID %", "Δ ID", "Aligned ratio", "Call", "Partner",
       "Partner Δ ID", "Partner call", "Reciprocal"]
with open(a.tsv, "w") as t, open(a.md, "w") as m:
    t.write("\t".join(hdr) + "\n")
    cells = []
    for r in sorted(rows, key=lambda r: (r["hap"], -r["length"])):
        cells.append([r["contig"], r["hap"], f"{r['length']:,}", fmt(r["gc"], 1), fmt(r["cov"], 1), r["phylum"],
                      r["dec_chrom"], r["vir_chrom"], fmt(r["dec_aligned"], 1), fmt(r["vir_aligned"], 1),
                      fmt(r["dec_id"]), fmt(r["vir_id"]), ("-" if r["delta"] is None else f"{r['delta']:+.2f}"),
                      fmt(r["ratio"]), r["call"], r["partner"],
                      ("-" if r["partner_delta"] is None else f"{r['partner_delta']:+.2f}"), r["partner_call"], r["reciprocal"]])
    for c in cells:
        t.write("\t".join(c) + "\n")
    w = [max(len(h), *(len(c[i]) for c in cells)) for i, h in enumerate(hdr)] if cells else [len(h) for h in hdr]
    right = {2, 3, 4, 8, 9, 10, 11, 12, 13, 16}
    line = lambda c: "| " + " | ".join(x.rjust(w[i]) if i in right else x.ljust(w[i]) for i, x in enumerate(c)) + " |"
    m.write(line(hdr) + "\n")
    m.write("| " + " | ".join(("-" * (w[i] - 1) + ":") if i in right else (":" + "-" * (w[i] - 1)) for i in range(len(hdr))) + " |\n")
    for c in cells:
        m.write(line(c) + "\n")

n = len(rows)
from collections import Counter
print(f"{n} contigs scored; calls: {dict(Counter(r['call'] for r in rows))}; "
      f"reciprocal partner pairs: {sum(1 for r in rows if r['reciprocal']=='yes')//2}", file=sys.stderr)

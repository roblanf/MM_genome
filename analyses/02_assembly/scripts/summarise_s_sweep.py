#!/usr/bin/env python3
"""
Summarise BUSCO (compleasm) and k-mer (Merqury) results across the hifiasm -s sweep.

Expected layout (run from 02_assembly/homcov55/):
    compleasm/<label>/summary.txt
    merqury/<label>/<label>_merq.qv
    merqury/<label>/<label>_merq.completeness.stats
    merqury/<label>/<label>_merq.spectra-cn.hist

where <label> looks like s035_hap1, s035_hap2, s035_union, ...

Usage:
    python scripts/summarise_s_sweep.py [--compleasm compleasm] [--merqury merqury] [--het 28] [--hom 55]

Prints three Markdown tables: BUSCO completeness, Merqury QV / k-mer completeness,
and a copy-number breakdown of the het-peak and hom-peak k-mers from spectra-cn.hist.
"""

import argparse
import glob
import os
import re
from collections import defaultdict

ASM_ORDER = {"hap1": 0, "hap2": 1, "p_ctg": 2, "union": 3}


def sort_key(label):
    """Sort s030_hap1, s030_hap2, s030_union, s035_hap1, ...; unknown labels go last."""
    m = re.match(r"s(\d+)_(\w+)$", label)
    if not m:
        return (9.0, 99, label)
    s_val = int(m.group(1)) / 100.0
    return (s_val, ASM_ORDER.get(m.group(2), 50), label)


def print_table(header, rows, align):
    """header: list of column names; rows: list of lists; align: list of '<' or '>'."""
    widths = [max(len(str(h)), *(len(str(r[i])) for r in rows)) for i, h in enumerate(header)]
    fmt = " | ".join("{:" + a + str(w) + "}" for a, w in zip(align, widths))
    sep = " | ".join(("-" * w) if a == "<" else ("-" * (w - 1) + ":") for a, w in zip(align, widths))
    print("| " + fmt.format(*header) + " |")
    print("| " + sep + " |")
    for r in rows:
        print("| " + fmt.format(*[str(x) for x in r]) + " |")
    print()


# ---------------------------------------------------------------- compleasm --
def parse_compleasm(summary_file):
    """Return dict with S/D/F/I/M as 'xx.xx% (N)' and the total N."""
    stats, total = {}, "0"
    with open(summary_file) as f:
        for line in f:
            m = re.search(r"([SDFIM]):([\d.]+)%,\s+(\d+)", line)
            if m:
                key, perc, count = m.groups()
                stats[key] = f"{perc}% ({count})"
            n = re.search(r"^[Nn]:(\d+)", line.strip())
            if n:
                total = n.group(1)
    return stats, total


def busco_table(compleasm_dir):
    files = sorted(glob.glob(os.path.join(compleasm_dir, "*", "summary.txt")),
                   key=lambda p: sort_key(os.path.basename(os.path.dirname(p))))
    if not files:
        print(f"No compleasm summaries found under {compleasm_dir}/*/summary.txt\n")
        return
    rows = []
    for fn in files:
        label = os.path.basename(os.path.dirname(fn))
        stats, total = parse_compleasm(fn)
        rows.append([label] + [stats.get(k, "0.00% (0)") for k in "SDFIM"] + [total])
    print("### BUSCO completeness (compleasm, eudicotyledons_odb12)\n")
    print_table(["Dataset", "S (%, (N))", "D (%, (N))", "F (%, (N))", "I (%, (N))", "M (%, (N))", "Total"],
                rows, ["<", ">", ">", ">", ">", ">", ">"])


# ------------------------------------------------------------------ merqury --
def merqury_table(merqury_dir):
    """QV and k-mer completeness. Handles multi-line files (two-assembly mode)."""
    qv_files = sorted(glob.glob(os.path.join(merqury_dir, "*", "*_merq.qv")),
                      key=lambda p: sort_key(os.path.basename(os.path.dirname(p))))
    if not qv_files:
        print(f"No Merqury .qv files found under {merqury_dir}/*/*_merq.qv\n")
        return
    rows = []
    for qv_fn in qv_files:
        label = os.path.basename(os.path.dirname(qv_fn))
        comp_fn = qv_fn.replace(".qv", ".completeness.stats")

        # <asm>  asm_only_kmers  total_kmers  QV  error_rate
        qv = {}
        with open(qv_fn) as f:
            for line in f:
                p = line.split()
                if len(p) >= 5:
                    qv[p[0]] = p
        # <asm>  all  found  total  completeness%
        comp = {}
        if os.path.exists(comp_fn):
            with open(comp_fn) as f:
                for line in f:
                    p = line.split()
                    if len(p) >= 5:
                        comp[p[0]] = p

        for asm, p in qv.items():
            c = comp.get(asm, comp.get(asm.lower(), None))
            rows.append([
                label if len(qv) == 1 else f"{label} ({asm})",
                p[3],                                   # QV
                f"{float(p[4]):.2e}",                   # error rate
                c[4] if c else "NA",                    # completeness %
                f"{int(c[2]):,} / {int(c[3]):,}" if c else "NA",
                f"{int(p[1]):,}",                       # assembly-only (erroneous) k-mers
            ])
    print("### k-mer completeness and QV (Merqury, k=21)\n")
    print_table(["Dataset", "QV", "Error rate", "Completeness (%)", "Solid k-mers found / total", "Asm-only k-mers"],
                rows, ["<", ">", ">", ">", ">", ">"])


def copy_number_table(merqury_dir, het_peak, hom_peak):
    """
    Where do the het-peak and hom-peak read k-mers land in the assembly?
    From <label>_merq.spectra-cn.hist (columns: Copies, kmer_multiplicity, Count).
    Windows: het = [0.5*het, 1.5*het), hom = [0.8*hom, 1.25*hom).
    """
    # Two-assembly runs write <label>_merq.spectra-cn.hist (the union); single-assembly
    # runs only write <label>_merq.<asm>.spectra-cn.hist, so fall back to that.
    hist_files = []
    for d in sorted(glob.glob(os.path.join(merqury_dir, "*")), key=lambda p: sort_key(os.path.basename(p))):
        label = os.path.basename(d)
        union_hist = os.path.join(d, f"{label}_merq.spectra-cn.hist")
        if os.path.exists(union_hist):
            hist_files.append(union_hist)
        else:
            hist_files.extend(sorted(glob.glob(os.path.join(d, f"{label}_merq.*.spectra-cn.hist"))))
    if not hist_files:
        print(f"No Merqury spectra-cn.hist files found under {merqury_dir}\n")
        return
    het_lo, het_hi = 0.5 * het_peak, 1.5 * het_peak
    hom_lo, hom_hi = 0.8 * hom_peak, 1.25 * hom_peak

    def bucket(copies):
        if copies == "read-only":
            return "0"
        if copies in ("1", "2"):
            return copies
        return "3+"   # 3, 4, >4

    rows = []
    for fn in hist_files:
        label = os.path.basename(os.path.dirname(fn))
        het, hom = defaultdict(int), defaultdict(int)
        with open(fn) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 3 or p[0] == "Copies":
                    continue
                mult, count = float(p[1]), int(float(p[2]))
                if het_lo <= mult < het_hi:
                    het[bucket(p[0])] += count
                elif hom_lo <= mult < hom_hi:
                    hom[bucket(p[0])] += count
        het_tot, hom_tot = sum(het.values()) or 1, sum(hom.values()) or 1
        pct = lambda d, k, t: f"{100.0 * d[k] / t:5.1f}"
        rows.append([
            label,
            pct(het, "0", het_tot), pct(het, "1", het_tot), f"{100.0 * (het['2'] + het['3+']) / het_tot:5.1f}",
            pct(hom, "0", hom_tot), pct(hom, "1", hom_tot), pct(hom, "2", hom_tot), pct(hom, "3+", hom_tot),
        ])
    print(f"### Copy number of het-peak ({het_lo:.0f}-{het_hi:.0f}x) and hom-peak ({hom_lo:.0f}-{hom_hi:.0f}x) k-mers in the assembly (%)\n")
    print_table(["Dataset", "het: 0x", "het: 1x", "het: >=2x", "hom: 0x", "hom: 1x", "hom: 2x", "hom: >=3x"],
                rows, ["<", ">", ">", ">", ">", ">", ">", ">"])
    print("Expected for one haplotype: het ~50% 0x / ~50% 1x, hom ~100% 1x. "
          "Expected for hap1+hap2 together: het ~100% 1x, hom ~100% 2x; "
          "het k-mers at >=2x or hom k-mers at >=3x are duplicated regions.\n")


# --------------------------------------------------------------------- main --
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compleasm", default="compleasm", help="directory holding <label>/summary.txt")
    ap.add_argument("--merqury", default="merqury", help="directory holding <label>/<label>_merq.*")
    ap.add_argument("--het", type=float, default=28, help="heterozygous k-mer peak (GenomeScope kcov)")
    ap.add_argument("--hom", type=float, default=55, help="homozygous k-mer peak (2 x kcov)")
    args = ap.parse_args()

    busco_table(args.compleasm)
    merqury_table(args.merqury)
    copy_number_table(args.merqury, args.het, args.hom)

#!/usr/bin/env python3
"""
oatk_circular_map.py - draw a circular gene map of an oatk organelle assembly
(mitochondrion or plastid).

Inputs are the two files written by oatk / pathfinder:
    prefix.mito.ctg.fasta + prefix.mito.ctg.bed     (mitochondrion)
    prefix.pltd.ctg.fasta + prefix.pltd.ctg.bed     (plastid)
FASTA headers carry length=, circular=, path=; the BED is BED6
(seq, start, end, gene, score, strand) with a '#' header line.

Examples:
    python oatk_circular_map.py x.mito.ctg.fasta x.mito.ctg.bed -o mito_map.png
    python oatk_circular_map.py x.pltd.ctg.fasta x.pltd.ctg.bed -o pltd_map.svg
    python oatk_circular_map.py x.mito.ctg.fasta x.mito.ctg.bed --contig ctg000001c --genbank mito.gb

Requires: pycirclize (pip install pycirclize), numpy, matplotlib.
"""

import argparse
import gzip
import math
import re
import sys
from collections import OrderedDict

import numpy as np
from matplotlib.patches import Patch
from pycirclize import Circos

# ---------------------------------------------------------------------------
# Gene categories (OatkDB embryophyta gene names) and colours
# ---------------------------------------------------------------------------
MITO_CATEGORIES = [
    # (legend label, colour, regex on gene name)
    ("NADH dehydrogenase (complex I)",            "#ffd166", r"^nad"),
    ("Succinate dehydrogenase (II)",              "#ef8354", r"^sdh"),
    ("Ubiquinol cytochrome c reductase (III)",    "#d62828", r"^cob"),
    ("Cytochrome c oxidase (IV)",                 "#c77dff", r"^cox"),
    ("ATP synthase (V)",                          "#06a77d", r"^atp"),
    ("Cytochrome c biogenesis",                   "#a0522d", r"^ccm"),
    ("Ribosomal proteins",                        "#118ab2", r"^rp[sl]"),
    ("Maturase / other proteins",                 "#6c757d", r"^(matR|mttB|rtl|orf)"),
    ("rRNA",                                      "#ff006e", r"^rrn"),
    ("tRNA",                                      "#2b2d42", r"^trn"),
]
PLTD_CATEGORIES = [
    ("Photosystem I",                             "#1b5e20", r"^psa"),
    ("Photosystem II",                            "#4caf50", r"^psb"),
    ("Cytochrome b/f complex",                    "#a5d6a7", r"^pet"),
    ("ATP synthase",                              "#00acc1", r"^atp"),
    ("NADH dehydrogenase",                        "#ffd166", r"^ndh"),
    ("RubisCO large subunit",                     "#5d4037", r"^rbcL"),
    ("RNA polymerase",                            "#d62828", r"^rpo"),
    ("Ribosomal proteins",                        "#118ab2", r"^rp[sl]"),
    ("Chlorophyll biosynthesis",                  "#8e24aa", r"^chl"),
    ("Other genes",                               "#6c757d", r"^(matK|clpP|accD|cemA|ccsA|infA|tufA|cys)"),
    ("Hypothetical ORFs (ycf)",                   "#adb5bd", r"^ycf"),
    ("rRNA",                                      "#ff006e", r"^rrn"),
    ("tRNA",                                      "#2b2d42", r"^trn"),
]
OTHER = ("Other", "#cccccc")
PLTD_MARKERS = re.compile(r"^(psa|psb|pet|ndh|rbcL|rpo|ycf|matK|clpP|accD|cemA|ccsA|rrn16|rrn23|rrn4\.5)", re.I)


def make_classifier(categories):
    def category(gene):
        for label, colour, pattern in categories:
            if re.match(pattern, gene, flags=re.IGNORECASE):
                return label, colour
        return OTHER
    return category


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------
def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    """Return OrderedDict name -> dict(seq, length, circular, path, header)."""
    contigs = OrderedDict()
    name, header, chunks = None, None, []

    def flush():
        if name is None:
            return
        seq = "".join(chunks).upper()
        tags = dict(re.findall(r"(\w+)=(\S+)", header))
        length = int(tags.get("length", len(seq)))
        if "circular" in tags:
            circular = tags["circular"].lower() == "true"
        else:
            circular = name.endswith("c")  # oatk naming: ctg000001c / ctg000001l
        contigs[name] = dict(seq=seq, length=length, circular=circular,
                             path=tags.get("path", ""), header=header)

    with open_any(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:]
                name = header.split()[0]
                chunks = []
            elif line:
                chunks.append(line.strip())
    flush()
    return contigs


def read_bed(path, min_score=0):
    """Return dict contig -> list of (start, end, gene, score, strand)."""
    feats = {}
    with open_any(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                f = line.split()
            if len(f) < 6:
                continue
            seq, start, end, gene, score, strand = f[:6]
            start, end = int(start), int(end)
            try:
                score = float(score)
            except ValueError:
                score = 0.0
            if score < min_score:
                continue
            if start > end:
                start, end = end, start
            feats.setdefault(seq, []).append((start, end, gene, score, strand))
    for v in feats.values():
        v.sort()
    return feats


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def nice_interval(length, target_ticks=20):
    """Pick a 1/2/5 x 10^k tick interval giving roughly target_ticks ticks."""
    raw = length / target_ticks
    exp = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 5, 10):
        if m * exp >= raw:
            return m * exp
    return 10 * exp


def gc_content(seq, window):
    """Sliding-window GC fraction; returns (window centres, gc values)."""
    n = len(seq)
    if n == 0:
        return np.array([]), np.array([])
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    is_gc = ((arr == ord("G")) | (arr == ord("C"))).astype(np.int32)
    is_base = np.isin(arr, [ord(b) for b in "ACGT"]).astype(np.int32)
    cs_gc = np.concatenate([[0], np.cumsum(is_gc)])
    cs_b = np.concatenate([[0], np.cumsum(is_base)])
    centres, values = [], []
    for s in range(0, n, window):
        e = min(s + window, n)
        bases = cs_b[e] - cs_b[s]
        if bases == 0:
            continue
        centres.append((s + e) / 2)
        values.append((cs_gc[e] - cs_gc[s]) / bases)
    return np.array(centres), np.array(values)


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def find_inverted_repeat(seq, circular=True, k=31, min_len=1000, max_gap=100):
    """Locate the largest pair of inverted repeats (plastid IRA/IRB).

    Returns ((a_start, a_end), (b_start, b_end)) as 0-based half-open intervals on
    the original sequence (an interval with start > end wraps around the origin),
    or None if no inverted repeat >= min_len is found.
    """
    L = len(seq)
    if L < 2 * min_len:
        return None
    s = seq + seq if circular else seq          # doubling lets a repeat cross the origin
    n = len(s)
    fwd = {}
    for p in range(n - k + 1):
        fwd.setdefault(s[p:p + k], []).append(p)
    rc = revcomp(s)
    diag = {}
    for i in range(n - k + 1):
        for p in fwd.get(rc[i:i + k], ()):
            q = n - i - k                       # forward position of the partner k-mer
            if abs(p - q) < k:                  # palindromic self-match
                continue
            diag.setdefault(p - i, []).append(p)

    best = None
    for d, ps in diag.items():
        ps.sort()
        run_start = prev = ps[0]
        for p in ps[1:] + [None]:
            if p is None or p - prev > max_gap:
                if prev + k - run_start >= min_len and (best is None or prev + k - run_start > best[0]):
                    best = (prev + k - run_start, run_start, prev + k, d)
                if p is not None:
                    run_start = p
            if p is not None:
                prev = p
    if best is None:
        return None
    length, a0, a1, d = best
    if length > L / 2 or (circular and length >= L - 1):
        return None
    # partner copy: i = p - d, forward position q = n - i - k = n - (p - d) - k
    b0, b1 = n - (a1 - k - d) - k, n - (a0 - d) - k + k
    if circular:
        a0, a1, b0, b1 = a0 % L, a1 % L, b0 % L, b1 % L
    return (a0, a1), (b0, b1)


def arc_len(start, end, L):
    return (end - start) % L if start != end else 0


def fmt_bp(x):
    if x >= 1_000_000:
        return f"{x/1e6:g} Mb"
    if x >= 1000:
        return f"{x/1000:g} kb"
    return f"{x:g} bp"


def detect_organelle(fasta_path, feats):
    name = fasta_path.lower()
    if "pltd" in name or "plastid" in name or "chloro" in name:
        return "pltd"
    if "mito" in name:
        return "mito"
    genes = [g for v in feats.values() for _, _, g, _, _ in v]
    return "pltd" if any(PLTD_MARKERS.match(g) for g in genes) else "mito"


# ---------------------------------------------------------------------------
# Optional GenBank export (e.g. for OGDRAW, Geneious, SnapGene)
# ---------------------------------------------------------------------------
def write_genbank(contigs, feats, path):
    import datetime
    date = datetime.date.today().strftime("%d-%b-%Y").upper()
    with open(path, "w") as out:
        for name, c in contigs.items():
            L, seq = c["length"], c["seq"]
            topo = "circular" if c["circular"] else "linear"
            out.write(f"LOCUS       {name:<23} {L:>7} bp    DNA     {topo:<8} PLN {date}\n")
            out.write(f"DEFINITION  {name} oatk organelle assembly, {topo}.\n")
            out.write(f"ACCESSION   {name}\nVERSION     {name}\n")
            out.write("KEYWORDS    .\nSOURCE      .\n  ORGANISM  .\n            .\n")
            out.write("FEATURES             Location/Qualifiers\n")
            out.write(f"     source          1..{L}\n")
            if c["path"]:
                out.write(f'                     /note="oatk path {c["path"]}"\n')
            for start, end, gene, score, strand in feats.get(name, []):
                # oatk writes nhmmscan's 1-based alignment coordinates into the BED
                # (not 0-based BED convention), so they are used as-is; ends can be +/-1 bp
                start, end = max(1, start), min(L, end)
                if end < start:
                    continue
                loc = f"{start}..{end}" if strand == "+" else f"complement({start}..{end})"
                g = gene.lower()
                key = "tRNA" if g.startswith("trn") else "rRNA" if g.startswith("rrn") else "CDS"
                out.write(f"     gene            {loc}\n                     /gene=\"{gene}\"\n")
                out.write(f"     {key:<16}{loc}\n                     /gene=\"{gene}\"\n")
                out.write(f'                     /note="oatk hmmannot score {score:g}"\n')
            out.write("ORIGIN\n")
            for i in range(0, len(seq), 60):
                chunk = seq[i:i + 60].lower()
                out.write(f"{i + 1:>9} " + " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10)) + "\n")
            out.write("//\n")
    print(f"wrote {path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def draw(contigs, feats, out, organelle="mito", title=None, label_size=7, no_trn=False,
         no_labels=False, gc_window=1000, dpi=300, figsize=9.0, show_ir=True):
    categories = PLTD_CATEGORIES if organelle == "pltd" else MITO_CATEGORIES
    category = make_classifier(categories)

    names = list(contigs)
    single = len(names) == 1
    all_circular = all(contigs[n]["circular"] for n in names)
    space = 0 if single else 4
    total = sum(contigs[n]["length"] for n in names)

    sectors = {n: contigs[n]["length"] for n in names}
    circos = Circos(sectors, space=space, endspace=not (single and all_circular))
    used_categories = OrderedDict()
    ir_found = False

    for sector in circos.sectors:
        name = sector.name
        c = contigs[name]
        L, seq = c["length"], c["seq"]
        sector_feats = feats.get(name, [])

        # --- outer scale ---------------------------------------------------
        scale = sector.add_track((88, 90))
        scale.axis(fc="none", ec="black", lw=0.8)
        major = nice_interval(L, target_ticks=20 if single else max(4, int(20 * L / total)))
        scale.xticks_by_interval(major, tick_length=2, label_size=6, label_margin=1.5,
                                 label_formatter=lambda v: f"{v/1000:g} kb" if v > 0 else "0",
                                 show_endlabel=not c["circular"] or not single)
        scale.xticks_by_interval(major / 5, tick_length=1, show_label=False)

        # --- gene tracks: + strand outer, - strand inner -------------------
        plus = sector.add_track((78, 87))
        minus = sector.add_track((68, 77))
        plus.axis(fc="#f7f7f7", ec="#cccccc", lw=0.5)
        minus.axis(fc="#f7f7f7", ec="#cccccc", lw=0.5)

        for start, end, gene, score, strand in sector_feats:
            start, end = max(0, min(start, L)), max(0, min(end, L))
            if end <= start:
                continue
            label, colour = category(gene)
            used_categories[label] = colour
            track = plus if strand == "+" else minus
            deg = 360 * (end - start) / total
            kw = dict(fc=colour, ec="black", lw=0.4)
            if deg < 1.0:                       # too short for an arrow head
                track.rect(start, end, **kw)
            elif strand == "+":
                track.arrow(start, end, head_length=min(1.5, deg / 2), shaft_ratio=0.6, **kw)
            else:
                track.arrow(end, start, head_length=min(1.5, deg / 2), shaft_ratio=0.6, **kw)

            if not no_labels and not (no_trn and gene.lower().startswith("trn")):
                plus.annotate((start + end) / 2, gene, min_r=90.5, max_r=97,
                              label_size=label_size,
                              line_kws=dict(color="#555555", lw=0.4),
                              text_kws=dict(fontstyle="italic" if strand == "-" else "normal"))

        # --- GC content ----------------------------------------------------
        if seq and gc_window > 0 and len(seq) > gc_window:
            gc = sector.add_track((52, 65))
            gc.axis(fc="none", ec="#999999", lw=0.5)
            x, y = gc_content(seq, gc_window)
            if len(y):
                mean = float(np.mean(y))
                pad = max(0.02, (y.max() - y.min()) * 0.1)
                vmin, vmax = max(0, y.min() - pad), min(1, y.max() + pad)
                gc.fill_between(x, np.maximum(y, mean), mean, vmin=vmin, vmax=vmax,
                                fc="#2a9d8f", ec="none", alpha=0.9)
                gc.fill_between(x, np.minimum(y, mean), mean, vmin=vmin, vmax=vmax,
                                fc="#9b2226", ec="none", alpha=0.9)
                gc.line(x, [mean] * len(x), vmin=vmin, vmax=vmax, color="#444444", lw=0.5)

        # --- plastid quadripartite structure (LSC / IRB / SSC / IRA) -------
        if organelle == "pltd" and show_ir and seq and L >= 20000:
            ir = find_inverted_repeat(seq, circular=c["circular"])
            if ir:
                ir_found = True
                (a0, a1), (b0, b1) = ir
                qt = sector.add_track((45, 50))
                qt.axis(fc="none", ec="#999999", lw=0.5)
                if c["circular"]:
                    # clockwise from the end of one copy you reach the start of the other:
                    # the two single-copy regions are (a1 -> b0) and (b1 -> a0)
                    sc1, sc2 = (a1, b0), (b1, a0)
                    if arc_len(*sc1, L) >= arc_len(*sc2, L):
                        lsc, ssc, irb, ira = sc1, sc2, (b0, b1), (a0, a1)
                    else:
                        lsc, ssc, irb, ira = sc2, sc1, (a0, a1), (b0, b1)
                    regions = [("LSC", lsc, "#e0e0e0"), ("IRB", irb, "#9e9e9e"),
                               ("SSC", ssc, "#e0e0e0"), ("IRA", ira, "#9e9e9e")]
                else:
                    if a0 > b0:
                        (a0, a1), (b0, b1) = (b0, b1), (a0, a1)
                    lsc, ssc, irb, ira = (0, a0), (a1, b0), (a0, a1), (b0, b1)
                    regions = [("IR", irb, "#9e9e9e"), ("IR", ira, "#9e9e9e")]
                for reg_name, (s0, s1), col in regions:
                    ln = arc_len(s0, s1, L)
                    if ln == 0:
                        continue
                    if s1 == 0:                             # region ends exactly at the origin
                        s1 = L
                    if s0 < s1:
                        qt.rect(s0, s1, fc=col, ec="black", lw=0.3)
                    else:                                   # wraps around the origin
                        qt.rect(s0, L, fc=col, ec="black", lw=0.3)
                        qt.rect(0, s1, fc=col, ec="black", lw=0.3)
                    mid = (s0 + ln / 2) % L
                    sector.text(f"{reg_name}\n{ln/1000:.1f} kb", x=mid, r=40, size=6)
                print(f"{name}\tinverted repeat: IRA {ira[0]+1}-{ira[1] or L} / IRB {irb[0]+1}-{irb[1] or L} "
                      f"({arc_len(*irb, L):,} bp), LSC {arc_len(*lsc, L):,} bp, SSC {arc_len(*ssc, L):,} bp",
                      file=sys.stderr)
            else:
                print(f"{name}\tno inverted repeat >= 1 kb found", file=sys.stderr)

        # --- sector label for multi-contig figures -------------------------
        if not single:
            topo = "circular" if c["circular"] else "linear"
            sector.text(f"{name}\n{fmt_bp(L)}, {topo}", r=33, size=7)

    # --- centre text and legend --------------------------------------------
    org_name = "plastid" if organelle == "pltd" else "mitochondrion"
    if single:
        c = contigs[names[0]]
        topo = "circular" if c["circular"] else "linear"
        centre = f"{org_name}\n{names[0]}\n{c['length']:,} bp\n{topo}"
    else:
        centre = f"{org_name}\n{len(names)} contigs\n{total:,} bp total"
    if title:
        centre = title + "\n\n" + centre
    circos.text(centre, r=0, size=10, ha="center", va="center")

    fig = circos.plotfig(figsize=(figsize, figsize), dpi=dpi)
    order = [lab for lab, _, _ in categories] + [OTHER[0]]
    handles = [Patch(facecolor=used_categories[lab], edgecolor="black", lw=0.4, label=lab)
               for lab in order if lab in used_categories]
    handles.append(Patch(facecolor="#2a9d8f", label=f"GC above mean ({gc_window} bp windows)"))
    handles.append(Patch(facecolor="#9b2226", label="GC below mean"))
    if ir_found:
        handles.append(Patch(facecolor="#9e9e9e", edgecolor="black", lw=0.3, label="Inverted repeats (IRA/IRB)"))
        handles.append(Patch(facecolor="#e0e0e0", edgecolor="black", lw=0.3, label="Single-copy regions (LSC/SSC)"))
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.97, 0.5),
               fontsize=7, frameon=False, title_fontsize=7, alignment="left",
               title="outer ring: + strand\ninner ring: - strand (italic labels)\nscale: kb")
    fig.savefig(out, dpi=dpi, bbox_inches="tight", pad_inches=0.3)
    print(f"wrote {out}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="prefix.{mito,pltd}.ctg.fasta from oatk/pathfinder")
    ap.add_argument("bed", help="prefix.{mito,pltd}.ctg.bed from oatk/pathfinder")
    ap.add_argument("-o", "--out", default="organelle_map.png",
                    help="output image (.png, .svg, .pdf) [organelle_map.png]")
    ap.add_argument("--organelle", choices=["auto", "mito", "pltd"], default="auto",
                    help="gene colour scheme; auto = guess from the file name, then from gene names [auto]")
    ap.add_argument("--contig", nargs="+", help="draw only these contig name(s) [all]")
    ap.add_argument("--min-len", type=int, default=0, help="skip contigs shorter than this [0]")
    ap.add_argument("--min-score", type=float, default=0, help="skip BED features below this HMM score [0]")
    ap.add_argument("--gc-window", type=int, default=1000, help="GC-content window in bp; 0 disables the track [1000]")
    ap.add_argument("--no-ir", action="store_true", help="plastid: do not detect/draw the LSC-IRB-SSC-IRA ring")
    ap.add_argument("--no-trn", action="store_true", help="do not label tRNA genes")
    ap.add_argument("--no-labels", action="store_true", help="draw genes without text labels")
    ap.add_argument("--label-size", type=float, default=7, help="gene label font size [7]")
    ap.add_argument("--title", help="title printed in the centre of the map")
    ap.add_argument("--figsize", type=float, default=9, help="figure width/height in inches [9]")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--genbank", metavar="FILE",
                    help="also write the selected contigs + annotation as GenBank (for OGDRAW etc.)")
    args = ap.parse_args()

    contigs = read_fasta(args.fasta)
    feats = read_bed(args.bed, args.min_score)
    organelle = args.organelle if args.organelle != "auto" else detect_organelle(args.fasta, feats)

    if args.contig:
        missing = [n for n in args.contig if n not in contigs]
        if missing:
            sys.exit(f"contig(s) not found in FASTA: {', '.join(missing)}\n"
                     f"available: {', '.join(contigs)}")
        contigs = OrderedDict((n, contigs[n]) for n in args.contig)
    if args.min_len:
        contigs = OrderedDict((n, c) for n, c in contigs.items() if c["length"] >= args.min_len)
    if not contigs:
        sys.exit("no contigs left to draw")

    print(f"organelle: {organelle}", file=sys.stderr)
    for n, c in contigs.items():
        ngenes = len(feats.get(n, []))
        print(f"{n}\t{c['length']:,} bp\t{'circular' if c['circular'] else 'linear'}\t{ngenes} features",
              file=sys.stderr)
    unknown = set(feats) - set(contigs)
    if unknown and not args.contig:
        print(f"note: BED features on contigs absent from the FASTA were ignored: {', '.join(sorted(unknown))}",
              file=sys.stderr)

    if args.genbank:
        write_genbank(contigs, feats, args.genbank)

    draw(contigs, feats, args.out, organelle=organelle, title=args.title, label_size=args.label_size,
         no_trn=args.no_trn, no_labels=args.no_labels, gc_window=args.gc_window,
         dpi=args.dpi, figsize=args.figsize, show_ir=not args.no_ir)


if __name__ == "__main__":
    main()

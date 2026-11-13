#!/usr/bin/env python3
"""Draw scaffolds as bars: contigs coloured by haplotype of origin, gaps and fills marked.
Before = RagTag AGP; after (optional) = TGS-GapCloser gap_fill_detail.

Usage:
  python3 plot_scaffolds.py --agp ragtag.scaffold.agp [--detail closed.gap_fill_detail --prefix vir_]
                            [--min-len 10000000] [--title "virginea-like"] --out scaffolds.png
"""
import argparse
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ap = argparse.ArgumentParser()
ap.add_argument("--agp", required=True)
ap.add_argument("--detail")
ap.add_argument("--fasta", help="the scaffold FASTA given to TGS-GapCloser; N runs inside contigs are then shown on the before track")
ap.add_argument("--prefix", default="", help="prefix added to scaffold names in the detail file (e.g. vir_)")
ap.add_argument("--min-len", type=int, default=10_000_000, help="draw scaffolds at least this long")
ap.add_argument("--title", default="")
ap.add_argument("--out", default="scaffolds.png")
a = ap.parse_args()

COL = {"h1": "#d95f02", "h2": "#1b9e77", "other": "#999999"}   # hap1 orange, hap2 green


def hap(contig):
    return "h1" if contig.startswith("h1") else "h2" if contig.startswith("h2") else "other"


# ---- before: AGP --------------------------------------------------------------------------
before = defaultdict(list)          # scaffold -> [(start, end, kind, contig)]
for line in open(a.agp):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    s, e = int(f[1]), int(f[2])
    if f[4] == "W":
        before[f[0]].append((s, e, "contig", f[5]))
    elif f[4] in ("U", "N"):
        before[f[0]].append((s, e, "gap", ""))
scaffolds = [s for s in before if max(e for _, e, _, _ in before[s]) >= a.min_len]
scaffolds.sort(key=lambda s: -max(e for _, e, _, _ in before[s]))


# ---- N runs inside the input scaffolds (optional): shown on the before track ---------------
nruns = defaultdict(list)
if a.fasta:
    name, pos, in_n, n_start = None, 0, False, 0
    def _flush():
        if in_n:
            nruns[name].append((n_start, pos))
    with open(a.fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                _flush(); name = line[1:].split()[0]
                if a.prefix and name.startswith(a.prefix):
                    name = name[len(a.prefix):]
                pos, in_n = 0, False
                continue
            for ch in line.strip():
                pos += 1
                if ch in "Nn":
                    if not in_n:
                        in_n, n_start = True, pos
                elif in_n:
                    nruns[name].append((n_start, pos - 1)); in_n = False
        _flush()
    for scaf in list(before):
        for s_, e_ in nruns.get(scaf, []):
            if e_ - s_ + 1 >= 10:
                before[scaf].append((s_, e_, "nrun", ""))


def contig_at(scaf, pos):
    """contig occupying position pos (1-based) of the original scaffold"""
    for s, e, kind, c in before[scaf]:
        if kind == "contig" and s <= pos <= e:
            return c
    return ""


# ---- after: gap_fill_detail -----------------------------------------------------------------
after = defaultdict(list)           # scaffold -> [(start, end, kind, contig)]
if a.detail:
    cur = None
    for line in open(a.detail):
        if line.startswith(">"):
            cur = line[1:].strip()
            if a.prefix and cur.startswith(a.prefix):
                cur = cur[len(a.prefix):]
            continue
        f = line.split()
        if cur is None or len(f) < 3:
            continue
        s, e, kind = int(f[0]), int(f[1]), f[2]
        if kind == "S":
            # split the kept segment at contig boundaries of the original scaffold
            os_, oe = int(f[3]), int(f[4])
            for cs, ce, k, c in before[cur]:
                if k != "contig" or ce < os_ or cs > oe:
                    continue
                lo, hi = max(cs, os_), min(ce, oe)
                after[cur].append((s + (lo - os_), s + (hi - os_), "contig", c))
        elif kind == "N":
            after[cur].append((s, e, "gap", ""))
        elif kind == "F":
            after[cur].append((s, e, "fill", ""))

# ---- draw -----------------------------------------------------------------------------------
rows = len(scaffolds)
fig_h = max(3, 0.9 * rows + 1.5)
fig, ax = plt.subplots(figsize=(14, fig_h))
L = max(e for s in scaffolds for _, e, _, _ in before[s])
bar_h, step = 0.32, 1.0
for i, scaf in enumerate(scaffolds):
    y0 = (rows - 1 - i) * step
    tracks = [("before", before[scaf], y0 + 0.36)] + ([("after", after.get(scaf, before[scaf]), y0)] if a.detail else [])
    for label, segs, y in tracks:
        for s, e, kind, c in segs:
            if kind == "contig":
                ax.add_patch(Rectangle((s, y), e - s + 1, bar_h, color=COL[hap(c)], lw=0))
            elif kind == "fill":   # drawn slightly taller than the bar and at least 0.5% of the axis wide, so kb-scale fills stay visible
                ax.add_patch(Rectangle((s, y - 0.05), max(e - s + 1, 0.005 * L), bar_h + 0.10, color="#7570b3", lw=0, zorder=3))
        for s, e, kind, c in segs:
            if kind == "nrun":       # N run inside a contig: blank it out so both tracks describe the sequence the same way
                ax.add_patch(Rectangle((s, y - 0.01), e - s + 1, bar_h + 0.02, color="white", lw=0, zorder=2))
        for s, e, kind, c in segs:
            if kind in ("gap", "nrun"):
                ax.plot([s, s], [y - 0.04, y + bar_h + 0.04], color="#d62728", lw=1.2, zorder=4)
        ax.text(-0.01 * L, y + bar_h / 2, label, ha="right", va="center", fontsize=8, color="#555555")
    ax.text(-0.09 * L, y0 + 0.36 + bar_h / 2, scaf.replace("_RagTag", ""), ha="right", va="center", fontsize=9)
ax.set_xlim(-0.14 * L, L * 1.01)
ax.set_ylim(-0.3, rows * step + 0.3)
ax.set_yticks([])
ax.set_xlabel("position (Mb)")
ax.set_xticks([x for x in range(0, int(L) + 1, 10_000_000)])
ax.set_xticklabels([f"{x // 1_000_000}" for x in ax.get_xticks()])
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
handles = [Rectangle((0, 0), 1, 1, color=COL["h1"]), Rectangle((0, 0), 1, 1, color=COL["h2"]),
           Rectangle((0, 0), 1, 1, color="#7570b3"), plt.Line2D([0], [0], color="#d62728", lw=1.2)]
ax.legend(handles, ["contig from hap1", "contig from hap2", "fill (TGS-GapCloser)", "gap"],
          loc="lower right", frameon=False, fontsize=8, ncol=4)
if a.title:
    ax.set_title(a.title, loc="left")
plt.tight_layout()
plt.savefig(a.out, dpi=200)
print(f"wrote {a.out}: {rows} scaffolds")

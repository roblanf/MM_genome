#!/usr/bin/env python3
"""Flag contigs that are (mostly) contained in a LONGER contig of the same haplotype.

Input: PAF from a self-alignment, e.g.
    minimap2 -cx asm20 -DP -m200 -t16 hap1.fa hap1.fa > h1_self.paf

Output (TSV, one row per flagged contig, most-covered first):
    contig  length  frac_covered  max_identity  contained_in

Usage:
    python3 flag_contained_contigs.py h1_self.paf > h1_dups.tsv
"""
import argparse
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("paf")
ap.add_argument("--min-aln", type=int, default=10000,
                help="min alignment block length in bp (default 10000)")
ap.add_argument("--min-id", type=float, default=0.90,
                help="min identity, 1 - de tag or matches/alen (default 0.90)")
ap.add_argument("--min-frac", type=float, default=0.80,
                help="flag contig if >= this fraction of its length is covered (default 0.80)")
a = ap.parse_args()

ivs = defaultdict(list)   # contig -> [(qstart, qend, identity, target)]
length = {}

with open(a.paf) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, ql, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
        t, tl, nmatch, alen = f[5], int(f[6]), int(f[9]), int(f[10])
        if q == t or alen < a.min_aln:
            continue
        # only count alignments to a longer contig, so the shorter copy of a
        # duplicated pair is flagged, not both (ties broken by name)
        if tl < ql or (tl == ql and t < q):
            continue
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        ident = 1.0 - de if de is not None else nmatch / alen
        if ident < a.min_id:
            continue
        length[q] = ql
        ivs[q].append((qs, qe, ident, t))

rows = []
for q, iv in ivs.items():
    iv.sort()
    covered, cs, ce = 0, iv[0][0], iv[0][1]
    for s, e, _, _ in iv[1:]:
        if s > ce:
            covered += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    covered += ce - cs
    frac = covered / length[q]
    if frac < a.min_frac:
        continue
    per_target = defaultdict(int)
    for s, e, _, t in iv:
        per_target[t] += e - s
    top = sorted(per_target, key=per_target.get, reverse=True)[:5]
    rows.append((q, length[q], frac, max(i for _, _, i, _ in iv), ",".join(top)))

print("contig\tlength\tfrac_covered\tmax_identity\tcontained_in")
for q, ln, frac, ident, tg in sorted(rows, key=lambda r: -r[2]):
    print(f"{q}\t{ln}\t{frac:.3f}\t{ident:.4f}\t{tg}")

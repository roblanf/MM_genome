#!/usr/bin/env python3
"""Cross-haplotype duplicate removal, after the parental split and before scaffolding.

A fragment (any contig that is not a chromosome-scale anchor) is a cross-haplotype duplicate if
its best alignment to another contig of the union is near-identical (>= --min-identity) and covers
>= --min-cov of it. Such a fragment is the same molecule as part of another contig. Read depth in
its own haplotype's BAM then decides: at single-parent depth (<= --max-depth) it is the other
parent's copy and is DROPPED; above that, both parents' reads support it, it is shared sequence
the chromosome contig lacks, and it is KEPT. Fragments at parental divergence from everything
(own-sequence) and fragments with no alignment (repeat-only) are kept.

Usage:
  python3 crosshap_dedup.py --union union.fa --anchors parent_of_chromosome.tsv \
      --bam-h1 03_qc/hap1.bam --bam-h2 03_qc/hap2.bam [--threads 32] [--max-depth 41.3] --out crosshap
Writes <out>.tsv (every fragment: best hit, identity, covered, depth, class, decision),
<out>.drop.txt, <out>.keep.txt. Alignment and depth are cached as <out>.paf / <out>.regions.bed.gz.
"""
import argparse, os, subprocess, sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--union", required=True, help="FASTA of both parental assemblies (all contigs)")
ap.add_argument("--anchors", required=True, help="list of chromosome-scale contigs (first column), excluded from the fragment pool")
ap.add_argument("--bam-h1", required=True, help="reads mapped to hap1 alone")
ap.add_argument("--bam-h2", required=True, help="reads mapped to hap2 alone")
ap.add_argument("--threads", type=int, default=32)
ap.add_argument("--min-aln", type=int, default=10000)
ap.add_argument("--min-identity", type=float, default=0.995)
ap.add_argument("--min-cov", type=float, default=0.5)
ap.add_argument("--max-depth", type=float, default=41.3, help="depth at or below which a duplicate is a single-parent copy")
ap.add_argument("--out", default="crosshap")
a = ap.parse_args()


def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"command failed: {cmd}\n{r.stderr[-2000:]}")


anchors = {l.split()[0] for l in open(a.anchors) if l.strip() and not l.startswith("#")}
if not os.path.exists(a.union + ".fai"):
    run(f"samtools faidx {a.union}")
length = {l.split()[0]: int(l.split()[1]) for l in open(a.union + ".fai")}
pool = [c for c in length if c not in anchors]
with open(a.out + ".pool.txt", "w") as f:
    f.write("\n".join(pool) + "\n")

# ---- alignment of the pool to the union -------------------------------------------------------
paf = a.out + ".paf"
if not os.path.exists(paf):
    print("aligning fragments to the union ...", file=sys.stderr)
    run(f"seqkit grep -f {a.out}.pool.txt {a.union} > {a.out}.pool.fa")
    run(f"minimap2 -cx asm20 -t {a.threads} {a.union} {a.out}.pool.fa > {paf} 2> {paf}.log")
best, cov = {}, defaultdict(lambda: defaultdict(list))
for line in open(paf):
    f = line.split("\t")
    q, t, qs, qe, nm, al = f[0], f[5], int(f[2]), int(f[3]), int(f[9]), int(f[10])
    if q == t or al < a.min_aln:
        continue
    de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
    ident = 1 - de if de is not None else nm / al
    cov[q][t].append((qs, qe))
    if q not in best or ident > best[q][1]:
        best[q] = (t, ident)


def covered(q, t):
    iv = sorted(cov[q][t]); tot, cs, ce = 0, iv[0][0], iv[0][1]
    for s, e in iv[1:]:
        if s > ce:
            tot += ce - cs; cs, ce = s, e
        else:
            ce = max(ce, e)
    return (tot + ce - cs) / length[q]


cls = {}
for q in pool:
    if q not in best:
        cls[q] = ("-", 0.0, 0.0, "no-hit")
    else:
        t, ident = best[q]; cv = covered(q, t)
        cls[q] = (t, ident, cv, "duplicate" if (ident >= a.min_identity and cv >= a.min_cov) else "own-sequence")

# ---- depth of the duplicates, in their own haplotype's BAM --------------------------------------
depth = {}
reg = a.out + ".regions.bed.gz"
if not os.path.exists(reg):
    for hap, bam in (("h1", a.bam_h1), ("h2", a.bam_h2)):
        dups = [q for q in pool if cls[q][3] == "duplicate" and q.startswith(hap)]
        if not dups:
            continue
        with open(f"{a.out}.{hap}.bed", "w") as f:
            for q in dups:
                f.write(f"{q}\t0\t{length[q]}\n")
        run(f"mosdepth -t {min(a.threads, 16)} --by {a.out}.{hap}.bed --no-per-base {a.out}.{hap} {bam}")
    run(f"zcat {a.out}.h1.regions.bed.gz {a.out}.h2.regions.bed.gz 2>/dev/null | gzip > {reg}")
import gzip
for line in gzip.open(reg, "rt"):
    f = line.split()
    depth[f[0]] = float(f[3])

# ---- decisions ---------------------------------------------------------------------------------
# A duplicate at single-parent depth is dropped - unless its best hit is itself a fragment that is
# being dropped (two copies of the same region, no chromosome copy): then this one is kept, so that
# a region never disappears entirely. Deciding shortest-first keeps the longer copy of such a pair.
decision = {}
for q in sorted(pool, key=lambda x: length[x]):
    t, ident, cv, c = cls[q]; d = depth.get(q)
    if c != "duplicate":
        decision[q] = "keep"
    elif d is None or d > a.max_depth:
        decision[q] = "keep (shared sequence)"
    elif t in anchors or decision.get(t, "keep").startswith("keep"):
        decision[q] = "drop (single-parent copy)"
    else:
        decision[q] = "keep (only remaining copy)"
rows, drop, keep = [], [], []
for q in sorted(pool, key=lambda x: -length[x]):
    t, ident, cv, c = cls[q]; d = depth.get(q)
    (drop if decision[q].startswith("drop") else keep).append(q)
    rows.append((q, length[q], t, f"{ident:.4f}", f"{cv:.2f}", "-" if d is None else f"{d:.1f}", c, decision[q]))
with open(a.out + ".tsv", "w") as f:
    f.write("fragment\tlength\tbest_hit\tidentity\tcovered\tdepth\tclass\tdecision\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
open(a.out + ".drop.txt", "w").write("".join(q + "\n" for q in drop))
open(a.out + ".keep.txt", "w").write("".join(q + "\n" for q in keep))

summ = defaultdict(lambda: [0, 0])
for r in rows:
    summ[(r[6], r[7])][0] += 1; summ[(r[6], r[7])][1] += r[1]
for (c, dec), (n, L) in sorted(summ.items()):
    print(f"{c:13s} {dec:26s} n={n:4d}  {L/1e6:6.2f} Mb")
print(f"drop: {len(drop)} fragments ({sum(length[q] for q in drop)/1e6:.2f} Mb) -> {a.out}.drop.txt;  keep: {len(keep)} -> {a.out}.keep.txt")

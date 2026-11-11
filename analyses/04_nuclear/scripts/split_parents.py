#!/usr/bin/env python3
"""Assign every contig of hap1 and hap2 to one of two parental genomes.

Rules, in priority order:
  chromosome : contig listed in --chrom (chromosome-scale contigs) -> parent given there
  group      : fragment whose best partner in the other haplotype is a listed chromosome contig
               covering >= --min-cov of it -> the OPPOSITE parent of that partner (the fragment
               belongs to its own haplotype's copy of that chromosome, which has a gap there)
  intrinsic  : ungrouped fragment -> its own evidence: Δ ID call from dv_scoring.tsv and/or
               hapmer majority >= --hapmer-majority from hapmers.count; if both exist they must agree
  unassigned : nothing confident

Also reports, for group fragments, whether their own evidence is consistent:
  consistent      : intrinsic call == group parent
  identical-filler: >= 99% identical to the partner chromosome (both parents carry it); resembling
                    the partner's parent is expected and not a conflict
  CHECK           : < 99% identical to the partner and intrinsic call == partner's parent -> possible
                    misphased fragment; verify with read depth (half depth = other parent's allele)

Usage:
  python3 split_parents.py --pairs both.pairs.tsv --chrom parent_of_chromosome.tsv \
      --fai hap1.purged.fa.fai hap2.purged.fa.fai [--dv dv_scoring.tsv] [--hapmers hapmers.count] \
      [--min-cov 0.5] [--hapmer-majority 0.9] --out parent_of_contig.tsv
Writes <out>, plus virginea.txt / decipiens.txt / unassigned.txt next to it.
"""
import argparse, csv, os
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--pairs", required=True)
ap.add_argument("--chrom", required=True, help="TSV: contig <tab> virginea|decipiens")
ap.add_argument("--fai", nargs="+", required=True, help="samtools faidx .fai of hap1 and hap2 (or the FASTA files)")
ap.add_argument("--dv"); ap.add_argument("--hapmers")
ap.add_argument("--min-cov", type=float, default=0.5)
ap.add_argument("--hapmer-majority", type=float, default=0.9)
ap.add_argument("--min-hapmers", type=int, default=1000)
ap.add_argument("--out", default="parent_of_contig.tsv")
a = ap.parse_args()

OPP = {"virginea": "decipiens", "decipiens": "virginea"}
parent_chrom = {}
for line in open(a.chrom):
    if line.strip() and not line.startswith("#"):
        c, p = line.split()[:2]
        if p not in OPP:
            raise SystemExit(f"parent must be virginea or decipiens: {line}")
        parent_chrom[c] = p

# contig lengths
length = {}
for f in a.fai:
    if f.endswith(".fai"):
        for line in open(f):
            x = line.split("\t"); length[x[0]] = int(x[1])
    else:
        name, n = None, 0
        for line in open(f):
            if line.startswith(">"):
                if name: length[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
        if name: length[name] = n

# best partner of every contig from the pairs table
best = {}   # contig -> (partner, cov, identity)
with open(a.pairs) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        for me, other, cov in ((r["h1"], r["h2"], float(r["h1_cov"])), (r["h2"], r["h1"], float(r["h2_cov"]))):
            if me not in best or cov > best[me][1]:
                best[me] = (other, cov, float(r["identity"]))

# intrinsic evidence
dv_call = {}
if a.dv:
    with open(a.dv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Call"] in OPP:
                dv_call[r["Contig"]] = r["Call"]
hap_call = {}
if a.hapmers:
    with open(a.hapmers) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v = int(r[next(k for k in r if "virginea" in k)]); d = int(r[next(k for k in r if "decipiens" in k)])
            if v + d >= a.min_hapmers:
                if v / (v + d) >= a.hapmer_majority: hap_call[r["Contig"]] = "virginea"
                elif d / (v + d) >= a.hapmer_majority: hap_call[r["Contig"]] = "decipiens"


def intrinsic(c):
    calls = {x for x in (dv_call.get(c), hap_call.get(c)) if x}
    return calls.pop() if len(calls) == 1 else None


rows = []
for c in sorted(length, key=lambda x: (x[:2], -length[x])):
    hap = "hap1" if c.startswith("h1") else "hap2"
    partner, cov, ident = best.get(c, ("-", 0.0, 0.0))
    intr = intrinsic(c)
    if c in parent_chrom:
        rule, parent, check = "chromosome", parent_chrom[c], "-"
    elif partner in parent_chrom and cov >= a.min_cov:
        rule, parent = "group", OPP[parent_chrom[partner]]
        if intr is None: check = "no intrinsic call"
        elif intr == parent: check = "consistent"
        elif ident >= 0.99: check = "identical-filler"
        else: check = "CHECK"
    elif intr:
        rule, parent, check = "intrinsic", intr, "-"
    else:
        rule, parent, check = "unassigned", "unassigned", "-"
    rows.append((c, hap, length[c], rule, parent, partner, f"{cov:.2f}", f"{ident:.3f}",
                 dv_call.get(c, "-"), hap_call.get(c, "-"), check))

outdir = os.path.dirname(os.path.abspath(a.out))
with open(a.out, "w") as out:
    out.write("contig\thap\tlength\trule\tparent\tpartner\tpartner_cov\tpartner_identity\tdv_call\thapmer_call\tcheck\n")
    for r in rows:
        out.write("\t".join(map(str, r)) + "\n")
for p in ("virginea", "decipiens", "unassigned"):
    with open(os.path.join(outdir, f"{p}.txt"), "w") as f:
        for r in rows:
            if r[4] == p: f.write(r[0] + "\n")

summ = defaultdict(lambda: [0, 0])
for r in rows:
    summ[(r[4], r[3])][0] += 1; summ[(r[4], r[3])][1] += r[2]
for (p, rule), (n, L) in sorted(summ.items()):
    print(f"{p:11s} {rule:11s} n={n:4d}  {L/1e6:7.2f} Mb")
n_check = sum(1 for r in rows if r[10] == "CHECK")
print(f"fragments flagged CHECK: {n_check}")

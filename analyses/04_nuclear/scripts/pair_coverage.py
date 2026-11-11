#!/usr/bin/env python3
"""Pairwise coverage between the contigs of two haplotypes from a minimap2 PAF.

Input: PAF of hap1 (query) aligned to hap2 (target), e.g.
    minimap2 -cx asm20 --secondary=no -t 128 hap2.purged.fa hap1.purged.fa > h1_vs_h2.paf

Outputs (TSV):
  <prefix>.pairs.tsv    one row per (hap1 contig, hap2 contig) pair:
        h1, h1_len, h2, h2_len, n_aln, aligned_bp, identity,
        h1_cov (fraction of the hap1 contig covered by alignments to this hap2 contig),
        h2_cov (fraction of the hap2 contig covered by alignments to this hap1 contig),
        best_for_h1, best_for_h2, reciprocal
  <prefix>.contigs.tsv  one row per contig of either haplotype:
        contig, hap, length, total_cov (fraction covered by the other haplotype, all partners),
        best_partner, best_partner_cov, reciprocal, partners (partner:cov for every partner >= --min-cov)

Usage:
    python3 pair_coverage.py h1_vs_h2.paf --prefix h1_vs_h2 [--min-aln 10000] [--min-id 0.90] [--min-cov 0.01]
"""
import argparse
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("paf")
ap.add_argument("--prefix", default="pairs")
ap.add_argument("--min-aln", type=int, default=10000, help="min alignment block length in bp (default 10000)")
ap.add_argument("--min-id", type=float, default=0.90, help="min identity (1 - de, or matches/alen) (default 0.90)")
ap.add_argument("--min-cov", type=float, default=0.01, help="report pairs covering >= this fraction of either contig (default 0.01)")
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


q_iv = defaultdict(list)      # (q,t) -> intervals on q
t_iv = defaultdict(list)      # (q,t) -> intervals on t
q_all = defaultdict(list)     # q -> intervals on q (all targets)
t_all = defaultdict(list)     # t -> intervals on t (all queries)
n_aln = defaultdict(int)
aln_bp = defaultdict(int)
match_bp = defaultdict(float)
qlen, tlen = {}, {}

with open(a.paf) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, ql, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
        t, tl, ts, te = f[5], int(f[6]), int(f[7]), int(f[8])
        nmatch, alen = int(f[9]), int(f[10])
        if alen < a.min_aln:
            continue
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        ident = 1.0 - de if de is not None else nmatch / alen
        if ident < a.min_id:
            continue
        qlen[q], tlen[t] = ql, tl
        key = (q, t)
        q_iv[key].append((qs, qe)); t_iv[key].append((ts, te))
        q_all[q].append((qs, qe));  t_all[t].append((ts, te))
        n_aln[key] += 1; aln_bp[key] += alen; match_bp[key] += ident * alen

pairs = {}
for key in q_iv:
    q, t = key
    pairs[key] = dict(h1_cov=merged_len(q_iv[key]) / qlen[q], h2_cov=merged_len(t_iv[key]) / tlen[t],
                      identity=match_bp[key] / aln_bp[key])

best_q = {}   # q -> t with highest h1_cov
best_t = {}   # t -> q with highest h2_cov
for (q, t), p in pairs.items():
    if q not in best_q or p["h1_cov"] > pairs[(q, best_q[q])]["h1_cov"]:
        best_q[q] = t
    if t not in best_t or p["h2_cov"] > pairs[(best_t[t], t)]["h2_cov"]:
        best_t[t] = q

with open(f"{a.prefix}.pairs.tsv", "w") as out:
    out.write("h1\th1_len\th2\th2_len\tn_aln\taligned_bp\tidentity\th1_cov\th2_cov\tbest_for_h1\tbest_for_h2\treciprocal\n")
    for (q, t), p in sorted(pairs.items(), key=lambda kv: -aln_bp[kv[0]]):
        if p["h1_cov"] < a.min_cov and p["h2_cov"] < a.min_cov:
            continue
        b1, b2 = best_q[q] == t, best_t[t] == q
        out.write(f"{q}\t{qlen[q]}\t{t}\t{tlen[t]}\t{n_aln[(q,t)]}\t{aln_bp[(q,t)]}\t{p['identity']:.4f}\t"
                  f"{p['h1_cov']:.3f}\t{p['h2_cov']:.3f}\t{'yes' if b1 else 'no'}\t{'yes' if b2 else 'no'}\t"
                  f"{'yes' if (b1 and b2) else 'no'}\n")

with open(f"{a.prefix}.contigs.tsv", "w") as out:
    out.write("contig\thap\tlength\ttotal_cov\tbest_partner\tbest_partner_cov\treciprocal\tpartners\n")
    rows = []
    for q in qlen:
        parts = sorted(((t, pairs[(q, t)]["h1_cov"]) for (qq, t) in pairs if qq == q), key=lambda x: -x[1])
        bt = best_q[q]
        rows.append((q, "hap1", qlen[q], merged_len(q_all[q]) / qlen[q], bt, pairs[(q, bt)]["h1_cov"],
                     best_t[bt] == q, parts))
    for t in tlen:
        parts = sorted(((q, pairs[(q, t)]["h2_cov"]) for (q, tt) in pairs if tt == t), key=lambda x: -x[1])
        bq = best_t[t]
        rows.append((t, "hap2", tlen[t], merged_len(t_all[t]) / tlen[t], bq, pairs[(bq, t)]["h2_cov"],
                     best_q[bq] == t, parts))
    for c, h, L, tc, bp, bpc, rec, parts in sorted(rows, key=lambda r: (r[1], -r[2])):
        plist = ",".join(f"{p}:{v:.2f}" for p, v in parts if v >= a.min_cov)
        out.write(f"{c}\t{h}\t{L}\t{tc:.3f}\t{bp}\t{bpc:.3f}\t{'yes' if rec else 'no'}\t{plist}\n")

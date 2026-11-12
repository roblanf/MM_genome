#!/usr/bin/env python3
"""Fragment audit after per-chromosome scaffolding.

For ONE parental genome (the WORK directory written by scaffold_per_chromosome.sh), using the
union of both haplotypes as the comparison set:

  1. every non-anchor contig of the union (the fragment pool, both parents) is aligned to the
     union and classified by its best non-self hit: 'duplicate' (>= --min-identity, i.e. the same
     molecule as a region of another contig) or 'own-sequence' (a distinct allele / unique), or
     'no-hit' (repeat-only)
  2. telomere arrays are searched in all fragments (tidk) and in the anchors of this parent
  3. for this parent: which fragments were placed in chromosome scaffolds, which lie beyond a
     telomere (from telomere_check), which are duplicates -> lists to exclude / drop
  4. for every anchor end of this parent WITHOUT a telomere: telomere-bearing fragments (from
     either parent) whose best hit is that anchor's end (same molecule) or the homologous
     chromosome's end (own allele) -> candidates to attach

Outputs in WORK/audit/: fragments.tsv, ends.tsv, drop_duplicates.txt, exclude_beyond_telomere.txt,
attach_candidates.tsv. External tools: seqkit, minimap2, tidk (results cached; rerun with --paf /
--tidk to reuse).

Usage:
  python3 fragment_audit.py --work 07c_scaffold_correct_bin/vir --union union.purged.fa \
      --chrom parent_of_chromosome.tsv --parent-of-contig parent_of_contig.tsv --label virginea [--threads 32]
"""
import argparse, os, subprocess, sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--work", required=True, help="WORK directory of scaffold_per_chromosome.sh for this parent")
ap.add_argument("--union", required=True, help="FASTA of both haplotypes (all contigs)")
ap.add_argument("--chrom", required=True, help="parent_of_chromosome.tsv (anchor <tab> parent)")
ap.add_argument("--parent-of-contig", required=True, help="parent_of_contig.tsv")
ap.add_argument("--label", required=True, help="this parent's label in the tables: virginea | decipiens")
ap.add_argument("--threads", type=int, default=32)
ap.add_argument("--min-aln", type=int, default=10000)
ap.add_argument("--min-identity", type=float, default=0.995)
ap.add_argument("--min-cov", type=float, default=0.5, help="fraction of a fragment covered by its best-hit contig to call it a duplicate")
ap.add_argument("--end-window", type=int, default=300000, help="a hit within this distance of an anchor end counts as 'at the end'")
ap.add_argument("--telomere-copies", type=int, default=500, help="repeat copies for a fragment/end to count as telomere-bearing")
ap.add_argument("--edge", type=int, default=100000)
a = ap.parse_args()

W = os.path.abspath(a.work); OUT = os.path.join(W, "audit"); os.makedirs(OUT, exist_ok=True)
def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"command failed: {cmd}\n{r.stderr[-2000:]}")
    return r.stdout

# ---- anchors, parents, pairs ----------------------------------------------------------------
anchor_parent = {l.split()[0]: l.split()[1] for l in open(a.chrom) if l.strip() and not l.startswith("#")}
anchors = set(anchor_parent)
parent_of = {}
with open(a.parent_of_contig) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ic, ip = hdr.index("contig"), hdr.index("parent")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) > max(ic, ip):
            parent_of[f[ic]] = f[ip]
pairs = {}
pf = os.path.join(W, "pairs.tsv")
if os.path.exists(pf):
    for line in open(pf):
        h1, h2 = line.split()[:2]
        pairs[h1] = h2; pairs.setdefault(h2, h1)        # homologue lookup both ways (first pairing wins for h2)
        if h2 not in pairs or pairs[h2] == h1: pairs[h2] = h1

# ---- lengths and the fragment pool ------------------------------------------------------------
fai = a.union + ".fai"
if not os.path.exists(fai):
    run(f"samtools faidx {a.union}")
length = {l.split()[0]: int(l.split()[1]) for l in open(fai)}
pool = [c for c in length if c not in anchors]
with open(os.path.join(OUT, "pool.txt"), "w") as f:
    f.write("\n".join(pool) + "\n")
pool_fa = os.path.join(OUT, "pool.fa")
if not os.path.exists(pool_fa):
    run(f"seqkit grep -f {OUT}/pool.txt {a.union} > {pool_fa}")

# ---- what this parent's scaffolding did with its fragments -------------------------------------
placed = {}            # fragment -> scaffold
for line in open(os.path.join(W, "combined.agp")):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) >= 9 and f[4] == "W" and f[5] not in anchors:
        placed[f[5]] = f[0]
beyond = set()
bt = os.path.join(W, "beyond_telomere.txt")
if os.path.exists(bt):
    beyond = {l.strip() for l in open(bt) if l.strip()}

# ---- alignment of the pool to the union (cached) --------------------------------------------
paf = os.path.join(OUT, "pool_vs_union.paf")
if not os.path.exists(paf):
    print("aligning fragment pool to the union ...", file=sys.stderr)
    run(f"minimap2 -cx asm20 -t {a.threads} {a.union} {pool_fa} > {paf} 2> {paf}.log")
best = {}              # fragment -> (hit, identity, hit_start, hit_end)
cov = defaultdict(lambda: defaultdict(list))   # fragment -> hit -> [query intervals]
for line in open(paf):
    f = line.split("\t")
    q, t, qs, qe, ts, te, nm, al = f[0], f[5], int(f[2]), int(f[3]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
    if q == t or al < a.min_aln: continue
    de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
    ident = 1 - de if de is not None else nm / al
    cov[q][t].append((qs, qe))
    if q not in best or ident > best[q][1] or (ident == best[q][1] and al > (best[q][3] - best[q][2])):
        best[q] = (t, ident, ts, te)
def covered(q, t):
    iv = sorted(cov[q][t]); tot, cs, ce = 0, iv[0][0], iv[0][1]
    for s, e in iv[1:]:
        if s > ce: tot += ce - cs; cs, ce = s, e
        else: ce = max(ce, e)
    return (tot + ce - cs) / length[q]

# ---- telomeres in fragments and anchors (cached) ------------------------------------------------
def tidk(fa, name):
    tsv = os.path.join(OUT, name, f"{name}_telomeric_repeat_windows.tsv")
    if not os.path.exists(tsv):
        run(f"tidk search -s AAACCCT -o {name} -d {OUT}/{name} {fa} > /dev/null 2>&1")
    win = defaultdict(list)
    with open(tsv) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4: win[f[0]].append((int(f[1]), int(f[2]) + int(f[3])))
    return win
pool_win = tidk(pool_fa, "tidk_pool")
telo = {c: sum(n for _, n in ws if n >= 100) for c, ws in pool_win.items()}
def telo_side(c):
    """where in the fragment the array sits: 'start' / 'end' / 'internal' / '-' """
    ws = [(w, n) for w, n in pool_win.get(c, []) if n >= 100]
    if not ws: return "-"
    L = length[c]; first, last = min(w for w, _ in ws), max(w for w, _ in ws)
    if first <= 20000 + 10000: return "start"
    if last >= L - 20000: return "end"
    return "internal"
my_anchors = [c for c in anchors if anchor_parent[c] == a.label]
anc_fa = os.path.join(OUT, "anchors.fa")
if not os.path.exists(anc_fa):
    with open(os.path.join(OUT, "anchors.txt"), "w") as f: f.write("\n".join(my_anchors) + "\n")
    run(f"seqkit grep -f {OUT}/anchors.txt {a.union} > {anc_fa}")
anc_win = tidk(anc_fa, "tidk_anchors")
def anchor_end_has_telomere(c, side):
    ws = anc_win.get(c, [])
    if side == "left":  return any(n >= 100 and w <= a.edge for w, n in ws)
    return any(n >= 100 and w >= length[c] - a.edge for w, n in ws)

# ---- fragments table ----------------------------------------------------------------------------
rows = []
for c in sorted(pool, key=lambda x: -length[x]):
    hit, ident, ts, te = best.get(c, ("-", 0.0, 0, 0))
    cv = covered(c, hit) if hit != "-" else 0.0
    cls = "no-hit" if hit == "-" else ("duplicate" if (ident >= a.min_identity and cv >= a.min_cov) else "own-sequence")
    hit_pos = "-"
    if hit in anchors:
        hit_pos = "left-end" if ts <= a.end_window else ("right-end" if te >= length[hit] - a.end_window else "internal")
    rows.append([c, parent_of.get(c, "-"), length[c], placed.get(c, "unplaced" if parent_of.get(c) == a.label else "-"),
                 "yes" if c in beyond else "no", hit, f"{ident:.4f}", f"{cv:.2f}", hit_pos, cls, telo.get(c, 0), telo_side(c)])
with open(os.path.join(OUT, "fragments.tsv"), "w") as f:
    f.write("fragment\tparent_file\tlength\tplaced_in\tbeyond_telomere\tbest_hit\tidentity\tcovered\thit_position\tclass\ttelomere_copies\ttelomere_side\n")
    for r in rows: f.write("\t".join(map(str, r)) + "\n")
mine = [r for r in rows if r[1] == a.label]
with open(os.path.join(OUT, "drop_duplicates.txt"), "w") as f:
    f.write("".join(r[0] + "\n" for r in mine if r[9] == "duplicate"))
with open(os.path.join(OUT, "exclude_beyond_telomere.txt"), "w") as f:
    f.write("".join(r[0] + "\n" for r in mine if r[4] == "yes"))

# ---- ends without a telomere: candidates --------------------------------------------------------
ends, cands = [], []
telo_frags = [r for r in rows if r[10] >= a.telomere_copies]
for c in sorted(my_anchors, key=lambda x: -length[x]):
    for side in ("left", "right"):
        has = anchor_end_has_telomere(c, side)
        found = []
        if not has:
            for r in telo_frags:
                frag, hit, pos, ident = r[0], r[5], r[8], float(r[6])
                if hit == c and pos == f"{side}-end":
                    found.append(f"{frag} (same molecule, id {ident:.3f}, from {r[1]} file)")
                elif hit == pairs.get(c) and pos not in ("-", "internal") and ident < a.min_identity:   # own allele: diverged from the homologue
                    found.append(f"{frag} (own allele, maps to homologue {hit} {pos}, id {ident:.3f}, from {r[1]} file)")
        ends.append((c, side, "yes" if has else "no", "; ".join(found) or ("-" if has else "none found")))
        for x in found: cands.append((c, side, x))
with open(os.path.join(OUT, "ends.tsv"), "w") as f:
    f.write("anchor\tside\ttelomere\tcandidate_fragments\n")
    for e in ends: f.write("\t".join(e) + "\n")
with open(os.path.join(OUT, "attach_candidates.tsv"), "w") as f:
    f.write("anchor\tside\tcandidate\n")
    for e in cands: f.write("\t".join(e) + "\n")

# ---- summary --------------------------------------------------------------------------------------
n_cls = defaultdict(lambda: [0, 0])
for r in mine: n_cls[r[9]][0] += 1; n_cls[r[9]][1] += r[2]
print(f"fragments in the {a.label} file: {len(mine)}")
for k, (n, L) in sorted(n_cls.items()): print(f"  {k:13s} {n:4d}  {L/1e6:6.2f} Mb")
print(f"  beyond a telomere: {sum(1 for r in mine if r[4]=='yes')}   placed in chromosomes: {sum(1 for r in mine if r[3] not in ('unplaced','-'))}")
print(f"  telomere-bearing (>= {a.telomere_copies} copies): {sum(1 for r in mine if r[10] >= a.telomere_copies)}")
print(f"anchor ends with a telomere: {sum(1 for e in ends if e[2]=='yes')}/{len(ends)}; ends with attachment candidates: {len({(e[0],e[1]) for e in cands})}")
for e in ends:
    if e[2] == "no": print(f"  {e[0]} {e[1]}: {e[3]}")
print(f"written to {OUT}: fragments.tsv, ends.tsv, attach_candidates.tsv, drop_duplicates.txt, exclude_beyond_telomere.txt")

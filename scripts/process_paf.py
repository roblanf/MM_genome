import os
import re

def parse_paf(paf_path):
    if not os.path.exists(paf_path):
        return {}
    
    stats = {}
    with open(paf_path, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 12: continue
            
            q_name = cols[0]
            q_len = int(cols[1])
            q_start = int(cols[2])
            q_end = int(cols[3])
            t_name = cols[5]
            matches = int(cols[9])
            block_len = int(cols[10])
            
            # Extract divergence (dv:f:) tag
            dv_match = re.search(r'dv:f:([\d.]+)', line)
            dv = float(dv_match.group(1)) if dv_match else 0.05
            
            if q_name not in stats:
                stats[q_name] = {
                    'q_len': q_len,
                    'aligned_bases': 0,
                    'weighted_id': 0.0,
                    'match_sum': 0,
                    'best_target': t_name
                }
            
            # Sum up stats
            stats[q_name]['aligned_bases'] += (q_end - q_start)
            stats[q_name]['match_sum'] += matches
            stats[q_name]['weighted_id'] += (1 - dv) * matches

    return stats

# File paths
files = {
    "H1_Dec": "04_parental_assignment/hap1_vs_decipiens.paf",
    "H1_Vir": "04_parental_assignment/hap1_vs_virginea.paf",
    "H2_Dec": "04_parental_assignment/hap2_vs_decipiens.paf",
    "H2_Vir": "04_parental_assignment/hap2_vs_virginea.paf"
}

# Load data
data = {k: parse_paf(v) for k, v in files.items()}
all_chroms_h1 = sorted(data["H1_Dec"].keys())
all_chroms_h2 = sorted(data["H2_Dec"].keys())

print("## Parental Assignment Summary (E. phylacis Haplotypes)")
print("\n### Haplotype 1 Assignment")
print("| Hybrid Chrom | Length (Mb) | Match E. decipiens (%) | Match E. virginea (%) | Assignment | Coverage |")
print("|:---|:---:|:---:|:---:|:---:|:---:|")

for c in all_chroms_h1:
    d = data["H1_Dec"].get(c, {'id':0, 'weighted_id':0, 'match_sum':1, 'aligned_bases':0, 'q_len':1})
    v = data["H1_Vir"].get(c, {'id':0, 'weighted_id':0, 'match_sum':1, 'aligned_bases':0, 'q_len':1})
    
    id_d = (d['weighted_id']/d['match_sum'])*100 if d['match_sum']>0 else 0
    id_v = (v['weighted_id']/v['match_sum'])*100 if v['match_sum']>0 else 0
    cov = (max(d['aligned_bases'], v['aligned_bases']) / d['q_len']) * 100
    winner = "**E. decipiens**" if id_d > id_v else "**E. virginea**"
    
    print(f"| {c} | {d['q_len']/1e6:.1f} | {id_d:.2f}% | {id_v:.2f}% | {winner} | {cov:.1f}% |")

print("\n### Haplotype 2 Assignment")
print("| Hybrid Chrom | Length (Mb) | Match E. decipiens (%) | Match E. virginea (%) | Assignment | Coverage |")
print("|:---|:---:|:---:|:---:|:---:|:---:|")

for c in all_chroms_h2:
    d = data["H2_Dec"].get(c, {'id':0, 'weighted_id':0, 'match_sum':1, 'aligned_bases':0, 'q_len':1})
    v = data["H2_Vir"].get(c, {'id':0, 'weighted_id':0, 'match_sum':1, 'aligned_bases':0, 'q_len':1})
    
    id_d = (d['weighted_id']/d['match_sum'])*100 if d['match_sum']>0 else 0
    id_v = (v['weighted_id']/v['match_sum'])*100 if v['match_sum']>0 else 0
    cov = (max(d['aligned_bases'], v['aligned_bases']) / d['q_len']) * 100
    winner = "**E. decipiens**" if id_d > id_v else "**E. virginea**"
    
    print(f"| {c} | {d['q_len']/1e6:.1f} | {id_d:.2f}% | {id_v:.2f}% | {winner} | {cov:.1f}% |")

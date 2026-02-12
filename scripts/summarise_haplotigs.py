import re
import os

# Configuration
input_dir = os.path.expanduser("~/Desktop/euc_contigs_ramdisk/") 
haplotig_pattern = r"h[12]tg\d+\w*" 

def get_gene_count(root_path):
    # ROADIES places num_gt.txt inside the statistics folder 
    gt_file = os.path.join(root_path, "statistics", "num_gt.txt")
    if os.path.exists(gt_file):
        with open(gt_file, 'r') as f:
            content = f.read().strip()
            # Grabs the integer count at the end of the line 
            match = re.search(r"(\d+)$", content)
            if match:
                return int(match.group(1))
    return 0

def parse_astral_node(newick_str, target_id):
    # Captures Target BLen, Sister Name, and Metadata Block
    # Case 1: (Target:BLen,Sister:Len)'[meta]'
    pattern = rf"\(({target_id}:([^,)]+),([^:]+):[^)]+)\)'\[(.*?)\]'"
    match = re.search(pattern, newick_str)
    
    if not match:
        # Case 2: (Sister:Len,Target:BLen)'[meta]'
        pattern = rf"\((([^:]+):[^,]+,{target_id}:([^)]+))\)'\[(.*?)\]'"
        match = re.search(pattern, newick_str)
        if match:
            sister = match.group(2)
            branch_len = match.group(3)
            meta_str = match.group(4)
        else:
            return None, None, None
    else:
        branch_len = match.group(2)
        sister = match.group(3)
        meta_str = match.group(4)
    
    try:
        # Extracts f1-f3 and q1-q3 from the ASTRAL metadata string
        meta = dict(item.split("=") for item in meta_str.split(";") if "=" in item)
        return sister, branch_len, meta
    except ValueError:
        return sister, branch_len, {"error": "meta_parse_fail"}

results = []

# Gather data from all directories
for root, dirs, files in os.walk(input_dir):
    for filename in files:
        if filename == "roadies_stats.nwk":
            filepath = os.path.join(root, filename)
            with open(filepath, 'r') as f:
                tree = f.read()
            
            hid_match = re.search(haplotig_pattern, tree)
            if hid_match:
                hid = hid_match.group(0)
                sister, blen, stats = parse_astral_node(tree, hid)
                num_genes = get_gene_count(root)
                
                if stats:
                    results.append({
                        'hid': hid,
                        'sister': sister.strip('('),
                        'blen': blen,
                        'genes': num_genes,
                        'stats': stats
                    })

# Natural sort (h1tg1, h1tg2, h1tg10, h2tg1...)
results.sort(key=lambda x: [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x['hid'])])

# Formatting Headers
header = f"{'Haplotig':<15} {'Sister_Group':<20} {'Genes':<6} {'BLen':<10} {'pp1':<6} {'f1':<8} {'f2':<8} {'f3':<8} {'q1':<8} {'q2':<8} {'q3':<8}"
print(header)
print("-" * len(header))

for r in results:
    s = r['stats']
    print(f"{r['hid']:<15} {r['sister']:<20} "
          f"{r['genes']:<6} "
          f"{r['blen']:<10.8} " # Limit BLen display to 8 chars
          f"{s.get('pp1','-'):<6} "
          f"{s.get('f1','-'):<8.6} " # f counts as strings can be long
          f"{s.get('f2','-'):<8.6} "
          f"{s.get('f3','-'):<8.6} "
          f"{s.get('q1','-'):<8} "
          f"{s.get('q2','-'):<8} "
          f"{s.get('q3','-'):<8}")

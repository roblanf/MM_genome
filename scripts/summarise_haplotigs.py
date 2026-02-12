# Call this to get a little summary table of each haplotig analysis from ROADIES
# It should tell you the sister group; number of genes used for the analysis, and the stats for the branch that unites the haplotig with the sister group
import re
import os

# Configuration
input_dir = os.path.expanduser("~/Desktop/euc_contigs_ramdisk/") 
haplotig_pattern = r"h[12]tg\d+\w*" 

def get_gene_count(root_path):
    gt_file = os.path.join(root_path, "statistics", "num_gt.txt")
    if os.path.exists(gt_file):
        with open(gt_file, 'r') as f:
            content = f.read().strip()
            match = re.search(r"(\d+)$", content)
            if match:
                return int(match.group(1))
    return 0

def parse_astral_node(newick_str, target_id):
    pattern = rf"\(({target_id}:[^,]+,([^:]+):[^)]+)\)'\[(.*?)\]'"
    match = re.search(pattern, newick_str)
    if not match:
        pattern = rf"\((([^:]+):[^,]+,{target_id}:[^)]+)\)'\[(.*?)\]'"
        match = re.search(pattern, newick_str)
    
    if match:
        sister = match.group(2)
        meta_str = match.group(3)
        try:
            meta = dict(item.split("=") for item in meta_str.split(";") if "=" in item)
            return sister, meta
        except ValueError:
            return sister, {"error": "meta_parse_fail"}
    return None, None

results = []

# Collect all data first
for root, dirs, files in os.walk(input_dir):
    for filename in files:
        if filename == "roadies_stats.nwk":
            filepath = os.path.join(root, filename)
            with open(filepath, 'r') as f:
                tree = f.read()
            
            hid_match = re.search(haplotig_pattern, tree)
            if hid_match:
                hid = hid_match.group(0)
                sister, stats = parse_astral_node(tree, hid)
                num_genes = get_gene_count(root)
                
                if stats:
                    results.append({
                        'hid': hid,
                        'sister': sister.strip('('),
                        'genes': num_genes,
                        'stats': stats
                    })

# Sort results by haplotig name (hid)
# This handles h1tg... then h2tg... and numerical order correctly
results.sort(key=lambda x: [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x['hid'])])

# Print Headers
header = f"{'Haplotig':<15} {'Sister_Group':<25} {'Genes':<8} {'pp1':<8} {'f1':<10} {'f2':<10} {'f3':<10} {'q1':<10} {'q2':<10} {'q3':<10}"
print(header)
print("-" * len(header))

# Print Sorted Rows
for r in results:
    s = r['stats']
    print(f"{r['hid']:<15} {r['sister']:<25} "
          f"{r['genes']:<8} "
          f"{s.get('pp1','-'):<8} "
          f"{s.get('f1','-'):<10} "
          f"{s.get('f2','-'):<10} "
          f"{s.get('f3','-'):<10} "
          f"{s.get('q1','-'):<10} "
          f"{s.get('q2','-'):<10} "
          f"{s.get('q3','-'):<10}")

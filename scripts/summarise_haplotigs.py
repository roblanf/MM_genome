import os
import re
from ete3 import Tree

# --- Configuration ---
# Path to your contigs folder
INPUT_DIR = os.path.expanduser("~/Desktop/euc_contigs_ramdisk/")
# Outgroup for consistent rooting
OUTGROUP_NAME = "Angophora_floribunda"
# Regex to find the haplotig name in the tree
HAPLOTIG_PATTERN = r"h[12]tg\d+\w*"

def get_gene_count(root_path):
    """Reads the gene tree count from statistics/num_gt.txt."""
    gt_file = os.path.join(root_path, "statistics", "num_gt.txt")
    if os.path.exists(gt_file):
        with open(gt_file, 'r') as f:
            content = f.read().strip()
            match = re.search(r"(\d+)$", content)
            if match:
                return int(match.group(1))
    return 0

def parse_astral_stats(name_string):
    """Extracts pp1, f1, f2, f3, q1, q2, q3 from the ASTRAL comment string."""
    # ASTRAL/ROADIES puts stats in '[...]' which ETE3 loads as the node name
    clean_str = name_string.strip("'[] ")
    pairs = re.findall(r"(\w+)=([^;]+)", clean_str)
    return dict(pairs)

def get_sister_label(node):
    """Identifies the sister leaf or clade name."""
    sisters = node.get_sisters()
    if not sisters:
        return "No_Sister"
    
    sister = sisters[0]
    if sister.is_leaf():
        return sister.name
    else:
        # For a clade, return the first leaf name as a representative
        leaves = sister.get_leaf_names()
        return f"Clade({leaves[0]},...)" if len(leaves) > 1 else leaves[0]

def process_contig(root):
    tree_path = os.path.join(root, "roadies_stats.nwk")
    if not os.path.exists(tree_path):
        return None
        
    with open(tree_path, 'r') as f:
        newick_str = f.read().strip()
        
    # Determine the haplotig ID for this folder
    hid_match = re.search(HAPLOTIG_PATTERN, newick_str)
    if not hid_match:
        return None
    hid = hid_match.group(0)
    
    try:
        # format=1 reads internal node names (where ROADIES stores metadata)
        # quoted_node_names=True handles the ASTRAL '[...]' syntax
        t = Tree(newick_str, format=1, quoted_node_names=True)
        
        # 1. Root the tree
        outgroups = t.search_nodes(name=OUTGROUP_NAME)
        if outgroups:
            t.set_outgroup(outgroups[0])
        else:
            t.set_midpoint_outgroup()
            
        # 2. Find target haplotig
        target_nodes = t.search_nodes(name=hid)
        if not target_nodes:
            return None
        target = target_nodes[0]
        
        # 3. Collect Data
        # Terminal branch length
        blen = target.dist
        # Sister group in the rooted tree
        sister = get_sister_label(target)
        # Stats are located on the parent node (representing the split)
        stats = parse_astral_stats(target.up.name) if target.up else {}
        # Gene count
        genes = get_gene_count(root)
        
        return {
            'hid': hid,
            'sister': sister,
            'genes': genes,
            'blen': blen,
            'stats': stats
        }
    except Exception:
        return None

def main():
    results = []
    
    # Iterate through all subdirectories
    for root, dirs, files in os.walk(INPUT_DIR):
        if "roadies_stats.nwk" in files:
            res = process_contig(root)
            if res:
                results.append(res)
                
    # Natural sort by haplotig name (h1tg1 before h1tg10)
    results.sort(key=lambda x: [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x['hid'])])
    
    # Print Table Header
    header = ["Haplotig", "Sister_Group", "Genes", "BLen", "pp1", "f1", "f2", "f3", "q1", "q2", "q3"]
    row_fmt = "{:<18} {:<22} {:<6} {:<10.8} {:<6} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8}"
    
    print(row_fmt.format(*header))
    print("-" * 120)
    
    # Print Rows
    for r in results:
        s = r['stats']
        print(row_fmt.format(
            r['hid'], r['sister'], r['genes'], r['blen'],
            s.get('pp1', '-'), s.get('f1', '-'), s.get('f2', '-'), 
            s.get('f3', '-'), s.get('q1', '-'), s.get('q2', '-'), s.get('q3', '-')
        ))

if __name__ == "__main__":
    main()
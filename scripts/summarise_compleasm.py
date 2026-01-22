import glob
import os
import re

# Use an absolute path or ensure you are in the directory above '03_hifiasm_assembly'
results_path = "03_hifiasm_assembly/QC/compleasm_results/*/summary.txt"
files = sorted(glob.glob(results_path))

if not files:
    print(f"Error: No summary files found at {results_path}")
    print("Check your current working directory with 'os.getcwd()'")
else:
    print(f"| {'Dataset':<20} | {'S (%)':>8} | {'D (%)':>8} | {'F (%)':>8} | {'M (%)':>8} |")
    print(f"| {'-'*20} | {'-'*8}: | {'-'*8}: | {'-'*8}: | {'-'*8}: |")

    for summary_file in files:
        # Get folder name (e.g., p_ctg_eudicot)
        label = os.path.basename(os.path.dirname(summary_file))
        
        stats = {}
        with open(summary_file, 'r') as f:
            for line in f:
                # Look for patterns like "S:96.50%, 1955"
                match = re.search(r'([SDFM]):([\d.]+)%', line)
                if match:
                    key, val = match.groups()
                    stats[key] = val
        
        # Only print if we actually found data
        if stats:
            print(f"| {label:<20} | {stats.get('S', '0.0'):>7}% | {stats.get('D', '0.0'):>7}% | {stats.get('F', '0.0'):>7}% | {stats.get('M', '0.0'):>7}% |")

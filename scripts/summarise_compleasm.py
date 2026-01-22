import glob
import os
import re

# Path to your results
results_path = "03_hifiasm_assembly/QC/compleasm_results/*/summary.txt"
files = sorted(glob.glob(results_path))

if not files:
    print(f"Error: No summary files found at {results_path}")
else:
    # Header with expanded columns for absolute numbers and Total
    header = f"| {'Dataset':<22} | {'S (%, (N))':>15} | {'D (%, (N))':>15} | {'F (%, (N))':>15} | {'M (%, (N))':>15} | {'Total':>8} |"
    separator = f"| {'-'*22} | {'-'*15}: | {'-'*15}: | {'-'*15}: | {'-'*15}: | {'-'*8}: |"
    print(header)
    print(separator)

    for summary_file in files:
        label = os.path.basename(os.path.dirname(summary_file))
        
        stats = {}
        total_n = "0"
        
        with open(summary_file, 'r') as f:
            for line in f:
                # Capture standard lines: "S:96.50%, 1955"
                match = re.search(r'([SDFM]):([\d.]+)%,\s+(\d+)', line)
                if match:
                    key, perc, count = match.groups()
                    stats[key] = f"{perc}% ({count})"
                
                # Capture the Total line: "N:2026" or "n:2026"
                n_match = re.search(r'[Nn]:(\d+)', line)
                if n_match:
                    total_n = n_match.group(1)
        
        if stats:
            print(f"| {label:<22} | "
                  f"{stats.get('S', '0.0% (0)'):>15} | "
                  f"{stats.get('D', '0.0% (0)'):>15} | "
                  f"{stats.get('F', '0.0% (0)'):>15} | "
                  f"{stats.get('M', '0.0% (0)'):>15} | "
                  f"{total_n:>8} |")

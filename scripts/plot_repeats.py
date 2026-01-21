import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np

# Configuration
data_dir = '03_hifiasm_assembly/QC/telomere_results'
assemblies = {
    "p_ctg": "E_phylacis_asm.bp.p_ctg.fa",
    "hap1": "E_phylacis_asm.bp.hap1.p_ctg.fa",
    "hap2": "E_phylacis_asm.bp.hap2.p_ctg.fa"
}

motif_colors = {
    'AAACCCT': '#D55E00', 
    'ACCCGTC': '#56B4E9', 
    'AAAAAAG': '#009E73', 
    'AAGACTC': '#CC79A7'  
}
motifs = ['AAACCCT', 'ACCCGTC', 'AAAAAAG', 'AAGACTC']

for label, assm_file in assemblies.items():
    print(f"Generating high-contrast plot (v3) for {label}...")
    
    fig, axes = plt.subplots(11, 4, figsize=(16, 10), sharey='col')
    
    # Adjusting top margin to accommodate title and subtitle
    plt.subplots_adjust(hspace=0.05, wspace=0.15, left=0.05, right=0.98, top=0.88, bottom=0.05)
    
    # Main Title and Subtitle
    fig.suptitle(f"Repeat distribution for {assm_file}", fontsize=16, fontweight='bold', y=0.96)
    subtitle_text = "Method: Motifs identified via 'tidk explore' and mapped across 10kb windows using 'tidk search'"
    fig.text(0.5, 0.92, subtitle_text, ha='center', fontsize=11, style='italic')
    
    prefix = 'h1' if "hap1" in label else ('h2' if "hap2" in label else 'p')

    for row_idx in range(11):
        chrom = f"{prefix}tg0000{row_idx+1:02}l"
        
        for col_idx, motif in enumerate(motifs):
            ax = axes[row_idx, col_idx]
            color = motif_colors[motif]
            
            file_path = f"{data_dir}/{assm_file}_{motif}_telomeric_repeat_windows.tsv"
            
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, sep='\t')
                chrom_data = df[df['id'] == chrom].copy()
                
                if not chrom_data.empty:
                    chrom_data['total_count'] = chrom_data['forward_repeat_number'] + chrom_data['reverse_repeat_number']
                    
                    # 1. Smoothed trend line
                    smoothed = chrom_data['total_count'].rolling(window=5, center=True, min_periods=1).mean()
                    ax.plot(chrom_data['window'] / 1e6, smoothed, color=color, linewidth=1.6, zorder=1)
                    
                    # 2. Black points (size 0.4, alpha 0.5)
                    ax.scatter(chrom_data['window'] / 1e6, chrom_data['total_count'], 
                               color='black', s=0.4, alpha=0.5, edgecolors='none', zorder=2)

            # Style Cleanup
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.tick_params(axis='both', which='major', labelsize=6, pad=0.5, length=2)
            
            if col_idx == 0:
                ax.set_ylabel(f"C{row_idx+1}", fontsize=8, fontweight='bold', rotation=0, labelpad=15, verticalalignment='center')
            else:
                ax.set_yticklabels([])

            if row_idx == 0:
                ax.set_title(motif, fontsize=10, fontweight='bold', pad=2, color=color)
            
            ax.set_xlabel("Mb", fontsize=7, labelpad=0.5)

    output_path = f"{data_dir}/{label}_repeat_fingerprint.png"
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"File saved: {output_path}")

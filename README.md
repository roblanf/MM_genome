# MM_genome
An assembly of the Meelup Mallee genome.

# Get the environment running

I'll try and do most things in one conda environment. The details of what that looks like are in the `environment.yml` file.

```bash
conda env create environment.yml
```

# Raw Data

The raw data are located here:

```bash
raw_data="/home/raw_data/MM/2025_long_reads/ONT_gDNA192_1071_RL/E_Phylacis/20251117_1227_2C_PBE85256_e040940c/"
```

# QC

First let's examine the raw long reads carefully.

```bash
# 1. Run NanoPlot
# Added --maxrows to keep the HTML report manageable given your high coverage
NanoPlot -t 64 --fastq ${raw_data}/*.fastq.gz --maxrows 1000000 -o 01_NanoPlot_Output

# 2. Long-read K-mer Counting
mkdir -p kmc_long_reads
mkdir -p tmp_kmc  # Ensure a local tmp directory exists

# Find files using the variable
find ${raw_data} -name "*.fastq.gz" > files.txt

# Execute KMC
# Increased -m to 256 as requested to utilize your server's RAM
kmc -k21 -t64 -m256 -ci1 -cs10000 @files.txt kmc_long_reads/ tmp_kmc/

# 3. Generate Histogram for GenomeScope
kmc_tools transform kmc_long_reads/ dump -s long_read_histogram.txt
```

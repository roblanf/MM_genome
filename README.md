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
# 1. Run NanoPlot, 64 threads
NanoPlot -t 64 --fastq raw_data/*.fastq.gz -o 01_NanoPlot_Output

# 2. Long-read K-mer Counting (using -ci1 -cs10000 to catch low frequency k-mers)
mkdir kmc_long_reads
find raw_data/ -name "*.fastq.gz" > files.txt
kmc -k21 -t64 -m256 -ci1 -cs10000 @files.txt kmc_long_reads/ tmp/
kmc_tools transform kmc_long_reads/ dump -s long_read_histogram.txt

# 3. GenomeScope 2.0 (Run this via R or the web tool)
# genomescope2 -i long_read_histogram.txt -o 02_GenomeScope_LongRead -k 21
```

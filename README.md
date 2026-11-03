# Summary of the Genome Assembly Trials

This branch is purposed to: (i) summarise the findings from the other two branches, and (ii) finalise the parameters and steps conducted to assemble the *Eucalyptus x phylacis* genome.

## Table of Contents
- <a href="#01">Prerequisites</a>
- <a href="#02">QC</a>
    - <a href="#02a">Running Chopper</a>
    - <a href="#02b">Running Kraken2</a>
    - <a href="#02c">Running GenomeScope2</a>
- <a href="#03">Assembly using Hifiasm</a>
    - <a href="#03a">Using default parameters</a>
    - <a href="#03b">Adding the estimated genome size</a>
    - <a href="#03c">Estimating the peak homozygous coverage</a>
    - <a href="#03d">Rerunning with `--hom-cov 56`</a>
    - <a href="#03e">Trying lower `-s` values</a>

## <a id="01">Prerequisites</a>
This analysis requires a number of software and R packages to run. Here, I am going to install them using `conda`:
```
conda create -n meelupmallee
conda activate meelupmallee
conda install -c conda-forge -c bioconda chopper compleasm genomescope2 gfatools hifiasm kmc kraken2 merqury meryl mosdepth nanoplot seqkit
```

> [!NOTE]
> The raw ONT (Oxford Nanopore Technologies) long reads are stored as 70 compressed FASTQ files at `00_raw_data/`.

## <a id="02">QC</a>

### <a id="02a">Running Chopper</a>
First, I filtered the raw ONT reads with Chopper with a minimum read length of **15kb** and read quality of **10**.
```
# run NanoPlot pre-filtering
NanoPlot -t 128 --fastq 00_raw_data/*.fastq.gz -o 01_qc/01_nanoplot --title "E_phylacis_Raw_ONT"

# pipe pigz (decompression) -> Chopper (filtering) -> bgzip (compression)
pigz -dc -p 16 00_raw_data/*.fastq.gz | \
chopper -q 10 -l 15000 | \
bgzip -@ 128 > 01_qc/02_chopper/E_phylacis_filtered.fastq.gz

# run NanoPlot post-filtering
NanoPlot -t 128 --fastq 01_qc/02_chopper/E_phylacis_filtered.fastq.gz -o 01_qc/03_nanoplot --title "E_phylacis_Filtered_ONT"
```

### OUTPUT
|        Metric       |         Raw       |      Filtered     |
| :------------------ | -----------------:| -----------------:|
| Total bases         | 44,326,974,907.0  | 33,967,234,187.0  |
| Number of reads     | 3,064,194.0       | 1,345,923.0       |
| Mean read length    | 14,466.1          | 25,237.1          |
| Median read length  | 12,334.0          | 22,752.0          |
| Read length N50     | 22,177.0          | 25,598.0          |
| STDEV read length   | 11,705.0          | 9,033.2           |
| Mean read quality   | 17.2              | 17.8              |
| Median read quality | 19.6              | 20.2              |
| >Q10                | 3051048 (99.6%)   | 1345923 (100.0%)  |
| >Q15                | 2547395 (83.1%)   | 1165231 (86.6%)   |
| >Q20                | 1396296 (45.6%)   | 696853 (51.8%)    |
| >Q25                | 66131 (2.2%)      | 32580 (2.4%)      |
| >Q30                | 992 (0.0%)        | 10 (0.0%)         |

### <a id="02b">Running Kraken2</a>
I performed an additional analysis using Kraken2 to detect foreign DNA on the filtered ONT reads. First, I have to download the PlusPFP database:
```
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20261022.tar.gz
tar -xzf k2_pluspfp_*.tar.gz
```

I then ran Kraken2 with two confidence thresholds: **0** (default) and **0.2**. At default, a read is assigned to a taxonomic group on the basis of a single matching k-mer and thus might capture an increased proportion of false positives. Setting this threshold higher should better capture true positives from the dataset.
```
# run Kraken2 with default parameter
kraken2 --db 01_qc/04_kraken/pluspfp --threads 128 \
    --gzip-compressed \
    --report 01_qc/04_kraken/default.kreport \
    --report-minimizer-data \
    --use-names \
    --output 01_qc/04_kraken/default.kraken \
    01_qc/02_chopper/E_phylacis_filtered.fastq.gz

# run Kraken2 with --conf 0.2
kraken2 --db 01_qc/04_kraken/pluspfp --threads 128 \
    --confidence 0.2 \
    --gzip-compressed \
    --report 01_qc/04_kraken/conf02.kreport \
    --report-minimizer-data \
    --use-names \
    --output 01_qc/04_kraken/conf02.kraken \
    01_qc/02_chopper/E_phylacis_filtered.fastq.gz
```

### OUTPUT
| Category                        |  `--confidence 0`  |  `--confidence 0.2`  |
| :------------------------------ | -----------------: | -------------------: |
| Total reads                     | 1,345,923          | 1,345,923            |
| **Streptophyta (land plants)**  | 1,336,060 (99.27%) | 1,222,485 (90.83%)   |
| – *Eucalyptus* (*E. grandis*)   | 1,332,485 (99.00%) | 1,118,097 (83.07%)   |
| – Other plant genera            | 3,165 (0.24%)      | 96 (0.01%)           |
| – Unresolved at genus level     | 410 (0.03%)        | 104,292 (7.75%)      |
| **Human (*Homo sapiens*)**      | 3,661 (0.27%)      | 116 (0.01%)          |
| **Fungi**                       | 3,496 (0.26%)      | 0                    |
| **Unresolved at kingdom level** | 180 (0.01%)        | 10,842 (0.81%)       |
| **Bacteria**                    | 1,704 (0.13%)      | 18 (<0.01%)          |
| **Archaea**                     | 6 (<0.01%)         | 0                    |
| **Viruses**                     | 0                  | 0                    |
| **Other sequences**             | 2 (<0.01%)         | 0                    |
| **Unresolved at domain level**  | 681 (0.05%)        | 4,009 (0.30%)        |
| **Unclassified**                | 133 (0.01%)        | 108,453 (8.06%)      |

The proportion of foreign reads in the dataset is very low and unlikely to pass the depth criteria during contig formation, so I think we can go directly to the assembly step.

### <a id="02c">Running GenomeScope2</a>
I ran KMC and GenomeScope2 to estimate the heterozygosity and genome size of *E. x phylacis* based on the filtered ONT reads.
```
# store the path to the filtered reads
find 01_qc/02_chopper/ -name "*.fastq.gz" > 01_qc/05_genomescope/files.txt

# run KMC and create the histogram
kmc -k21 -t128 -ci1 -cs10000 @01_qc/05_genomescope/files.txt 01_qc/05_genomescope/kmc_db 01_qc/05_genomescope/
kmc_tools transform 01_qc/05_genomescope/kmc_db histogram 01_qc/05_genomescope/histogram.txt -ci1 -cs10000

# run GenomeScope2
genomescope2 -i 01_qc/05_genomescope/histogram.txt -o 01_qc/05_genomescope/ -k 21 -p 2
```

### OUTPUT
![Figure 01. GenomeScope2 transformed linear plot on filtered ONT reads.](./figures/fig01.png)

## <a id="03">Assembly using Hifiasm</a>

### <a id="03a">Using default parameters</a>
```
hifiasm \
    -o 02_assembly/default/E_phylacis_asm \
    -t 128 \
    --ont \
    -l 3 \
    --telo-m AAACCCT \
    --dual-scaf \
    01_qc/02_chopper/E_phylacis_filtered.fastq.gz \
    2>&1 | tee 02_assembly/default/hifiasm_restart.log
```

After the third correction round, Hifiasm automatically sets the homozygous read coverage threshold to be **52**, even though the peak heterozygous coverage is **29**:

![Figure 02. Corrected distribution of read coverage from Hifiasm with default parameters.](./figures/fig02.png)

### <a id="03b">Adding the estimated genome size</a>
```
hifiasm \
    -o 02_assembly/hgsize533/E_phylacis_asm \
    -t 128 \
    --ont \
    -l 3 \
    --telo-m AAACCCT \
    --hg-size 533m \
    --dual-scaf \
    01_qc/02_chopper/E_phylacis_filtered.fastq.gz \
    2>&1 | tee 02_assembly/hgsize533/hifiasm_restart.log
```

The homozygous read coverage threshold is now **63** (i.e., total read length / estimated genome size), while the peak heterozygous coverage is still **29**. It seems that we need a more objective way to set the threshold...

### <a id="03c">Estimating the peak homozygous coverage</a>
To empirically set the homozygous peak coverage, I calculated the median BUSCO coverage across the two haplotypes. First, I concatenated the two haplotypes inferred from the default run and remapped the filtered ONT reads.
```
# change the working directory
cd 02_assembly/default/

# convert GFA to FASTA format
gfatools gfa2fa E_phylacis_asm.bp.hap1.p_ctg.gfa > hap1.fa
gfatools gfa2fa E_phylacis_asm.bp.hap2.p_ctg.gfa > hap2.fa

# concatenate both haplotypes
cat hap1.fa hap2.fa > both_haps.fa

# map the filtered ONT reads to both haplotypes
minimap2 \
    -t 128 \
    -ax lr:hq both_haps.fa \
    ../../01_qc/02_chopper/E_phylacis_filtered.fastq.gz | \
    samtools sort -@16 -o both_haps.bam
samtools index both_haps.bam
```

I then ran compleasm on each haplotype, extracted single-copy BUSCO, and calculated the depth statistics across loci.
```
# download the ODB database
compleasm download eudicotyledons --odb odb12 -L busco_lineages/

# run compleasm on each haplotype
compleasm run -a hap1.fa -o compleasm/compleasm_hap1 -l eudicotyledons --odb odb12 -t 128 -L busco_lineages/
compleasm run -a hap2.fa -o compleasm/compleasm_hap2 -l eudicotyledons --odb odb12 -t 128 -L busco_lineages/

# extract all single-copy BUSCO
for h in hap1 hap2; do
  awk -F'\t' 'NR>1 && $2=="Single"{print $3"\t"$4-1"\t"$5"\t"$1}' compleasm/compleasm_$h/eudicotyledons_odb12/full_table.tsv
done | sort -k1,1 -k2,2n > compleasm/busco_single.bed

# calculate the depth across BUSCO
mosdepth -t 16 --by compleasm/busco_single.bed --no-per-base compleasm/busco both_haps.bam

# extract the Q1, Q2, and Q3 depth across BUSCO
zcat compleasm/busco.regions.bed.gz | \
    cut -f5 | sort -n | \
    awk '{a[NR]=$1} END{printf "n=%d  Q1=%.1f  median=%.1f  Q3=%.1f\n", NR, a[int(NR*0.25)], a[int(NR*0.5)], a[int(NR*0.75)]}'

# change the working directory
cd ../..
```

This command prints the following: `n=5173  Q1=23.6  median=27.9  Q3=32.2`. The median closely matches the heterozygosity peak from <a href="#03a">Hifiasm</a> and <a href="#03b">GenomeScope2</a>. Thus, we should redo the Hifiasm run with **`--hom-cov 55`**.

### <a id="03d">Rerunning with `--hom-cov 55`</a>
```
hifiasm \
    -o 02_assembly/homcov55/E_phylacis_asm \
    -t 128 \
    --ont \
    -l 3 \
    --telo-m AAACCCT \
    --hom-cov 55 \
    --dual-scaf \
    01_qc/02_chopper/E_phylacis_filtered.fastq.gz \
    2>&1 | tee 02_assembly/homcov55/hifiasm_restart.log
```

The next question is whether we need to adjust the `-s` flag, which reflects the *"similarity threshold for duplicate haplotigs that should be purged"*. At `-l 3`, the default threshold is 0.55, meaning that a pair of contigs need to have at least 55% cross-parent overlap to be considered as homologous. This might be too strict for an F1 hybrid genome.

### <a id="03e">Trying lower `-s` values</a>
Let's try five different `-s` values (0.50, 0.45, 0.40, 0.35, 0.30) based on the corrected reads from the Hifiasm run. Here's the example code for one Hifiasm run:
```
# change the working directory
cd 02_assembly/homcov55

# set path to the original run
ORIG=/home/jeremias/meelupmallee/long_reads/02_assembly/homcov55/E_phylacis_asm

# set output path
NEW=E_phylacis_s050.asm

# create a new directory and shortcuts
mkdir -p hifiasm_s050 && cd hifiasm_s050
for f in ec.bin ovlp.source.bin ovlp.reverse.bin; do
    ln -s ${ORIG}.${f} ${NEW}.${f}
done

# re-run Hifiasm under the new parameter
hifiasm --ont -t 128 -o ${NEW} --hom-cov 55 -l 3 -s 0.50 --dual-scaf --telo-m AAACCCT /dev/null 2> ${NEW}.log

# change the working directory
cd ../../..
```

I converted all GFA outputs to FASTA format, then calculated the sequence statistics:
```
# iterate over all GFA files
for g in 02_assembly/homcov55/hifiasm_s0*/*.p_ctg.gfa; do
    gfatools gfa2fa "$g" > "${g%.gfa}.fa"
done

# calculate the statistics of all FASTA files
seqkit stats 02_assembly/homcov55/hifiasm_s0*/*.p_ctg.fa > hifiasm_s0XX_seqkit.txt
```

### OUTPUT

### Haplotype 1
| `-s` | Contigs | Total length (bp) | Min (bp) | Mean (bp) | Max (bp)   |
| :--: | ------: | ----------------: | -------: | --------: | ---------: |
| 0.30 |     449 |       583,208,380 |   30,917 | 1,298,905 | 65,682,487 |
| 0.35 |     405 |       579,832,984 |   24,809 | 1,431,686 | 65,682,487 |
| 0.40 |     421 |       573,216,200 |   24,809 | 1,361,559 | 65,682,487 |
| 0.45 |     431 |       576,202,291 |   24,809 | 1,336,896 | 65,682,487 |
| 0.50 |     402 |       574,198,743 |   24,809 | 1,428,355 | 65,682,487 |
| 0.55 |     440 |       575,869,290 |   24,809 | 1,308,794 | 65,682,487 |

### Haplotype 2
| `-s` | Contigs | Total length (bp) | Min (bp) | Mean (bp) | Max (bp)   |
| :--: | ------: | ----------------: | -------: | --------: | ---------: |
| 0.30 |     308 |       563,732,156 |   24,809 | 1,830,299 | 61,702,390 |
| 0.35 |     348 |       565,571,274 |   36,811 | 1,625,205 | 61,702,390 |
| 0.40 |     319 |       571,160,799 |   36,811 | 1,790,473 | 61,702,390 |
| 0.45 |     320 |       566,439,604 |   36,811 | 1,770,124 | 61,702,390 |
| 0.50 |     305 |       560,735,642 |   36,811 | 1,838,478 | 61,702,390 |
| 0.55 |     293 |       559,542,586 |   36,811 | 1,909,702 | 61,702,390 |

### Primary
| `-s` | Contigs | Total length (bp) | Min (bp) | Mean (bp)  | Max (bp)   |
| :--: | ------: | ----------------: | -------: | ---------: | ---------: |
| 0.30 |      55 |       538,510,568 |   30,917 |  9,791,101 | 61,704,473 |
| 0.35 |      54 |       538,043,526 |   30,917 |  9,963,769 | 61,704,473 |
| 0.40 |      53 |       538,126,685 |   30,917 | 10,153,334 | 61,704,473 |
| 0.45 |      52 |       538,070,186 |   30,917 | 10,347,504 | 61,704,473 |
| 0.50 |      56 |       538,678,952 |   30,917 |  9,619,267 | 61,704,473 |
| 0.55 |      58 |       545,448,847 |   30,917 |  9,404,291 | 61,704,473 |

The results seem to be pretty consistent across runs, so let's check each haplotype's BUSCO completeness and k-mer completeness. Here's the example code for one set of haplotypes:
```
# build the read k-mer database once
meryl count k=21 threads=128 01_qc/02_chopper/E_phylacis_filtered.fastq.gz output 02_assembly/homcov55/merqury/reads.meryl

# change the working directory
cd 02_assembly/homcov55/

# run compleasm
compleasm run -a hifiasm_s030/E_phylacis_s030.asm.bp.hap1.p_ctg.fa -o compleasm/s030_hap1 -l eudicotyledons --odb odb12 -t 128 -L ../default/busco_lineages/
compleasm run -a hifiasm_s030/E_phylacis_s030.asm.bp.hap2.p_ctg.fa -o compleasm/s030_hap2 -l eudicotyledons --odb odb12 -t 128 -L ../default/busco_lineages/
compleasm run -a hifiasm_s030/E_phylacis_s030.union.fa -o compleasm/s030_union -l eudicotyledons --odb odb12 -t 128 -L ../default/busco_lineages/

# run merqury for haplotype 1
mkdir -p merqury/s030_hap1 && cd merqury/s030_hap1
merqury.sh ../reads.meryl ../../hifiasm_s030/E_phylacis_s030.asm.bp.hap1.p_ctg.fa s030_hap1_merq
cd ../..

# run merqury for haplotype 2
mkdir -p merqury/s030_hap2 && cd merqury/s030_hap2
merqury.sh ../reads.meryl ../../hifiasm_s030/E_phylacis_s030.asm.bp.hap2.p_ctg.fa s030_hap2_merq
cd ../..

# run merqury for the combined haplotypes
mkdir -p merqury/s030_union && cd merqury/s030_union
merqury.sh ../reads.meryl ../../hifiasm_s030/E_phylacis_s030.union.fa s030_union_merq
cd ../..

# change the working directory
cd ../..
```

I asked Claude to generate a python <a href="./scripts/summarise_s_sweep.py">script</a> to summarise the results. The following tables show the statistics across 2,805 BUSCOs (for Compleasm) and 601,145,508 solid read k-mers (for Merqury) per run:

### Haplotype 1
| Metric                 |        0.30 |        0.35 |        0.40 |        0.45 |        0.50 |        0.55 |
| ---------------------- | ----------: | ----------: | ----------: | ----------: | ----------: | ----------: |
| BUSCO S (%)            |       92.16 |       93.05 |       92.55 |       92.73 |       93.44 |       92.69 |
| BUSCO D (%)            |        7.42 |        6.56 |        7.09 |        6.88 |        6.20 |        6.95 |
| BUSCO F (%)            |        0.39 |        0.36 |        0.32 |        0.36 |        0.32 |        0.32 |
| BUSCO M (%)            |        0.04 |        0.04 |        0.04 |        0.04 |        0.04 |        0.04 |
| k-mer completeness (%) |     66.6775 |     66.2871 |     65.9580 |     66.1271 |     65.9356 |     66.1318 |
| QV                     |     64.3232 |     64.5597 |     64.2477 |     64.4424 |     64.4342 |     64.4088 |

### Haplotype 2
| Metric                 |        0.30 |        0.35 |        0.40 |        0.45 |        0.50 |        0.55 |
| ---------------------- | ----------: | ----------: | ----------: | ----------: | ----------: | ----------: |
| BUSCO S (%)            |       94.01 |       93.26 |       93.94 |       94.08 |       94.40 |       93.90 |
| BUSCO D (%)            |        5.70 |        6.45 |        5.74 |        5.63 |        5.28 |        5.17 |
| BUSCO F (%)            |        0.25 |        0.25 |        0.29 |        0.25 |        0.29 |        0.32 |
| BUSCO M (%)            |        0.04 |        0.04 |        0.04 |        0.04 |        0.04 |        0.61 |
| k-mer completeness (%) |     65.3983 |     65.6864 |     65.9052 |     65.6105 |     65.2542 |     64.9363 |
| QV                     |     65.5825 |     65.3600 |     65.7464 |     65.4516 |     65.6166 |     65.9861 |

### Union (hap1 + hap2)
| Metric                 |        0.30 |        0.35 |        0.40 |        0.45 |        0.50 |        0.55 |
| ---------------------- | ----------: | ----------: | ----------: | ----------: | ----------: | ----------: |
| BUSCO S (%)            |        2.96 |        2.96 |        2.89 |        2.92 |        2.92 |        2.85 |
| BUSCO D (%)            |       96.83 |       96.83 |       96.90 |       96.86 |       96.86 |       96.93 |
| BUSCO F (%)            |        0.18 |        0.18 |        0.18 |        0.18 |        0.18 |        0.18 |
| BUSCO M (%)            |        0.04 |        0.04 |        0.04 |        0.04 |        0.04 |        0.04 |
| k-mer completeness (%) |     99.3074 |     99.3073 |     99.3073 |     99.3073 |     99.3071 |     99.3071 |
| QV                     |     64.8963 |     64.9374 |     64.9329 |     64.9145 |     64.9792 |     65.1172 |

We also calculated the copy number of heterozygous peak (14–44x) and homozygous peak (44–69x) read k-mers in the assembly (from Merqury `spectra-cn.hist`):
- For one haplotype, `het k-mers: 0x` and `het k-mers: 1x` should be ~50% each, while `hom k-mers: 1x` should be ~100%
- For the union, `het k-mers: 1x` and `hom k-mers: 2x` should be ~100% each
- `het k-mers: ≥2x` and `hom k-mers: ≥3x` mark regions present as both parental copies in both haplotypes

### Haplotype 1
| Metric           | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 | 0.55 |
| ---------------- | ---: | ---: | ---: | ---: | ---: | ---: |
| het k-mers: 0x   | 44.4 | 44.8 | 45.2 | 45.0 | 45.3 | 45.0 |
| het k-mers: 1x   | 55.2 | 54.7 | 54.4 | 54.6 | 54.3 | 54.6 |
| het k-mers: ≥2x  |  0.4 |  0.4 |  0.4 |  0.4 |  0.4 |  0.4 |
| hom k-mers: 0x   |  8.5 |  8.6 |  8.9 |  8.8 |  8.7 |  8.9 |
| hom k-mers: 1x   | 77.6 | 78.2 | 77.6 | 77.7 | 78.5 | 77.5 |
| hom k-mers: 2x   | 13.7 | 13.0 | 13.2 | 13.2 | 12.6 | 13.3 |
| hom k-mers: ≥3x  |  0.3 |  0.3 |  0.2 |  0.3 |  0.3 |  0.3 |

### Haplotype 2
| Metric           | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 | 0.55 |
| ---------------- | ---: | ---: | ---: | ---: | ---: | ---: |
| het k-mers: 0x   | 45.8 | 45.4 | 45.2 | 45.6 | 46.0 | 46.3 |
| het k-mers: 1x   | 53.8 | 54.2 | 54.4 | 54.0 | 53.7 | 53.3 |
| het k-mers: ≥2x  |  0.4 |  0.4 |  0.4 |  0.4 |  0.4 |  0.4 |
| hom k-mers: 0x   |  9.3 |  9.2 |  8.9 |  9.1 |  9.4 |  9.6 |
| hom k-mers: 1x   | 78.3 | 77.9 | 78.5 | 78.5 | 78.7 | 78.8 |
| hom k-mers: 2x   | 12.2 | 12.7 | 12.3 | 12.1 | 11.6 | 11.3 |
| hom k-mers: ≥3x  |  0.3 |  0.3 |  0.3 |  0.3 |  0.3 |  0.2 |

### Union (hap1 + hap2)
| Metric (%)       | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 | 0.55 |
| ---------------- | ---: | ---: | ---: | ---: | ---: | ---: |
| het k-mers: 0x   |  0.6 |  0.6 |  0.6 |  0.6 |  0.6 |  0.6 |
| het k-mers: 1x   | 88.5 | 88.6 | 88.8 | 88.9 | 89.6 | 89.6 |
| het k-mers: ≥2x  | 10.9 | 10.8 | 10.7 | 10.5 |  9.8 |  9.8 |
| hom k-mers: 0x   |  0.2 |  0.2 |  0.2 |  0.2 |  0.2 |  0.2 |
| hom k-mers: 1x   |  4.1 |  4.1 |  4.2 |  4.2 |  4.2 |  4.3 |
| hom k-mers: 2x   | 83.0 | 83.2 | 83.3 | 83.5 | 84.7 | 84.6 |
| hom k-mers: ≥3x  | 12.7 | 12.5 | 12.3 | 12.1 | 10.9 | 11.0 |

Based on these results, roughly 10% of the genome was duplicated within the first haplotype, and another 10% on the other haplotype. `-s 0.50` seems to be the most optimal parameter, with the highest number of single-copy BUSCO recovered for both haplotypes and comparable k-mer completeness. 

---
*Last update: 03 November 2026 by Jeremias Ivan*
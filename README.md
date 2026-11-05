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
    - <a href="#03b">Estimating the peak homozygous coverage</a>
    - <a href="#03c">Rerunning with `--hom-cov 55`</a>
    - <a href="#03d">Trying other `-s` values</a>
    - <a href="#03e">Conclusion</a>
- <a href="#04">Scaffolding</a>

## <a id="01">Prerequisites</a>
This analysis requires a number of software and R packages to run. Here, I am going to install them using `conda`:
```
conda create -n meelupmallee
conda activate meelupmallee
conda install -c conda-forge -c bioconda \
    blobtoolkit \
    chopper \
    compleasm \
    quast \
    genomescope2 \
    gfatools \
    hifiasm \
    kmc \
    kraken2 \
    merqury \
    meryl \
    mosdepth \
    nanoplot \
    r-optparse \
    r-plotly \
    r-tidyverse \
    seqkit
```

> [!NOTE]
> The raw ONT (Oxford Nanopore Technologies) long reads are stored as 70 compressed FASTQ files at `00_raw_data/`.

## <a id="02">QC</a>

### <a id="02a">Running Chopper</a>
First, I filtered the raw ONT reads with Chopper with a minimum read length of **15kb** and read quality of **10**.
```
# run NanoPlot pre-filtering
NanoPlot \
    -t 128 \
    --fastq 00_raw_data/*.fastq.gz \
    -o 01_qc/01_nanoplot \
    --title "E_phylacis_Raw_ONT"

# pipe pigz (decompression) -> Chopper (filtering) -> bgzip (compression)
pigz -dc -p 16 00_raw_data/*.fastq.gz | \
chopper -q 10 -l 15000 | \
bgzip -@ 128 > 01_qc/02_chopper/E_phylacis_filtered.fastq.gz

# run NanoPlot post-filtering
NanoPlot \
    -t 128 \
    --fastq 01_qc/02_chopper/E_phylacis_filtered.fastq.gz \
    -o 01_qc/03_nanoplot \
    --title "E_phylacis_Filtered_ONT"
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
kraken2 \
    --db 01_qc/04_kraken/pluspfp \
    --threads 128 \
    --gzip-compressed \
    --report 01_qc/04_kraken/default.kreport \
    --report-minimizer-data \
    --use-names \
    --output 01_qc/04_kraken/default.kraken \
    01_qc/02_chopper/E_phylacis_filtered.fastq.gz

# run Kraken2 with --conf 0.2
kraken2 \
    --db 01_qc/04_kraken/pluspfp \
    --threads 128 \
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
kmc \
    -k21 \
    -t128 \
    -ci1 -cs10000 \
    @01_qc/05_genomescope/files.txt \
    01_qc/05_genomescope/kmc_db \
    01_qc/05_genomescope/

kmc_tools transform \
    01_qc/05_genomescope/kmc_db \
    histogram \
    01_qc/05_genomescope/histogram.txt \
    -ci1 -cs10000

# run GenomeScope2
genomescope2 \
    -i 01_qc/05_genomescope/histogram.txt \
    -o 01_qc/05_genomescope/ \
    -k 21 -p 2
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

After the third correction round, Hifiasm automatically sets the homozygous read coverage threshold to be **52**, even though the peak homozygous coverage ranges between 54 to 56 based on the histogram:

![Figure 02. Corrected distribution of read coverage from Hifiasm with default parameters.](./figures/fig02.png)

### <a id="03b">Estimating the peak homozygous coverage</a>
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

I then ran compleasm on each haplotype, extracted all single-copy BUSCO, and calculated the depth statistics across loci.
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

### <a id="03c">Rerunning with `--hom-cov 55`</a>
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

### <a id="03d">Trying other `-s` values</a>
Let's try different `-s` values (0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.65, 0.70, 0.75) based on the corrected reads from the Hifiasm run. I also tried a stricter `-l 2` with four different `-s` values (0.45, 0.55, 0.65, 0.75). Here's the example code for one Hifiasm run:
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

I then converted all GFA outputs to FASTA format.
```
for g in 02_assembly/homcov55/hifiasm_s0*/*.p_ctg.gfa; do
    gfatools gfa2fa "$g" > "${g%.gfa}.fa"
done
```

### SeqKit
```
seqkit stats 02_assembly/homcov55/hifiasm_s0*/*.p_ctg.fa > hifiasm_all_seqkit.txt
```

- <b><i>Haplotype 1</b></i>

    | `-l` | `-s` | Contigs | Total length (bp) | Min (bp) | Mean (bp) |   Max (bp) |
    | :--: | :--: | ------: | ----------------: | -------: | --------: | ---------: |
    |  3   | 0.30 |     449 |       583,208,380 |   30,917 | 1,298,905 | 65,682,487 |
    |  3   | 0.35 |     405 |       579,832,984 |   24,809 | 1,431,686 | 65,682,487 |
    |  3   | 0.40 |     421 |       573,216,200 |   24,809 | 1,361,559 | 65,682,487 |
    |  3   | 0.45 |     431 |       576,202,291 |   24,809 | 1,336,896 | 65,682,487 |
    |  3   | 0.50 |     402 |       574,198,743 |   24,809 | 1,428,355 | 65,682,487 |
    |  3   | 0.55 |     440 |       575,869,290 |   24,809 | 1,308,794 | 65,682,487 |
    |  3   | 0.60 |     402 |       570,069,076 |   24,809 | 1,418,082 | 65,682,487 |
    |  3   | 0.65 |     414 |       546,652,395 |   24,809 | 1,320,416 | 65,682,487 |
    |  3   | 0.70 |     413 |       601,523,016 |   30,917 | 1,456,472 | 61,713,612 |
    |  3   | 0.75 |     329 |       697,630,498 |   30,917 | 2,120,457 | 65,682,489 |
    |  2   | 0.45 |     423 |       581,017,343 |   24,809 | 1,373,564 | 65,682,487 |
    |  2   | 0.55 |     406 |       570,533,238 |   30,917 | 1,405,254 | 65,682,487 |
    |  2   | 0.65 |     418 |       596,630,555 |   30,917 | 1,427,346 | 65,682,487 |
    |  2   | 0.75 |     410 |       620,392,369 |   30,917 | 1,513,152 | 65,682,489 |

- <b><i>Haplotype 2</b></i>

    | `-l` | `-s` | Contigs | Total length (bp) | Min (bp) | Mean (bp) |   Max (bp) |
    | :--: | :--: | ------: | ----------------: | -------: | --------: | ---------: |
    |  3   | 0.30 |     308 |       563,732,156 |   24,809 | 1,830,299 | 61,702,390 |
    |  3   | 0.35 |     348 |       565,571,274 |   36,811 | 1,625,205 | 61,702,390 |
    |  3   | 0.40 |     319 |       571,160,799 |   36,811 | 1,790,473 | 61,702,390 |
    |  3   | 0.45 |     320 |       566,439,604 |   36,811 | 1,770,124 | 61,702,390 |
    |  3   | 0.50 |     305 |       560,735,642 |   36,811 | 1,838,478 | 61,702,390 |
    |  3   | 0.55 |     293 |       559,542,586 |   36,811 | 1,909,702 | 61,702,390 |
    |  3   | 0.60 |     351 |       567,535,942 |   43,064 | 1,616,912 | 61,702,390 |
    |  3   | 0.65 |     329 |       589,961,923 |   41,601 | 1,793,197 | 61,702,390 |
    |  3   | 0.70 |     285 |       638,319,662 |   24,809 | 2,239,718 | 65,543,494 |
    |  3   | 0.75 |     325 |       587,982,618 |   35,019 | 1,809,177 | 65,543,479 |
    |  2   | 0.45 |     328 |       560,484,069 |   36,811 | 1,708,793 | 61,702,390 |
    |  2   | 0.55 |     321 |       564,637,389 |   42,060 | 1,758,995 | 61,702,390 |
    |  2   | 0.65 |     330 |       539,894,173 |   36,811 | 1,636,043 | 61,702,390 |
    |  2   | 0.75 |     266 |       669,882,354 |   35,019 | 2,518,355 | 65,543,479 |

- <b><i>Primary</b></i>

    | `-l` | `-s` | Contigs | Total length (bp) | Min (bp) |  Mean (bp) |   Max (bp) |
    | :--: | :--: | ------: | ----------------: | -------: | ---------: | ---------: |
    |  3   | 0.30 |      55 |       538,510,568 |   30,917 |  9,791,101 | 61,704,473 |
    |  3   | 0.35 |      54 |       538,043,526 |   30,917 |  9,963,769 | 61,704,473 |
    |  3   | 0.40 |      53 |       538,126,685 |   30,917 | 10,153,334 | 61,704,473 |
    |  3   | 0.45 |      52 |       538,070,186 |   30,917 | 10,347,504 | 61,704,473 |
    |  3   | 0.50 |      56 |       538,678,952 |   30,917 |  9,619,267 | 61,704,473 |
    |  3   | 0.55 |      58 |       545,448,847 |   30,917 |  9,404,290 | 61,704,473 |
    |  3   | 0.60 |      57 |       608,733,002 |   30,917 | 10,679,526 | 61,704,473 |
    |  3   | 0.65 |      58 |       633,708,634 |   30,917 | 10,926,011 | 61,704,473 |
    |  3   | 0.70 |      63 |       844,781,570 |   30,917 | 13,409,231 | 61,713,612 |
    |  3   | 0.75 |      64 |       960,983,357 |   30,917 | 15,015,365 | 65,682,489 |
    |  2   | 0.45 |      56 |       544,895,664 |   30,917 |  9,730,280 | 61,704,473 |
    |  2   | 0.55 |      56 |       545,452,896 |   30,917 |  9,740,230 | 61,704,473 |
    |  2   | 0.65 |      58 |       635,269,094 |   30,917 | 10,952,915 | 61,704,473 |
    |  2   | 0.75 |      64 |       961,221,210 |   30,917 | 15,019,081 | 65,682,489 |

### Compleasm + Merqury
First, we need to create a k-mer database from the filtered ONT reads.
```
meryl count k=21 threads=128 01_qc/02_chopper/E_phylacis_filtered.fastq.gz output 02_assembly/homcov55/merqury/reads.meryl
```

Now, let's check each haplotype's BUSCO completeness and k-mer completeness. Here's the example code for one set of haplotypes:
```
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

- <b><i>Haplotype 1</b></i>

    | Metric                      |  0.30 |  0.35 |  0.40 |  0.45 |  0.50 |  0.55 |  0.60 |  0.65 |  0.70 |  0.75 | 0.45 (l2) | 0.55 (l2) | 0.65 (l2) | 0.75 (l2) |
    | --------------------------- | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | --------: | --------: | --------: | --------: |
    | BUSCO S (%)                 | 92.16 | 93.05 | 92.55 | 92.73 | 93.44 | 92.69 | 93.58 | 91.30 | 63.28 | 48.70 |     92.30 |     93.19 |     85.78 |     50.41 |
    | BUSCO D (%)                 |  7.42 |  6.56 |  7.09 |  6.88 |  6.20 |  6.95 |  5.49 |  5.42 | 22.07 | 37.58 |      7.31 |      5.88 |     10.98 |     31.41 |
    | BUSCO F (%)                 |  0.39 |  0.36 |  0.32 |  0.36 |  0.32 |  0.32 |  0.39 |  0.57 |  1.25 |  0.86 |      0.36 |      0.39 |      0.53 |      1.18 |
    | BUSCO&nbsp;M&nbsp;(%)       |  0.04 |  0.04 |  0.04 |  0.04 |  0.04 |  0.04 |  0.53 |  2.71 | 13.40 | 12.87 |      0.04 |      0.53 |      2.71 |     17.01 |
    | k-mer completeness&nbsp;(%) | 66.68 | 66.29 | 65.96 | 66.13 | 65.94 | 66.13 | 65.55 | 63.06 | 63.65 | 71.98 |     66.46 |     65.46 |     66.57 |     64.97 |
    | Solid k-mers found (M)      | 400.8 | 398.5 | 396.5 | 397.5 | 396.4 | 397.5 | 394.1 | 379.1 | 382.6 | 432.7 |     399.5 |     393.5 |     400.2 |     390.6 |
    | QV                          | 64.32 | 64.56 | 64.25 | 64.44 | 64.43 | 64.41 | 64.55 | 64.47 | 64.35 | 64.93 |     64.60 |     64.58 |     64.40 |     64.62 |
    | Asm-only k-mers             | 4,521 | 4,243 | 4,507 | 4,332 | 4,326 | 4,357 | 4,197 | 4,093 | 4,614 | 4,702 |     4,203 |     4,156 |     4,546 |     4,478 |

- <b><i>Haplotype 2</b></i>

    | Metric                      |  0.30 |  0.35 |  0.40 |  0.45 |  0.50 |  0.55 |  0.60 |  0.65 |  0.70 |  0.75 | 0.45 (l2) | 0.55 (l2) | 0.65 (l2) | 0.75 (l2) |
    | --------------------------- | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | --------: | --------: | --------: | --------: |
    | BUSCO S (%)                 | 94.01 | 93.26 | 93.94 | 94.08 | 94.40 | 93.90 | 93.19 | 91.19 | 63.53 | 56.08 |     93.69 |     93.58 |     87.31 |     56.08 |
    | BUSCO D (%)                 |  5.70 |  6.45 |  5.74 |  5.63 |  5.28 |  5.17 |  6.52 |  8.52 | 27.59 | 28.98 |      5.42 |      6.13 |      7.81 |     34.33 |
    | BUSCO F (%)                 |  0.25 |  0.25 |  0.29 |  0.25 |  0.29 |  0.32 |  0.25 |  0.25 |  0.82 |  1.28 |      0.29 |      0.25 |      0.71 |      1.00 |
    | BUSCO&nbsp;M&nbsp;(%)       |  0.04 |  0.04 |  0.04 |  0.04 |  0.04 |  0.61 |  0.04 |  0.04 |  8.06 | 13.65 |      0.61 |      0.04 |      4.17 |      8.59 |
    | k-mer completeness&nbsp;(%) | 65.40 | 65.69 | 65.91 | 65.61 | 65.25 | 64.94 | 65.64 | 67.22 | 67.91 | 63.26 |     65.05 |     65.57 |     62.04 |     70.06 |
    | Solid k-mers found (M)      | 393.1 | 394.9 | 396.2 | 394.4 | 392.3 | 390.4 | 394.6 | 404.1 | 408.3 | 380.3 |     391.0 |     394.2 |     373.0 |     421.2 |
    | QV                          | 65.58 | 65.36 | 65.75 | 65.45 | 65.62 | 65.99 | 65.43 | 65.61 | 65.56 | 65.15 |     65.32 |     65.53 |     65.67 |     65.33 |
    | Asm-only k-mers             | 3,265 | 3,457 | 3,193 | 3,389 | 3,227 | 2,961 | 3,401 | 3,387 | 3,695 | 3,775 |     3,455 |     3,321 |     3,058 |     4,122 |

- <b><i>Union (hap1 + hap2)</b></i>

    | Metric                      |  0.30 |  0.35 |  0.40 |  0.45 |  0.50 |  0.55 |  0.60 |  0.65 |  0.70 |  0.75 | 0.45 (l2) | 0.55 (l2) | 0.65 (l2) | 0.75 (l2) |
    | --------------------------- | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | --------: | --------: | --------: | --------: |
    | BUSCO S (%)                 |  2.96 |  2.96 |  2.89 |  2.92 |  2.92 |  2.85 |  2.92 |  2.92 |  2.78 |  2.57 |      2.89 |      2.96 |      2.89 |      2.60 |
    | BUSCO D (%)                 | 96.83 | 96.83 | 96.90 | 96.86 | 96.86 | 96.93 | 96.86 | 96.86 | 97.01 | 97.22 |     96.90 |     96.83 |     96.90 |     97.18 |
    | BUSCO F (%)                 |  0.18 |  0.18 |  0.18 |  0.18 |  0.18 |  0.18 |  0.18 |  0.18 |  0.18 |  0.18 |      0.18 |      0.18 |      0.18 |      0.18 |
    | BUSCO&nbsp;M&nbsp;(%)       |  0.04 |  0.04 |  0.04 |  0.04 |  0.04 |  0.04 |  0.04 |  0.04 |  0.04 |  0.04 |      0.04 |      0.04 |      0.04 |      0.04 |
    | k-mer completeness&nbsp;(%) | 99.31 | 99.31 | 99.31 | 99.31 | 99.31 | 99.31 | 99.31 | 99.31 | 99.31 | 99.31 |     99.31 |     99.31 |     99.31 |     99.31 |
    | Solid k-mers found (M)      | 597.0 | 597.0 | 597.0 | 597.0 | 597.0 | 597.0 | 597.0 | 597.0 | 597.0 | 597.0 |     597.0 |     597.0 |     597.0 |     597.0 |
    | QV                          | 64.90 | 64.94 | 64.93 | 64.91 | 64.98 | 65.12 | 64.96 | 65.02 | 64.93 | 65.03 |     64.94 |     65.02 |     64.95 |     64.98 |
    | Asm-only k-mers             | 7,786 | 7,700 | 7,700 | 7,721 | 7,553 | 7,318 | 7,598 | 7,480 | 8,309 | 8,477 |     7,658 |     7,477 |     7,604 |     8,600 |

We also calculated the copy number of heterozygous peak (14–44x) and homozygous peak (44–69x) read k-mers in the assembly (from Merqury `spectra-cn.hist`).

> [!NOTE]
> - For one haplotype, `het k-mers: 0x` and `het k-mers: 1x` should be ~50% each, while `hom k-mers: 1x` should be ~100%
> - For the union, `het k-mers: 1x` and `hom k-mers: 2x` should be ~100% each
> - `het k-mers: ≥2x` and `hom k-mers: ≥3x` mark regions present as both parental copies in both haplotypes

- <b><i>Haplotype 1</b></i>

    | Metric               | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.45 (l2) | 0.55 (l2) | 0.65 (l2) | 0.75 (l2) |
    | -------------------- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --------: | --------: | --------: | --------: |
    | het k-mers: 0x       | 44.4 | 44.8 | 45.2 | 45.0 | 45.3 | 45.0 | 45.6 | 47.9 | 43.9 | 34.2 |      44.5 |      45.7 |      43.4 |      41.5 |
    | het k-mers: 1x       | 55.2 | 54.7 | 54.4 | 54.6 | 54.3 | 54.6 | 54.0 | 51.8 | 55.1 | 63.9 |      55.0 |      53.9 |      55.9 |      57.0 |
    | het k-mers:&nbsp;≥2x |  0.4 |  0.4 |  0.4 |  0.4 |  0.4 |  0.4 |  0.4 |  0.4 |  1.0 |  1.8 |       0.4 |       0.4 |       0.6 |       1.5 |
    | hom k-mers: 0x       |  8.5 |  8.6 |  8.9 |  8.8 |  8.7 |  8.9 |  9.4 | 12.5 | 21.0 | 15.3 |       8.8 |       9.5 |      11.0 |      22.2 |
    | hom k-mers: 1x       | 77.6 | 78.2 | 77.6 | 77.7 | 78.5 | 77.5 | 78.0 | 75.6 | 50.6 | 43.1 |      77.3 |      78.1 |      70.3 |      42.6 |
    | hom k-mers: 2x       | 13.7 | 13.0 | 13.2 | 13.2 | 12.6 | 13.3 | 12.4 | 11.7 | 27.8 | 40.8 |      13.6 |      12.2 |      18.4 |      34.5 |
    | hom k-mers:&nbsp;≥3x |  0.3 |  0.3 |  0.2 |  0.3 |  0.3 |  0.3 |  0.2 |  0.2 |  0.5 |  0.8 |       0.3 |       0.2 |       0.4 |       0.7 |

- <b><i>Haplotype 2</b></i>

    | Metric               | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.45 (l2) | 0.55 (l2) | 0.65 (l2) | 0.75 (l2) |
    | -------------------- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --------: | --------: | --------: | --------: |
    | het k-mers: 0x       | 45.8 | 45.4 | 45.2 | 45.6 | 46.0 | 46.3 | 45.5 | 43.5 | 39.2 | 44.2 |      46.2 |      45.6 |      48.2 |      36.7 |
    | het k-mers: 1x       | 53.8 | 54.2 | 54.4 | 54.0 | 53.7 | 53.3 | 54.1 | 56.0 | 59.4 | 54.6 |      53.5 |      54.0 |      51.3 |      61.8 |
    | het k-mers:&nbsp;≥2x |  0.4 |  0.4 |  0.4 |  0.4 |  0.4 |  0.4 |  0.4 |  0.5 |  1.4 |  1.2 |       0.4 |       0.4 |       0.5 |       1.5 |
    | hom k-mers: 0x       |  9.3 |  9.2 |  8.9 |  9.1 |  9.4 |  9.6 |  9.1 |  8.4 | 16.4 | 21.3 |       9.6 |       9.1 |      15.1 |      15.3 |
    | hom k-mers: 1x       | 78.3 | 77.9 | 78.5 | 78.5 | 78.7 | 78.8 | 78.2 | 76.0 | 50.4 | 49.0 |      78.2 |      78.1 |      71.1 |      47.4 |
    | hom k-mers: 2x       | 12.2 | 12.7 | 12.3 | 12.1 | 11.6 | 11.3 | 12.5 | 15.2 | 32.5 | 29.2 |      11.9 |      12.5 |      13.5 |      36.7 |
    | hom k-mers:&nbsp;≥3x |  0.3 |  0.3 |  0.3 |  0.3 |  0.3 |  0.2 |  0.3 |  0.3 |  0.6 |  0.5 |       0.2 |       0.3 |       0.3 |       0.7 |

- <b><i>Union (hap1 + hap2)</b></i>

    | Metric               | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.45 (l2) | 0.55 (l2) | 0.65 (l2) | 0.75 (l2) |
    | -------------------- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --------: | --------: | --------: | --------: |
    | het k-mers: 0x       |  0.6 |  0.6 |  0.6 |  0.6 |  0.6 |  0.6 |  0.6 |  0.6 |  0.6 |  0.6 |       0.6 |       0.6 |       0.6 |       0.6 |
    | het k-mers: 1x       | 88.5 | 88.6 | 88.8 | 88.9 | 89.6 | 89.6 | 89.4 | 89.6 | 80.6 | 75.8 |      89.0 |      89.6 |      89.6 |      75.5 |
    | het k-mers:&nbsp;≥2x | 10.9 | 10.8 | 10.7 | 10.5 |  9.8 |  9.8 | 10.0 |  9.8 | 18.9 | 23.7 |      10.4 |       9.8 |       9.8 |      23.9 |
    | hom k-mers: 0x       |  0.2 |  0.2 |  0.2 |  0.2 |  0.2 |  0.2 |  0.2 |  0.2 |  0.2 |  0.2 |       0.2 |       0.2 |       0.2 |       0.2 |
    | hom k-mers: 1x       |  4.1 |  4.1 |  4.2 |  4.2 |  4.2 |  4.3 |  4.2 |  4.3 |  3.7 |  3.6 |       4.2 |       4.2 |       4.3 |       3.5 |
    | hom k-mers: 2x       | 83.0 | 83.2 | 83.3 | 83.5 | 84.7 | 84.6 | 84.4 | 84.6 | 76.3 | 68.9 |      83.7 |      84.6 |      84.6 |      68.5 |
    | hom k-mers:&nbsp;≥3x | 12.7 | 12.5 | 12.3 | 12.1 | 10.9 | 11.0 | 11.2 | 11.0 | 19.8 | 27.4 |      11.9 |      11.0 |      11.0 |      27.8 |

### Quast
```
# change the working directory
cd 02_assembly/homcov55/

# run Quast
mkdir -p quast
for part in hap1 hap2; do
    files=(); labels=()
    for d in hifiasm_s0*/ hifiasm_l2_s0*/; do
        d=${d%/}; label=${d#hifiasm_}
        files+=("${d}/E_phylacis_${label}.asm.bp.${part}.p_ctg.fa")
        labels+=("$label")
    done
    quast.py --large --est-ref-size 534000000 \
        -o quast/${part} -l "$(IFS=,; echo "${labels[*]}")" "${files[@]}"
done

# change the working directory
cd ../..
```

- <b><i>Haplotype 1</b></i>

    ![Figure 03a. Quast contig viewer for haplotype 1.](./figures/fig03_hap1.png)

    | `l`  | `s`  |   N50 (bp) | L50 |  NG50 (bp) |   N90 (bp) | L90 |  NG90 (bp) |   auN (bp) | N's&nbsp;per&nbsp;100&nbsp;kb |
    | :--: | :--: | ---------: | --: | ---------: | ---------: | --: | ---------: | ---------: | ----------------------------: |
    |  3   | 0.30 | 43,772,942 |   6 | 51,317,246 | 17,090,617 |  12 | 33,590,266 | 43,842,519 |                         110.3 |
    |  3   | 0.35 | 51,317,246 |   5 | 51,317,246 | 20,641,050 |  12 | 34,632,737 | 45,079,273 |                         428.9 |
    |  3   | 0.40 | 51,330,377 |   5 | 51,330,377 | 17,090,617 |  12 | 33,387,292 | 44,105,774 |                         429.1 |
    |  3   | 0.45 | 41,957,526 |   6 | 51,317,246 | 17,090,617 |  12 | 33,387,292 | 44,313,252 |                         426.8 |
    |  3   | 0.50 | 51,330,377 |   5 | 51,330,377 | 20,449,059 |  12 | 34,632,737 | 45,320,782 |                         406.7 |
    |  3   | 0.55 | 41,957,526 |   6 | 51,359,987 | 18,186,819 |  12 | 33,662,322 | 44,075,504 |                         567.1 |
    |  3   | 0.60 | 39,841,891 |   6 | 41,957,526 | 20,394,891 |  13 | 28,751,543 | 41,292,460 |                         134.2 |
    |  3   | 0.65 | 39,841,891 |   6 | 41,957,526 | 20,394,891 |  12 | 28,751,543 | 41,997,972 |                         236.6 |
    |  3   | 0.70 | 51,351,936 |   6 | 54,056,018 | 29,899,960 |  12 | 36,525,735 | 45,920,878 |                         471.2 |
    |  3   | 0.75 | 51,317,246 |   6 | 56,292,571 | 20,394,891 |  14 | 36,525,735 | 45,842,407 |                         109.6 |
    |  2   | 0.45 | 41,957,526 |   6 | 51,317,246 | 18,177,790 |  12 | 33,387,292 | 44,085,185 |                         539.0 |
    |  2   | 0.55 | 39,841,891 |   6 | 43,799,538 | 18,352,263 |  13 | 28,751,543 | 41,751,729 |                         528.5 |
    |  2   | 0.65 | 41,957,526 |   6 | 51,359,987 | 28,751,543 |  12 | 33,549,482 | 44,094,347 |                         154.1 |
    |  2   | 0.75 | 51,144,740 |   6 | 51,317,246 | 17,320,306 |  13 | 34,632,737 | 44,526,071 |                         325.9 |

- <b><i>Haplotype 2</b></i>

    ![Figure 03b. Quast contig viewer for haplotype 2.](./figures/fig03_hap2.png)

    | `l`  | `s`  |   N50 (bp) | L50 |  NG50 (bp) |   N90 (bp) | L90 |  NG90 (bp) |   auN (bp) | N's&nbsp;per&nbsp;100&nbsp;kb |
    | :--: | :--: | ---------: | --: | ---------: | ---------: | --: | ---------: | ---------: | ----------------------------: |
    |  3   | 0.30 | 42,197,018 |   6 | 45,190,855 | 20,289,095 |  12 | 30,438,125 | 44,221,886 |                         264.7 |
    |  3   | 0.35 | 43,786,684 |   6 | 45,190,855 | 17,142,888 |  12 | 30,327,973 | 43,077,801 |                           0.0 |
    |  3   | 0.40 | 43,786,684 |   6 | 45,187,080 | 20,289,095 |  12 | 30,311,269 | 44,128,093 |                          30.3 |
    |  3   | 0.45 | 43,786,684 |   6 | 45,190,855 | 20,289,095 |  12 | 30,311,269 | 44,029,941 |                          30.6 |
    |  3   | 0.50 | 42,197,018 |   6 | 45,193,510 | 18,310,092 |  12 | 30,255,385 | 43,679,519 |                         117.5 |
    |  3   | 0.55 | 42,197,018 |   6 | 42,197,018 | 24,839,503 |  12 | 24,839,503 | 41,131,917 |                           0.0 |
    |  3   | 0.60 | 43,786,684 |   6 | 45,190,855 | 18,318,392 |  12 | 33,387,292 | 43,416,672 |                         447.7 |
    |  3   | 0.65 | 43,786,684 |   6 | 46,025,662 | 31,562,876 |  12 | 33,387,292 | 43,494,742 |                         525.0 |
    |  3   | 0.70 | 45,187,080 |   6 | 51,144,740 | 17,291,453 |  14 | 35,385,602 | 43,126,439 |                         729.4 |
    |  3   | 0.75 | 42,197,018 |   6 | 45,190,855 | 29,899,960 |  12 | 38,016,435 | 45,700,235 |                          11.1 |
    |  2   | 0.45 | 42,197,018 |   6 | 42,197,018 | 20,442,954 |  13 | 24,839,503 | 40,614,434 |                           0.0 |
    |  2   | 0.55 | 42,197,018 |   6 | 45,190,855 | 20,627,072 |  12 | 33,387,292 | 43,479,232 |                           0.0 |
    |  2   | 0.65 | 42,197,018 |   6 | 42,197,018 | 20,289,095 |  12 | 20,289,095 | 41,264,420 |                         437.7 |
    |  2   | 0.75 | 45,190,855 |   6 | 54,056,018 | 30,311,269 |  13 | 38,519,572 | 46,639,921 |                          40.7 |

### <a id="03e">Conclusion</a>
From all these tests, it seems that `-l 3 -s 0.50` has the most consistent results and is never an outlier compared to other parameters.

### Dotplot
```
# change the working directory
cd 02_assembly/homcov55/s050

# map both haplotypes
minimap2 -x asm20 -t 128 -N 1000 --secondary=no \
    E_phylacis_s050.asm.bp.hap1.p_ctg.fa \
    E_phylacis_s050.asm.bp.hap2.p_ctg.fa \
    > hap1_vs_hap2.paf

# generate the dotplot
Rscript scripts/pafCoordsDotPlotly.R \
    -i hap1_vs_hap2.paf \
    -o hap1_vs_hap2_dotplot2 -s -t

# change the working directory
cd ../../..
```
![Figure 04. Dotplot between hap1 and hap2 from s050.](./figures/fig04.png)

### Contamination check
```
# map the filtered ONT reads to the union
minimap2 -ax map-ont -t 128 --secondary=no \
    02_assembly/homcov55/hifiasm_s050/E_phylacis_s050.union.fa \
    01_qc/02_chopper/E_phylacis_filtered.fastq.gz \
    | samtools sort -@8 -o 03_qc/diploid.bam - && samtools index 03_qc/diploid.bam

# index each haplotype
samtools faidx 02_assembly/homcov55/hifiasm_s050/E_phylacis_s050.asm.bp.hap1.p_ctg.fa
samtools faidx 02_assembly/homcov55/hifiasm_s050/E_phylacis_s050.asm.bp.hap2.p_ctg.fa

# split the BAM file
samtools view -b 03_qc/diploid.bam $(cut -f1 02_assembly/homcov55/hifiasm_s050/E_phylacis_s050.asm.bp.hap1.p_ctg.fa.fai) > 03_qc/hap1.bam && samtools index 03_qc/hap1.bam
samtools view -b 03_qc/diploid.bam $(cut -f1 02_assembly/homcov55/hifiasm_s050/E_phylacis_s050.asm.bp.hap2.p_ctg.fa.fai) > 03_qc/hap2.bam && samtools index 03_qc/hap2.bam

# download the taxonomy metadata
mkdir -p 03_qc/taxdump && cd 03_qc/taxdump
curl -L https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -
cd ../..

# download core_nt database
mkdir -p 03_qc/blastdb && cd 03_qc/blastdb
update_blastdb.pl --decompress --num_threads 128 --source ncbi core_nt
update_blastdb.pl --decompress --source ncbi taxdb
cd ../..

# download UniProt / Swiss-Prot database
mkdir -p 03_qc/uniprot && cd 03_qc/uniprot
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
diamond makedb --in uniprot_sprot.fasta.gz -d uniprot_sprot \
    --taxonmap prot.accession2taxid.gz \
    --taxonnodes ../taxdump/nodes.dmp
cd ../..

# generate 100kb sliding-windows
seqkit sliding -s 100000 -W 100000 -g 02_assembly/homcov55/hifiasm_s050/E_phylacis_s050.asm.bp.hap1.p_ctg.fa > 03_qc/hap1.chunks.fa
seqkit sliding -s 100000 -W 100000 -g 02_assembly/homcov55/hifiasm_s050/E_phylacis_s050.asm.bp.hap2.p_ctg.fa > 03_qc/hap2.chunks.fa

# check windows against core_nt
blastn -task megablast -db 03_qc/blastdb/core_nt -query 03_qc/hap1.chunks.fa \
    -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 \
    -evalue 1e-25 -num_threads 128 \
    | sed 's/_sliding:[0-9]*-[0-9]*//' > 03_qc/hap1.blastn.out

blastn -task megablast -db 03_qc/blastdb/core_nt -query 03_qc/hap2.chunks.fa \
    -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 \
    -evalue 1e-25 -num_threads 128 \
    | sed 's/_sliding:[0-9]*-[0-9]*//' > 03_qc/hap2.blastn.out

# check windows against UniProt / Swiss-Prot
diamond blastx --query 03_qc/hap1.chunks.fa --db 03_qc/uniprot/uniprot_sprot.dmnd \
    --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 128 \
    | sed 's/_sliding:[0-9]*-[0-9]*//' > 03_qc/hap1.diamond.out

diamond blastx --query 03_qc/hap2.chunks.fa --db 03_qc/uniprot/uniprot_sprot.dmnd \
    --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 128 \
    | sed 's/_sliding:[0-9]*-[0-9]*//' > 03_qc/hap2.diamond.out

# build the Blobdir
blobtools create --fasta E_phylacis_s050.asm.bp.hap1.p_ctg.fa --taxid 3932 --taxdump 03_qc/taxdump/ 03_qc/hap1_btk
blobtools add --cov 03_qc/hap1.bam --threads 16 03_qc/hap1_btk
blobtools add --hits 03_qc/hap1.blastn.out --hits 03_qc/hap1.diamond.out \
    --taxrule bestsumorder --taxdump taxdump/ 03_qc/hap1_btk
blobtools add --busco 02_assembly/homcov55/compleasm/s050_hap1/eudicotyledons_odb12/full_table.tsv 03_qc/hap1_btk

blobtools create --fasta E_phylacis_s050.asm.bp.hap2.p_ctg.fa --taxid 3932 --taxdump 03_qc/taxdump/ 03_qc/hap2_btk
blobtools add --cov 03_qc/hap2.bam --threads 16 03_qc/hap2_btk
blobtools add --hits 03_qc/hap2.blastn.out --hits 03_qc/hap2.diamond.out \
    --taxrule bestsumorder --taxdump taxdump/ 03_qc/hap2_btk
blobtools add --busco 02_assembly/homcov55/compleasm/s050_hap2/eudicotyledons_odb12/full_table.tsv 03_qc/hap2_btk

# export the output
blobtools filter --table 03_qc/hap1_btk.tsv hap1_btk
blobtools filter --table 03_qc/hap2_btk.tsv hap2_btk

blobtools view \
    --view blob \
    --format png \
    --out 03_qc/fig_hap1 \
    03_qc/hap1_btk

blobtools view \
    --view blob \
    --format png \
    --out 03_qc/fig_hap2 \
    03_qc/hap2_btk
```

## <a id="04">Scaffolding</a>

---
*Last update: 05 November 2026 by Jeremias Ivan*
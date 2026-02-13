# An assembly of the Meelup Mallee genome.

# Get the environment running

I'll try and do most things in one conda environment. The details of what that looks like are in the `environment.yml` file.

```bash
conda env create environment.yml
```


# Raw Data

The raw data are located here:

```bash
raw_data="/home/raw_data/MM/2025_long_reads/ONT_gDNA192_1071_RL/E_Phylacis/20251117_1227_2C_PBE85256_e040940c/fastq_pass"
```

Let's get a first impression of how much data there is:

```bash
seqkit stats -j 64 -T ${raw_data}/*.fastq.gz > raw_data_seqkit_stats.tsv
grep -v "^file" raw_data_seqkit_stats.tsv | sed 's/,//g' | awk -F'\t' '{r+=$4; b+=$5} END {printf "Reads: %'\''d | Bases: %'\''d | Coverage: %.2fx\n", r, b, b/500000000}'
```

* **Reads**: 3,064,194 
* **Bases**: 44,326,974,907 
* **Coverage**: 88.65x

This shows that we have ~90x coverage (~45 of each haplotype) before QC and filtering, so a good place to start. This is based on an estiamted 500MB genome size.

# 01 QC

First let's examine the raw long reads carefully.

## Basic QC

```bash
qc_dir="01_QC"
mkdir ${qc_dir}

# 1. Run NanoPlot
NanoPlot -t 128 \
         --fastq ${raw_data}/*.fastq.gz \
         --downsample 100000 \
         -o ${qc_dir}/01_NanoPlot_Raw \
         --title "E_phylacis_Raw_ONT"

# 2. Long-read K-mer Counting
## set up directories
mkdir -p ${qc_dir}/02_Kmer_distribution
tmp_dir="tmp_processing"
mkdir -p ${tmp_dir}

## list the raw data files
find ${raw_data} -name "*.fastq.gz" > ${tmp_dir}/files.txt

## run KMC
# -k21: Standard k-mer length for GenomeScope 2.0
# -t128: Using all 128 threads for speed
# -m256: 256GB RAM limit
# -ci1: Include k-mers that occur at least once
kmc -k21 -t128 -m256 -ci1 -cs10000 \
    @${tmp_dir}/files.txt \
    ${tmp_dir}/kmc_db \
    ${tmp_dir}/

## create histogram
# keep everything down to 1 occurrence; bunch up the stuff that occurrs >10K times
kmc_tools transform ${tmp_dir}/kmc_db histogram ${qc_dir}/02_Kmer_distribution/long_read_histogram.txt -ci1 -cs1000000

# Run GenomeScope2
# -p : diploid
genomescope2 \
    -i ${qc_dir}/02_Kmer_distribution/long_read_histogram.txt \
    -o ${qc_dir}/02_Kmer_distribution/genomescope_results \
    -k 21 \
    -p 2

rm -rf ${tmp_dir}

```

Nanoplot basics

| Metric | Value |
| :--- | :--- |
| **Total Yield** | 44.3 Gb |
| **Number of Reads** | 3.06 M |
| **Read Length N50** | 22.2 kb |
| **Mean Read Length** | 14.5 kb |
| **Median Read Length** | 12.3 kb |
| **Mean Read Quality** | Q17.2 |
| **Median Read Quality** | Q19.6 (~99%) |
| **> Q10** (90.0% accuracy) | 44.3 Gb |
| **> Q15** (96.8% accuracy) | 38.2 Gb |
| **> Q20** (99.0% accuracy) | 22.2 Gb |

![Read Length vs Quality](01_QC/01_NanoPlot_Raw/LengthvsQualityScatterPlot_kde.png)

KMC basics output
```
Stats:
   No. of k-mers below min. threshold :            0
   No. of k-mers above max. threshold :            0
   No. of unique k-mers               :   5216090955
   No. of unique counted k-mers       :   5216090955
   Total no. of k-mers                :  44265691027
   Total no. of reads                 :      3064194
   Total no. of super-k-mers          :   6430181221
```

Genomescope:

![GenomeScope K-mer Plot](01_QC/02_Kmer_distribution/genomescope_results/transformed_linear_plot.png)


So we have about 4% heterozygosity, a genome size of ~523MB, and a low error rate of 0.6%. Quite nice!

## Contamaination checking

A quick check of the GC content suggests no serious contamination. 

```bash
zcat ${raw_data}/*.fastq.gz | head -n 400000 | \
seqkit fx2tab -n -g | cut -f 2 | \
awk '{printf "%.0f\n", $1}' | sort | uniq -c | \
awk '{printf "%s%%\t%s\n", $2, $1}' | sort -n
```

```
0%      1
9%      1
10%     1
11%     1
13%     1
15%     1
16%     1
17%     5
18%     5
19%     7
20%     6
21%     10
22%     17
23%     22
24%     32
25%     38
26%     67
27%     136
28%     145
29%     238
30%     298
31%     470
32%     690
33%     1117
34%     1895
35%     3027
36%     5375
37%     8983
38%     13473
39%     15529
40%     15264
41%     11650
42%     7815
43%     4459
44%     2612
45%     1632
46%     980
47%     636
48%     516
49%     387
50%     342
51%     228
52%     219
53%     257
54%     303
55%     442
56%     203
57%     85
58%     64
59%     43
60%     45
61%     22
62%     31
63%     24
64%     25
65%     13
66%     25
67%     17
68%     12
69%     11
70%     12
71%     11
72%     4
73%     6
74%     3
75%     3
76%     2
77%     3
79%     1
100%    1
```

This shows a unimodal distribution around 38-41% GC, so it's best to just assemble the genome first and decontaminate contigs later.



# 02 Read filtering

Let's filter the data with Chopper.

```bash
# 1. Setup Directories
filter_dir="02_filtering"
mkdir -p ${filter_dir}/qc_filtered

# let's not get ourselves in trouble
echo "02_filtering/E_phylacis_filtered.fastq.gz" >> .gitignore


# 2. Run Chopper
# Pipe pigz (decompression) -> chopper (filtering) -> bgzip (compression)
echo "Starting Chopper: Filtering for Length > 15kb and Quality > Q10..."
pigz -dc -p 128 ${raw_data}/*.fastq.gz | \
chopper -q 10 -l 15000 | \
bgzip -@ 128 > ${filter_dir}/E_phylacis_filtered.fastq.gz
echo "Filtering complete. Output saved to: ${filter_dir}/E_phylacis_filtered.fastq.gz"

# 3. Post-Filter QC (NanoPlot)
NanoPlot -t 128 \
         --fastq ${filter_dir}/E_phylacis_filtered.fastq.gz \
         --downsample 100000 \
         -o ${filter_dir}/qc_filtered \
         --title "E_phylacis_Filtered_ONT"
```

Post filtering QC:

| Metric | Pre-Filtering | Post-Filtering |
| :--- | :--- | :--- |
| **Total Yield** | 44.3 Gb | **34.0 Gb** |
| **Number of Reads** | 3.06 M | **1.35 M** |
| **Read Length N50** | 22.2 kb | **25.6 kb** |
| **Mean Read Length** | 14.5 kb | **25.2 kb** |
| **Median Read Length** | 12.3 kb | **22.8 kb** |
| **Mean Read Quality** | Q17.2 | **Q17.8** |
| **Median Read Quality** | Q19.6 (~99%) | **Q20.2** |
| **> Q10** (90.0% accuracy) | 44.3 Gb | **34.0 Gb** |
| **> Q15** (96.8% accuracy) | 38.2 Gb | **29.3 Gb** |
| **> Q20** (99.0% accuracy) | 22.2 Gb | **17.2 Gb** |

![Read Length vs Quality](02_filtering/qc_filtered/LengthvsQualityScatterPlot_kde.png)


This is great. We still have >30x coverage per haplotype, median quality is >Q20, and average read length is well over 20KB. Over half the data is >Q20 reads too. 

# 03 basic hifiasm assembly

Now I'll assemble the filtered reads with hifiasm. to do this I'll make a Ramdisk to do everything in RAM. I have 2.2TB so this should work OK...


```bash
mkdir -p 03_hifiasm_assembly

# 1. Create a mount point
sudo mkdir -p /mnt/ramdisk

# 2. Mount 1500GB of the 2.2TB RAM as a disk
sudo mount -t tmpfs -o size=1500G tmpfs /mnt/ramdisk
sudo chown $USER /mnt/ramdisk

# 3. Copy your filtered data THERE
cp 02_filtering/E_phylacis_filtered.fastq.gz /mnt/ramdisk/

# 4. Run hifiasm inside the ramdisk
cd /mnt/ramdisk

hifiasm \
    -o E_phylacis_asm \
    -t 160 \
    --ont \
    -l 3 \
    --telo-m AAACCCT \
    --dual-scaf \
    E_phylacis_filtered.fastq.gz \
    2>&1 | tee hifiasm_restart.log

# 5. VERY IMPORTANT: Move results back to /home before rebooting!
rsync -av --progress --exclude='*.fastq.gz' ./ ~/MM_genome/03_hifiasm_assembly/

cd ~/MM_genome

# ignore this stuff for git, mostly
echo "03_hifiasm_assembly/" >> .gitignore

# but keep these guys
git add -f 03_hifiasm_assembly/hifiasm_restart.log
git add -f 03_hifiasm_assembly/*.noseq.gfa
git add -f 03_hifiasm_assembly/*.lowQ.bed
git commit -m "Ignore raw assembly binaries but keep logs and noseq graphs"

```


## Assembly summary

Let's summarise the graphs with gfastats

```bash
mkdir -p 03_hifiasm_assembly/QC
gfastats 03_hifiasm_assembly/E_phylacis_asm.bp.hap1.p_ctg.gfa --discover-paths --segment-report > 03_hifiasm_assembly/QC/stats_hap1_segments.txt
gfastats 03_hifiasm_assembly/E_phylacis_asm.bp.hap1.p_ctg.gfa --discover-paths > 03_hifiasm_assembly/QC/stats_hap1.txt

gfastats 03_hifiasm_assembly/E_phylacis_asm.bp.hap2.p_ctg.gfa --discover-paths --segment-report > 03_hifiasm_assembly/QC/stats_hap2_segments.txt
gfastats 03_hifiasm_assembly/E_phylacis_asm.bp.hap2.p_ctg.gfa --discover-paths > 03_hifiasm_assembly/QC/stats_hap2.txt

gfastats 03_hifiasm_assembly/E_phylacis_asm.bp.p_ctg.gfa --discover-paths --segment-report > 03_hifiasm_assembly/QC/stats_primary_segments.txt
gfastats 03_hifiasm_assembly/E_phylacis_asm.bp.p_ctg.gfa --discover-paths > 03_hifiasm_assembly/QC/stats_primary.txt

git add -f 03_hifiasm_assembly/QC/
```

Hap1:
* Length: 579570136
* N50: 39833816
* Sum of top 11: 479,816,686 (82.7% of total)

Hap2:
* Length: 580987427
* N50: 43786601
* Sum of top 11: 466,703,285 (80.3% of total)

Primary contig assembly:
* Length: 545423967
* N50: 45651274
* Sum of top 11: 512,180,968 (93.9% of total)

Hap1 top 11 scaffolds:
```
Seq     Header  Comment Total segment length    A       C       G       T       GC content %    # soft-masked bases     Is circular: 
1       h1tg000001l             51317275        15581784        10088548        10060793        15586150        39.26   0       N
2       h1tg000002l             34632735        10403692        6918163 6905508 10405372        39.92   0       N
3       h1tg000003l             59377836        18014668        11731310        11638807        17993051        39.36   0       N
4       h1tg000004l             38027376        11524382        7473569 7506629 11522796        39.39   0       N
5       h1tg000005l             39657523        11954888        7853657 7859159 11989819        39.62   0       N
6       h1tg000006l             28751543        8692468 5642168 5679416 8737491 39.38   0       N
7       h1tg000007l             41957522        12681753        8316763 8322106 12636900        39.66   0       N
8       h1tg000008l             61728905        18653045        12184965        12191786        18699109        39.49   0       N
9       h1tg000009l             39833816        12005760        7924497 7892044 12011515        39.71   0       N
10      h1tg000010l             33387284        10075060        6608064 6625293 10078867        39.64   0       N
11      h1tg000011l             51144871        15434505        10210279        10161001        15339086        39.83   0       N
```

Hap2 top 11 scaffolds:

```
Seq     Header  Comment Total segment length    A       C       G       T       GC content %    # soft-masked bases     Is circular: 
1       h2tg000001l             54055996        16246510        10739895        10777096        16292495        39.81   0       N
2       h2tg000002l             45187061        13709758        8862292 8891714 13723297        39.29   0       N
3       h2tg000003l             65679588        19845392        12967282        12971991        19894923        39.49   0       N
4       h2tg000004l             38512735        11580937        7688737 7644481 11598580        39.81   0       N
5       h2tg000005l             18271583        5236722 3476365 3493002 5220136 39.99   0       N
6       h2tg000006l             36524851        11017298        7238859 7254846 11013848        39.68   0       N
7       h2tg000007l             42203136        12767072        8362376 8356582 12717106        39.62   0       N
8       h2tg000008l             30255392        9129461 5979442 5996696 9149793 39.58   0       N
9       h2tg000009l             56292537        17057888        11068887        11069636        17096126        39.33   0       N
10      h2tg000010l             43786601        13237045        8685725 8676979 13186752        39.65   0       N
11      h2tg000011l             35933805        10838134        7021993 7082002 10855327        39.40   0       N
```

Primary contig (bp.p_ctg) top 11 scaffolds
```
Seq     Header  Comment Total segment length    A       C       G       T       GC content %    # soft-masked bases     Is circular: 
1       ptg000001l              45651274        13837378        8968599 8979051 13866246        39.31   0       N
2       ptg000002l              54055996        16246510        10739895        10777096        16292495        39.81   0       N
3       ptg000003l              59271212        17985060        11714462        11612997        17958693        39.36   0       N
4       ptg000004l              38027376        11524382        7473569 7506629 11522796        39.39   0       N
5       ptg000005l              61704420        18645782        12179659        12186397        18692582        39.49   0       N
6       ptg000006l              38521672        11582632        7689359 7647684 11601997        39.81   0       N
7       ptg000007l              42285308        12737509        8368553 8385101 12794145        39.62   0       N
8       ptg000008l              41957522        12681753        8316763 8322106 12636900        39.66   0       N
9       ptg000009l              57093733        17310108        11249947        11235059        17298619        39.38   0       N
10      ptg000010l              39824047        12001580        7920315 7892040 12010112        39.71   0       N
11      ptg000011l              33788408        10173370        6719776 6724548 10170714        39.79   0       N
```

## Telomere checks

Let's check for common repeats, and see where they are.

First we explore for common repeats:

```bash
tidk explore --length 7 --minimum 5 --maximum 12 03_hifiasm_assembly/E_phylacis_asm.bp.p_ctg.fa
tidk explore --length 7 --minimum 5 --maximum 12 03_hifiasm_assembly/E_phylacis_asm.bp.hap1.p_ctg.fa
tidk explore --length 7 --minimum 5 --maximum 12 03_hifiasm_assembly/E_phylacis_asm.bp.hap1.p_ctg.fa
```

All looks good:

```
canonical_repeat_unit   count_repeat_runs_gt_100
AAACCCT 8816
ACCCGTC 2829
[+]     Exploring genome for potential telomeric repeats of length: 7
[+]     Finished searching genome
[+]     Generating output
canonical_repeat_unit   count_repeat_runs_gt_100
AAACCCT 8816
AAAAAAG 1157
AAGACTC 635
[+]     Exploring genome for potential telomeric repeats of length: 7
[+]     Finished searching genome
[+]     Generating output
canonical_repeat_unit   count_repeat_runs_gt_100
AAACCCT 8816
AAAAAAG 1157
AAGACTC 635
```

This is good! Exactly what we expect. The AAACCCT is the telomere. The others are likely centromeres, LINEs, or similar. Let's see where they are on the scaffolds (I'll focus on the first 11 scaffolds of each assembly, becuase that's most of it). 


First we'll use `tidk search` to figure out where they are in each of the three assemblies:

```bash
# Define motifs and assemblies
MOTIFS=("AAACCCT" "ACCCGTC" "AAAAAAG" "AAGACTC")
FILES=("p_ctg" "hap1.p_ctg" "hap2.p_ctg")

mkdir -p 03_hifiasm_assembly/QC/tidk_plots

for f in "${FILES[@]}"; do
  for m in "${MOTIFS[@]}"; do
    echo "Searching for $m in $f..."
    tidk search --string $m \
      --dir 03_hifiasm_assembly/QC/tidk_plots \
      --output "${f}_${m}" \
      --extension tsv \
      "03_hifiasm_assembly/E_phylacis_asm.bp.${f}.fa"
  done
done
```

| Primary Assembly | Haplotype 1 | Haplotype 2 |
| :---: | :---: | :---: |
| [![Primary](03_hifiasm_assembly/QC/telomere_results/p_ctg_repeat_fingerprint.png)](03_hifiasm_assembly/QC/telomere_results/p_ctg_repeat_fingerprint.png) | [![Hap1](03_hifiasm_assembly/QC/telomere_results/hap1_repeat_fingerprint.png)](03_hifiasm_assembly/QC/telomere_results/hap1_repeat_fingerprint.png) | [![Hap2](03_hifiasm_assembly/QC/telomere_results/hap2_repeat_fingerprint.png)](03_hifiasm_assembly/QC/telomere_results/hap2_repeat_fingerprint.png) |
| *Click to enlarge* | *Click to enlarge* | *Click to enlarge* |

**Legend:**
* **Red (AAACCCT):** Canonical plant telomere motif.
* **Sky Blue (ACCCGTC):** Putative centromeric satellite.
* **Green (AAAAAAG):** Transposon-associated / Poly-A repeats.
* **Purple (AAGACTC):** Secondary satellite motif.
* Black points represent raw 10kb window counts (alpha 0.5); colored lines indicate smoothed trends.

This is pretty good. Most of the 11 big scaffolds are T2T in the primary and the two haplotype assemblies.

## BUSCO with Compleasm

First I'll make versions of the assemblies with just the 11 longest contigs, since the telomere analysis suggests that this might be quite complete.

```bash
# Define your assemblies
ASSEMBLIES=(
    "E_phylacis_asm.bp.p_ctg.fa"
    "E_phylacis_asm.bp.hap1.p_ctg.fa"
    "E_phylacis_asm.bp.hap2.p_ctg.fa"
)

THREADS=64

for FASTA in "${ASSEMBLIES[@]}"; do
    # Create a base name (e.g., p_ctg, hap1, hap2)
    BASE=$(echo "$FASTA" | sed 's/E_phylacis_asm.bp.//; s/.p_ctg.fa//; s/.fa//')
    
    echo "-------------------------------------------------------"
    echo "Processing: $BASE"
    echo "-------------------------------------------------------"

    # 1. Index the file
    samtools faidx 03_hifiasm_assembly/$FASTA

    # 2. Get the names of the top 11 longest contigs
    cut -f1,2 03_hifiasm_assembly/${FASTA}.fai | \
        sort -k2,2rn | \
        head -n 11 | \
        cut -f1 > ${BASE}_top11_list.txt

    # 3. Extract these contigs
    samtools faidx 03_hifiasm_assembly/$FASTA \
        -r ${BASE}_top11_list.txt > 03_hifiasm_assembly/E_phylacis_${BASE}_top11.fa

    rm ${BASE}_top11_list.txt
done
```


Let's look at the three assemblies completeness with BUSCOs. We'll do general and specific, using embryophyta_odb12 and eudicots_odb12 respectively. We'll repeat it for the full assembly and the 11 biggest contigs.

```bash
# Configuration
OUT_BASE="03_hifiasm_assembly/QC/compleasm_results"
THREADS=128

mkdir -p "$OUT_DIR"

# Download the databases (compleasm handles the URLs automatically)
compleasm download embryophyta
compleasm download eudicotyledons

# --- PRIMARY ASSEMBLY ---
compleasm run -a 03_hifiasm_assembly/E_phylacis_asm.bp.p_ctg.fa -o ${OUT_BASE}/p_ctg_embryo -l embryophyta -t $THREADS
compleasm run -a 03_hifiasm_assembly/E_phylacis_asm.bp.p_ctg.fa -o ${OUT_BASE}/p_ctg_eudicot -l eudicotyledons -t $THREADS

compleasm run -a 03_hifiasm_assembly/E_phylacis_p_ctg_top11.fa -o ${OUT_BASE}/p_ctg_top11_embryo -l embryophyta -t $THREADS
compleasm run -a 03_hifiasm_assembly/E_phylacis_p_ctg_top11.fa -o ${OUT_BASE}/p_ctg_top11_eudicot -l eudicotyledons -t $THREADS

# --- HAPLOTYPE 1 ---
compleasm run -a 03_hifiasm_assembly/E_phylacis_asm.bp.hap1.p_ctg.fa -o ${OUT_BASE}/hap1_embryo -l embryophyta -t $THREADS
compleasm run -a 03_hifiasm_assembly/E_phylacis_asm.bp.hap1.p_ctg.fa -o ${OUT_BASE}/hap1_eudicot -l eudicotyledons -t $THREADS

compleasm run -a 03_hifiasm_assembly/E_phylacis_hap1_top11.fa -o ${OUT_BASE}/hap1_top11_embryo -l embryophyta -t $THREADS
compleasm run -a 03_hifiasm_assembly/E_phylacis_hap1_top11.fa -o ${OUT_BASE}/hap1_top11_eudicot -l eudicotyledons -t $THREADS

# --- HAPLOTYPE 2 ---
compleasm run -a 03_hifiasm_assembly/E_phylacis_asm.bp.hap2.p_ctg.fa -o ${OUT_BASE}/hap2_embryo -l embryophyta -t $THREADS
compleasm run -a 03_hifiasm_assembly/E_phylacis_asm.bp.hap2.p_ctg.fa -o ${OUT_BASE}/hap2_eudicot -l eudicotyledons -t $THREADS

compleasm run -a 03_hifiasm_assembly/E_phylacis_hap2_top11.fa -o ${OUT_BASE}/hap2_top11_embryo -l embryophyta -t $THREADS
compleasm run -a 03_hifiasm_assembly/E_phylacis_hap2_top11.fa -o ${OUT_BASE}/hap2_top11_eudicot -l eudicotyledons -t $THREADS

```

Once that's done, we can get the table below using

```bash
python scripts/summarise_compleasm.py
```

| Dataset with embryophyta |      S (%, (N)) |      D (%, (N)) |      F (%, (N)) |      M (%, (N)) |    Total |
| ---------------------- | ---------------: | ---------------: | ---------------: | ---------------: | --------: |
| p_ctg_embryo           |   96.50% (1955) |      3.11% (63) |       0.39% (8) |       0.00% (0) |     2026 |
| p_ctg_top11_embryo     |   94.32% (1911) |      2.62% (53) |      0.59% (12) |      2.47% (50) |     2026 |
| hap1_embryo            |   91.31% (1850) |     7.90% (160) |       0.39% (8) |       0.39% (8) |     2026 |
| hap1_top11_embryo      |   91.51% (1854) |      2.42% (49) |      0.79% (16) |     5.28% (107) |     2026 |
| hap2_embryo            |   93.14% (1887) |     6.52% (132) |       0.35% (7) |       0.00% (0) |     2026 |
| hap2_top11_embryo      |   93.48% (1894) |      1.97% (40) |      0.59% (12) |      3.95% (80) |     2026 |


| Dataset with eudicots  |      S (%, (N)) |      D (%, (N)) |      F (%, (N)) |      M (%, (N)) |    Total |
| ---------------------- | ---------------: | ---------------: | ---------------: | ---------------: | --------: |
| p_ctg_eudicot          |   96.79% (2715) |      2.89% (81) |       0.29% (8) |       0.04% (1) |     2805 |
| p_ctg_top11_eudicot    |   95.04% (2666) |      2.14% (60) |      0.53% (15) |      2.28% (64) |     2805 |
| hap1_eudicot           |   91.76% (2574) |     7.31% (205) |      0.39% (11) |      0.53% (15) |     2805 |
| hap1_top11_eudicot     |   91.94% (2579) |      2.00% (56) |      0.75% (21) |     5.31% (149) |     2805 |
| hap2_eudicot           |   92.66% (2599) |     7.06% (198) |       0.25% (7) |       0.04% (1) |     2805 |
| hap2_top11_eudicot     |   93.69% (2628) |      2.03% (57) |      0.57% (16) |     3.71% (104) |     2805 |

A few observations. The primary is exceptionally complete, and even has 97.5% of BUSCOs (either set) when you just take the top 11 contigs. The smaller contigs contain quite a few duplicates.

The two haplotype assemblies are also very complete, and dropping the small contigs massively drops the duplication rate, at the cost of putting the missing % up to ~4-5%. 

In terms of the top11 assemblies which have just the 11 'chromosome' contigs, completeness is still ~95% or higher. The duplication rate is ~2%, and there's ~1% fragmented.

This is good, showing we can stick with these assemblies from here. 

## Comparing hap1 and hap2

First let's figure out which contigs are which in the haplotypes. We'll align it and dot-plot it.

```bash
mkdir 03_hifiasm_assembly/compare_h1h2

# Align Hap1 to Hap2
minimap2 -x asm5 -t 128 -N 1000 --secondary=no \
    03_hifiasm_assembly/E_phylacis_hap1_top11.fa \
    03_hifiasm_assembly/E_phylacis_hap2_top11.fa \
    > 03_hifiasm_assembly/compare_h1h2/hap1_vs_hap2.paf

Rscript scripts/pafCoordsDotPlotly.R   -i 03_hifiasm_assembly/compare_h1h2/hap1_vs_hap2.paf   -o hap1_vs_hap2_dotplot2 -s -t -m 2000 -q 500000

mv hap1_vs_hap2_dotplot2.* 03_hifiasm_assembly/compare_h1h2/.
```

The dotplot shows good contiguity between them, but the % identity is low (gaps and repeats??)
![Hap1 vs Hap2 Dotplot](03_hifiasm_assembly/compare_h1h2/hap1_vs_hap2_dotplot2.png)

We can look at coverage vs. identity using the `dv` column of the PAF, which is similar to P distance.


```bash
awk 'BEGIN {
    print "| H2 ID | H1 ID | Total Matches | H2 Coverage | Identity (1-dv) |";
    print "| :--- | :--- | :--- | :--- | :--- |";
}
{
    # Extract dv tag
    dv=0; for(i=13;i<=NF;i++) if($i ~ /^dv:f:/) dv=substr($i,6);
    
    # Key is Query + Target pair
    key=$1"|"$6;
    matches[key]+=$10;
    # Store target length (Col 7) to calculate coverage later
    t_lengths[key]=$7;
    # Store the dv for the pair
    dvs[key]=dv; 
}
END {
    for (k in matches) {
        split(k, names, "|");
        coverage = (matches[k] / t_lengths[k]) * 100;
        identity = (1 - dvs[k]) * 100;
        
        # Filter for significant matches only (optional, e.g., > 1% coverage)
        # to avoid printing the repetitive noise again
        if (coverage > 1) {
            printf "| %s | %s | %d | %.2f%% | %.2f%% |\n", \
            names[1], names[2], matches[k], coverage, identity;
        }
    }
}' 03_hifiasm_assembly/compare_h1h2/hap1_vs_hap2.paf | sort -t'|' -k2,2 -k4,4nr | awk -F'|' '!seen[$2]++'
```

| H2 ID | H1 ID | Total Matches | H2 Coverage | Identity (1-dv) |
| :--- | :--- | :--- | :--- | :--- |
| h2tg000001l | h1tg000011l | 18993427 | 37.14% | 97.35% |
| h2tg000002l | h1tg000001l | 13955987 | 27.20% | 94.97% |
| h2tg000003l | h1tg000008l | 20029949 | 32.45% | 96.42% |
| h2tg000004l | h1tg000002l | 13101698 | 37.83% | 89.83% |
| h2tg000006l | h1tg000009l | 12991176 | 32.61% | 92.96% |
| h2tg000007l | h1tg000005l | 14464822 | 36.47% | 96.42% |
| h2tg000008l | h1tg000010l | 10033092 | 30.05% | 93.99% |
| h2tg000009l | h1tg000006l | 8494638 | 29.54% | 97.89% |
| h2tg000010l | h1tg000007l | 14004298 | 33.38% | 92.96% |
| h2tg000011l | h1tg000004l | 11697283 | 30.76% | 94.38% |
| h2tg000307l | h1tg000003l | 15749157 | 26.52% | 94.97% |

This table gives the correspondence between the chromosomes in h1 and h2. It shows that the coverage is low (as shown in the dotplot) due to indels and structural variants. But where the sequences align, they are very similar (89-98% identity). Chr 4l from h2 and 2l from h1 is interesting - they match by the lowest amount of 89%, but the dotplot for these looks the cleanest! So they are highly sytenic, but with high sequence divergence. 


# 04 Comparison with decipiens and virginea

Next we want to compare the new genomes to those of the putative parental species (not individuals, we don't have the individual parents). These are E. decipiens and E. virginea.

First we get the genomes:

```bash
# Create the root folder
mkdir -p parental_spp_genomes
cd parental_spp_genomes

# Download E. decipiens (GCA_014182575.1)
datasets download genome accession GCA_014182575.1 --filename e_decipiens.zip
unzip e_decipiens.zip -d e_decipiens_dir
mv e_decipiens_dir/ncbi_dataset/data/GCA_014182575.1/*.fna E_decipiens.fa

# Download E. virginea (GCA_014182375.1)
datasets download genome accession GCA_014182375.1 --filename e_virginea.zip
unzip e_virginea.zip -d e_virginea_dir
mv e_virginea_dir/ncbi_dataset/data/GCA_014182375.1/*.fna E_virginea.fa

# Cleanup
rm -rf *_dir *.zip
cd ..
```
## Via genome alignment

Now we align them all to the two parents. I'll cut out secondary alignments for now.

```bash
# Create output directory for alignments
mkdir -p 04_parental_assignment

# Alignment of Hap1 vs both parents
minimap2 -x asm5 -N 1000 --secondary=no -t 128 parental_spp_genomes/E_decipiens.fa 03_hifiasm_assembly/E_phylacis_hap1_top11.fa > 04_parental_assignment/hap1_vs_decipiens.paf
minimap2 -x asm5 -N 1000 --secondary=no -t 128 parental_spp_genomes/E_virginea.fa 03_hifiasm_assembly/E_phylacis_hap1_top11.fa > 04_parental_assignment/hap1_vs_virginea.paf

# Alignment of Hap2 vs both parents
minimap2 -x asm5 -N 1000 --secondary=no -t 128 parental_spp_genomes/E_decipiens.fa 03_hifiasm_assembly/E_phylacis_hap2_top11.fa > 04_parental_assignment/hap2_vs_decipiens.paf
minimap2 -x asm5 -N 1000 --secondary=no -t 128 parental_spp_genomes/E_virginea.fa 03_hifiasm_assembly/E_phylacis_hap2_top11.fa > 04_parental_assignment/hap2_vs_virginea.paf
```

A rough look at these alignments is to get the identity to decipiens and virginea across the whole chromosomes of hap1 and hap2, using the approach we did above for the hap1 v hap2 paf. 

```bash
awk 'BEGIN {
    print "| S1 ID | S2 ID | Total Matches | S1 Coverage | Identity |";
    print "| :--- | :--- | :--- | :--- | :--- |";
}
{
    # Extract dv tag
    dv=0; for(i=13;i<=NF;i++) if($i ~ /^dv:f:/) dv=substr($i,6);
    
    # Key is Query + Target pair
    key=$1"|"$6;
    matches[key]+=$10;
    # Store target length (Col 7) to calculate coverage later
    t_lengths[key]=$7;
    # Store the dv for the pair
    dvs[key]=dv; 
}
END {
    for (k in matches) {
        split(k, names, "|");
        coverage = (matches[k] / t_lengths[k]) * 100;
        identity = (1 - dvs[k]) * 100;
        
        # Filter for significant matches only (optional, e.g., > 1% coverage)
        # to avoid printing the repetitive noise again
        if (coverage > 1) {
            printf "| %s | %s | %d | %.2f%% | %.2f%% |\n", \
            names[1], names[2], matches[k], coverage, identity;
        }
    }
}' 03_hifiasm_assembly/compare_h1h2/hap1_vs_hap2.paf | sort -t'|' -k2,2 -k4,4nr | awk -F'|' '!seen[$2]++'
```

In the following tables - I ran the above twice for each of h1 and h2 (vs the two parents). I then organised the rows according to the decipiens chromosome, so row 1 corresponds to homologous chromosomes in both tables, etc. I bolded the parent with the highest identity for each row.

### Haplotype 1 Assignment
| S1 ID | decip ID | virgin ID | decip Matches | virgin Matches | decip Coverage | virgin Coverage | decip Identity | virgin Identity |
| :---  | :---     | :---      | :---          | :---           | :---           | :---            | :---           | :---            |
| h1tg000001l | CM024615.1 |  CM024520.1 | 24507427 |  14259343 | **40.71%** | 24.90% | **97.35%** |  92.96% |
| h1tg000002l | CM024618.1 |  CM024523.1 | 12497885 |  21930450 | 32.71% | **59.56%** | **97.35%** |  91.27% |
| h1tg000003l | CM024613.1 |  CM024518.1 | 15221461 |  26814450 | 24.60% | **47.00%** | 93.86% |  **97.35%** |
| h1tg000004l | CM024612.1 |  CM024517.1 | 21161217 |  10238986 | **49.09%** | 31.61% | 88.95% |  **93.39%** |
| h1tg000005l | CM024619.1 |  CM024524.1 | 14717449 |  25884302 | 31.37% | **54.94%** | 97.35% |  **97.89%** |
| h1tg000006l | CM024611.1 |  CM024516.1 | 14070625 |  7971233 |  **21.01%** | 13.35% | **97.35%** |  94.97% |
| h1tg000007l | CM024609.1 |  CM024514.1 | 14192176 |  25459023 | 29.80% | **61.75%** | **97.35%** |  96.42% |
| h1tg000008l | CM024616.1 |  CM024521.1 | 19421761 |  32622160 | 29.67% | **50.27%** | **97.35%** |  **97.35%** |
| h1tg000009l | CM024617.1 |  CM024522.1 | 22115718 |  11535835 | **49.93%** | 35.01% | 95.64% |  **97.35%** |
| h1tg000010l | CM024610.1 |  CM024515.1 | 11611686 |  18607522 | 22.00% | **39.55%** | **97.35%** |  96.42% |
| h1tg000011l | CM024614.1 |  CM024519.1 | 18324338 |  33323862 | 30.25% | **59.35%** | **96.42%** |  95.64% |

### Haplotype 2 Assignment
| S1 ID | decip ID | virgin ID | decip Matches | virgin Matches | decip Coverage | virgin Coverage | decip Identity | virgin Identity |
| :---  | :---     | :---      | :---          | :---           | :---           | :---            | :---           | :---            |
| h2tg000002l | CM024615.1 |  CM024520.1 | 13540418 | 23333631 | 22.49% | **40.74%** | **97.35%** |   92.96% |
| h2tg000004l | CM024618.1 |  CM024523.1 | 21985119 | 12297418 | **57.54%** | 33.40% | **97.35%** |   **97.35%** |
| h2tg000307l | CM024613.1 |  CM024518.1 | 25985950 | 14930779 | **42.00%** | 26.17% | 93.48% |   **97.35%** |
| h2tg000011l | CM024612.1 |  CM024517.1 | 12352436 | 19166451 | 28.65% | **59.17%** | **97.35%** |   **97.35%** |
| h2tg000007l | CM024619.1 |  CM024524.1 | 26865746 | 13888371 | **57.26%** | 29.48% | 94.82% |   **95.64%** |
| h2tg000009l | CM024611.1 |  CM024516.1 | 16223298 | 26266346 | 24.23% | **44.00%** | **95.64%** |   91.57% |
| h2tg000010l | CM024609.1 |  CM024514.1 | 24752428 | 13217580 | **51.97%** | 32.06% | 88.59% |   **94.38%** |
| h2tg000003l | CM024616.1 |  CM024521.1 | 34587010 | 19185595 | **52.84%** | 29.57% | **97.89%** |   **97.89%** |
| h2tg000006l | CM024617.1 |  CM024522.1 | 12763271 | 20615297 | 28.82% | **62.57%** | 84.85% |   **93.39%** |
| h2tg000008l | CM024610.1 |  CM024515.1 | 16518915 | 9262525 | **31.30%** | 19.69% | 95.64% |   **96.42%** |
| h2tg000001l | CM024614.1 |  CM024519.1 | 32535393 | 18373938 | **53.72%** | 32.72% | **97.35%** |   94.38% |


Based on these results, the identity looks a bit confusing, and maybe has a bug. Generously, one would say that it doesn't help distinguish much between the two parents. But the coverage tells a different story. We get perfect assignment between the two parents - h1 prefers one parent, h2 always prefers the other. You can read these off. In the following table I rearranged them into order based on the E decipiens genome from NCBI.


| Chr# | decipiens-original | virginea-original | decipiens-parent | virginea-parent |
| :--- | :---             | :---     	     | :---      		| :---            |
| Chr1 | CM024609.1 |  CM024514.1 | h2tg000010l | h1tg000007l |
| Chr2 | CM024610.1 |  CM024515.1 | h2tg000008l | h1tg000010l |
| Chr3 | CM024611.1 |  CM024516.1 | h1tg000006l | h2tg000009l |
| Chr4 | CM024612.1 |  CM024517.1 | h1tg000004l | h2tg000011l |
| Chr5 | CM024613.1 |  CM024518.1 | h2tg000307l | h1tg000003l |
| Chr6 | CM024614.1 |  CM024519.1 | h2tg000001l | h1tg000011l |
| Chr7 | CM024615.1 |  CM024520.1 | h1tg000001l | h2tg000002l |
| Chr8 | CM024616.1 |  CM024521.1 | h2tg000003l | h1tg000008l |
| Chr9 | CM024617.1 |  CM024522.1 | h1tg000009l | h2tg000006l |
| Chr10 | CM024618.1 |  CM024523.1 | h2tg000004l | h1tg000002l |
| Chr11 | CM024619.1 |  CM024524.1 | h2tg000007l | h1tg000005l |


## Dotplots vs. putative parents

Let's further examine the putative assignments above with dotplots.

First we'll get the 11 main chromosomes from the two parental genomes. This script also names the chromosomes in the parental genomes with the chromosome numbers, e.g. `CM024609.1_Chr1`

```bash
# ---------------------------------------------------------
# 1. Process E. decipiens
# ---------------------------------------------------------
# IDs: CM024609.1 (Chr 1) through CM024619.1 (Chr 11)
for i in {0..10}; do
    # Calculate the ID (09, 10, 11...) and the Chr number (1, 2, 3...)
    ACC=$(printf "CM0246%02d.1" $((9 + i)))
    CHR="Chr$((1 + i))"
    
    # Extract and rename on the fly
    samtools faidx parental_spp_genomes/E_decipiens.fa $ACC | \
    sed "s/>$ACC/>${ACC}_${CHR}/" >> parental_spp_genomes/E_decipiens_top11.fa
done

# ---------------------------------------------------------
# 2. Process E. virginea
# ---------------------------------------------------------
# IDs: CM024514.1 (Chr 1) through CM024524.1 (Chr 11)
for i in {0..10}; do
    ACC=$(printf "CM0245%02d.1" $((14 + i)))
    CHR="Chr$((1 + i))"
    
    samtools faidx parental_spp_genomes/E_virginea.fa $ACC | \
    sed "s/>$ACC/>${ACC}_${CHR}/" >> parental_spp_genomes/E_virginea_top11.fa
done
```

Now we can make the dot plots:

```bash
#!/bin/bash

# --- File and Directory Names ---
H1="03_hifiasm_assembly/E_phylacis_hap1_top11.fa"
H2="03_hifiasm_assembly/E_phylacis_hap2_top11.fa"
DEC="parental_spp_genomes/E_decipiens_top11_named.fa"
VIR="parental_spp_genomes/E_virginea_top11_named.fa"
OUT_DIR="04_parental_assignment/all_vs_all_dotplots"

mkdir -p $OUT_DIR

# Define the pairs to compare: "Reference,Query,OutputName"
PAIRS=(
    "$H1,$H2,h1_vs_h2"
    "$DEC,$H1,h1_vs_dec"
    "$VIR,$H1,h1_vs_vir"
    "$DEC,$H2,h2_vs_dec"
    "$VIR,$H2,h2_vs_vir"
    "$DEC,$VIR,dec_vs_vir"
)

          
echo "Done! All 6 comparisons are in $OUT_DIR"
```


## All-vs-All Haplotype & Parental Comparison

This table provides a comprehensive overview of the synteny and divergence between the assembled haplotypes (H1, H2) and the two parental reference genomes (*E. decipiens* and *E. virginea*). 

| | **Haplotype 1 (H1)** | **Haplotype 2 (H2)** | **E. decipiens (Ref)** | **E. virginea (Ref)** |
| :--- | :---: | :---: | :---: | :---: |
| **H1** | — | — | — | — |
| **H2** | [![H2 vs H1](04_parental_assignment/all_vs_all_dotplots/h1_vs_h2_plot.png)](04_parental_assignment/all_vs_all_dotplots/h1_vs_h2_plot.png) | — | — | — |
| **Decipiens** | [![Dec vs H1](04_parental_assignment/all_vs_all_dotplots/h1_vs_dec_plot.png)](04_parental_assignment/all_vs_all_dotplots/h1_vs_dec_plot.png) | [![Dec vs H2](04_parental_assignment/all_vs_all_dotplots/h2_vs_dec_plot.png)](04_parental_assignment/all_vs_all_dotplots/h2_vs_dec_plot.png) | — | — |
| **Virginea** | [![Vir vs H1](04_parental_assignment/all_vs_all_dotplots/h1_vs_vir_plot.png)](04_parental_assignment/all_vs_all_dotplots/h1_vs_vir_plot.png) | [![Vir vs H2](04_parental_assignment/all_vs_all_dotplots/h2_vs_vir_plot.png)](04_parental_assignment/all_vs_all_dotplots/h2_vs_vir_plot.png) | [![Vir vs Dec](04_parental_assignment/all_vs_all_dotplots/dec_vs_vir_plot.png)](04_parental_assignment/all_vs_all_dotplots/dec_vs_vir_plot.png) | — |

---
*Note: Click on any thumbnail to view the full-resolution interactive/static dotplot.*

Visual analysis of the dotplots suggest the following.

First, it looks to me like there could be a bit of phase switching. E.g. chr3 for h1 maps well to the start of virginea, but to the middle of decipiens. Merqury should help to sort out if this is real.

Second, we can try to allocate chromosomes to parents. However, doing this reveals that there's often not a clear distinction, and even if there is, one would often assign both H1 and H2 to the same parents on visual inspection (e.g. Chr 5 looks to be best aligned to decipiens for both H1 and H2). This could just be because the particular accumulation of structural variants (which are common) in the particular samples we have access to for virginea and decipiens. 

Everyone here is clearly closely related, and there are clearly a lot of structural variants between all comparisons. 

Another notable difference is that mapping E. decipiens vs. E. virginea looks a lot more contiguous than H1 vs H2, or either haplotypes vs decipiens OR virginea. This could well be because the two parental genomes were scaffolded against E. grandis, while H1 and H2 were not. This is important, because this is likely to reduce apparent variation between decipiens and virginea, as both will be affected by reference bias. Of course, my assembly may be full of errors as well!

All of this suggests that Merqury might be a useful approach to paint chromosomes. 

## Merqury

First let's download the short reads for the two parent species. They are here: 

* decipiens: https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR25119757
* virginea: https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR25119786

```bash
cd parental_spp_genomes

# decipiens
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR025/25119/SRR25119757/SRR25119757.lite.1
fasterq-dump SRR25119757.lite.1 --split-files --threads 16 --progress && pigz -p 16 *.fastq
mv SRR25119757.lite.1_1.fastq.gz E_decipiens_R1.fastq.gz
mv SRR25119757.lite.1_2.fastq.gz E_decipiens_R2.fastq.gz

# virginea
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR025/25119/SRR25119786/SRR25119786.lite.1
fasterq-dump SRR25119786.lite.1 --split-files --threads 16 --progress && pigz -p 16 *.fastq
mv SRR25119786.lite.1_1.fastq.gz E_virginea_R1.fastq.gz
mv SRR25119786.lite.1_2.fastq.gz E_virginea_R2.fastq.gz
```

Now we'll build the kmer databases. note that these are not actual parents, but parental species, so I won't bother using the kmers to check the assembly itself. I'm only intersted in mapping parental species kmers to the assembly to look for switching.

```bash
export k=31
cd parental_spp_genomes/

# E virginea
meryl k=$k count threads=128 memory=1200G \
    output E_virginea.meryl \
    E_virginea_R1.fastq.gz E_virginea_R2.fastq.gz

# E decipiens
meryl k=$k count threads=128 memory=1200G \
    output E_decipiens.meryl \
    E_decipiens_R1.fastq.gz E_decipiens_R2.fastq.gz

meryl statistics E_virginea.meryl
meryl statistics E_decipiens.meryl

```

```
decipiens
Found 1 command tree.
Number of 31-mers that are:
  unique             1290874493  (exactly one instance of the kmer is in the input)
  distinct           2030909683  (non-redundant kmer sequences in the input)
  present           18810934773  (...)
  missing   4611686016396478221  (non-redundant kmer sequences not in the input)

             number of   cumulative   cumulative     presence
              distinct     fraction     fraction   in dataset
frequency        kmers     distinct        total       (1e-6)
--------- ------------ ------------ ------------ ------------
        1   1290874493       0.6356       0.0686     0.000053
        2     89038578       0.6795       0.0781     0.000106
        3     34417617       0.6964       0.0836     0.000159
        4     24398445       0.7084       0.0888     0.000213
        5     21071260       0.7188       0.0944     0.000266
        6     20750383       0.7290       0.1010     0.000319
        7     21256492       0.7395       0.1089     0.000372
        8     21951538       0.7503       0.1182     0.000425
        9     22547668       0.7614       0.1290     0.000478
       10     22978776       0.7727       0.1412     0.000532

virginea
Found 1 command tree.
Number of 31-mers that are:
  unique             1225112802  (exactly one instance of the kmer is in the input)
  distinct           1897334044  (non-redundant kmer sequences in the input)
  present           17081200798  (...)
  missing   4611686016530053860  (non-redundant kmer sequences not in the input)

             number of   cumulative   cumulative     presence
              distinct     fraction     fraction   in dataset
frequency        kmers     distinct        total       (1e-6)
--------- ------------ ------------ ------------ ------------
        1   1225112802       0.6457       0.0717     0.000059
        2     94219756       0.6954       0.0828     0.000117
        3     32918583       0.7127       0.0885     0.000176
        4     21986254       0.7243       0.0937     0.000234
        5     18523886       0.7341       0.0991     0.000293
        6     18102892       0.7436       0.1055     0.000351
        7     18383380       0.7533       0.1130     0.000410
        8     18817127       0.7632       0.1218     0.000468
        9     19189702       0.7733       0.1319     0.000527
       10     19389542       0.7835       0.1433     0.000585
```

The trough in both of these distributions is at 6. So we'll get unique 31mers by filtering out everything that appears <6 times.

```bash
meryl threads=128 greater-than 5 \
    output E_virginea.filtered.meryl \
    E_virginea.meryl

meryl threads=128 greater-than 5 \
    output E_decipiens.filtered.meryl \
    E_decipiens.meryl
```

Now we can create species-spefic 31mers via subtraction

```bash
meryl threads=128 subtract \
    output E_virginea.only.meryl \
    E_virginea.filtered.meryl E_decipiens.filtered.meryl

meryl threads=128 subtract \
    output E_decipiens.only.meryl \
    E_decipiens.filtered.meryl E_virginea.filtered.meryl

meryl statistics E_virginea.only.meryl | head -n 20
meryl statistics E_decipiens.only.meryl | head -n 20
```

Both of these have a similar distribution:

```
decipiens
Found 1 command tree.
Number of 31-mers that are:
  unique                4979883  (exactly one instance of the kmer is in the input)
  distinct            486744832  (non-redundant kmer sequences in the input)
  present           10560074052  (...)
  missing   4611686017940643072  (non-redundant kmer sequences not in the input)

             number of   cumulative   cumulative     presence
              distinct     fraction     fraction   in dataset
frequency        kmers     distinct        total       (1e-6)
--------- ------------ ------------ ------------ ------------
        1      4979883       0.0102       0.0005     0.000095
        2      4832476       0.0202       0.0014     0.000189
        3      4666418       0.0297       0.0027     0.000284
        4      4474050       0.0389       0.0044     0.000379
        5      4263121       0.0477       0.0064     0.000473
        6     22101931       0.0931       0.0190     0.000568
        7     22034781       0.1384       0.0336     0.000663
        8     22173871       0.1839       0.0504     0.000758
        9     22250304       0.2296       0.0694     0.000852
       10     22199813       0.2753       0.0904     0.000947

virginea
Found 1 command tree.
Number of 31-mers that are:
  unique                4954189  (exactly one instance of the kmer is in the input)
  distinct            414255629  (non-redundant kmer sequences in the input)
  present            8912622148  (...)
  missing   4611686018013132275  (non-redundant kmer sequences not in the input)

             number of   cumulative   cumulative     presence
              distinct     fraction     fraction   in dataset
frequency        kmers     distinct        total       (1e-6)
--------- ------------ ------------ ------------ ------------
        1      4954189       0.0120       0.0006     0.000112
        2      4770616       0.0235       0.0016     0.000224
        3      4571300       0.0345       0.0032     0.000337
        4      4346218       0.0450       0.0051     0.000449
        5      4116279       0.0549       0.0074     0.000561
        6     19090626       0.1010       0.0203     0.000673
        7     18795727       0.1464       0.0350     0.000785
        8     18693943       0.1915       0.0518     0.000898
        9     18568529       0.2363       0.0706     0.001010
       10     18318337       0.2806       0.0911     0.001122

```

Now we filter out the rare unique kmers again, and use merqury to map these to the two haplotypes.

```bash

meryl threads=128 greater-than 5 \
    output E_virginea.final_probes.meryl \
    E_virginea.only.meryl

meryl threads=128 greater-than 5 \
    output E_decipiens.final_probes.meryl \
    E_decipiens.only.meryl
```

Now we run Merqury. I generate a dummy child db (I don't have short reads for the child, and regardless it's not actually a child...)

```bash
# 1. Create a placeholder 'read-db' by merging the two probe sets
meryl threads=128 union-sum \
    output species_union.meryl \
    E_virginea.final_probes.meryl E_decipiens.final_probes.meryl

# copy the files with a .fasta extension, because merqury seems to reall want that
cp ../03_hifiasm_assembly/E_phylacis_hap1_top11.fa E_phylacis_hap1_top11.fasta
cp ../03_hifiasm_assembly/E_phylacis_hap2_top11.fa E_phylacis_hap2_top11.fasta

merqury.sh \
    species_union.meryl \
    E_virginea.final_probes.meryl \
    E_decipiens.final_probes.meryl \
    E_phylacis_hap1_top11.fasta \
    E_phylacis_hap2_top11.fasta \
    E_phylacis_species_map
```


Frustratingly, this fails at the last steps, likely because using parental species not individuals leads to a lot of noise and ~infinite switching, so no blocks long enough to call a real block by merqury's approach. Instead, let's map the kmers directly and plot them out. 

Of note, the blob plot and the histogram (same data) show a consistent 2-3x preference for one parent spp over the other:

| Assembly | Contig | E_virginea Probes | E_decipiens Probes | Size (bp) | Dominant Parent | Dominant Ratio |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| E_phylacis_hap1_top11 | h1tg000008l | 29,574,230 | 11,485,726 | 61,728,875 | Virginea | 2.57x |
| E_phylacis_hap1_top11 | h1tg000003l | 24,166,149 | 12,032,714 | 59,377,806 | Virginea | 2.01x |
| E_phylacis_hap1_top11 | h1tg000001l | 7,233,487 | 25,796,697 | 51,317,245 | Decipiens | 3.57x |
| E_phylacis_hap1_top11 | h1tg000011l | 26,173,487 | 9,579,941 | 51,144,841 | Virginea | 2.73x |
| E_phylacis_hap1_top11 | h1tg000007l | 21,214,760 | 7,752,565 | 41,957,492 | Virginea | 2.74x |
| E_phylacis_hap1_top11 | h1tg000009l | 5,729,451 | 20,405,351 | 39,833,786 | Decipiens | 3.56x |
| E_phylacis_hap1_top11 | h1tg000005l | 20,178,769 | 7,482,591 | 39,657,493 | Virginea | 2.70x |
| E_phylacis_hap1_top11 | h1tg000004l | 5,222,786 | 19,422,560 | 38,027,346 | Decipiens | 3.72x |
| E_phylacis_hap1_top11 | h1tg000002l | 17,617,635 | 6,408,337 | 34,632,705 | Virginea | 2.75x |
| E_phylacis_hap1_top11 | h1tg000010l | 15,674,330 | 6,474,976 | 33,387,254 | Virginea | 2.42x |
| E_phylacis_hap1_top11 | h1tg000006l | 4,178,238 | 13,450,401 | 28,751,513 | Decipiens | 3.22x |
| E_phylacis_hap2_top11 | h2tg000003l | 10,174,689 | 33,286,426 | 65,679,558 | Decipiens | 3.27x |
| E_phylacis_hap2_top11 | h2tg000307l | 8,408,538 | 28,955,231 | 59,218,057 | Decipiens | 3.44x |
| E_phylacis_hap2_top11 | h2tg000009l | 23,820,241 | 11,758,129 | 56,292,507 | Virginea | 2.03x |
| E_phylacis_hap2_top11 | h2tg000001l | 7,696,497 | 29,401,315 | 54,055,966 | Decipiens | 3.82x |
| E_phylacis_hap2_top11 | h2tg000002l | 18,573,496 | 10,028,448 | 45,187,031 | Virginea | 1.85x |
| E_phylacis_hap2_top11 | h2tg000010l | 6,147,171 | 23,104,469 | 43,786,441 | Decipiens | 3.76x |
| E_phylacis_hap2_top11 | h2tg000007l | 5,967,530 | 23,348,547 | 42,203,106 | Decipiens | 3.91x |
| E_phylacis_hap2_top11 | h2tg000004l | 5,520,519 | 19,783,505 | 38,512,705 | Decipiens | 3.58x |
| E_phylacis_hap2_top11 | h2tg000006l | 17,538,735 | 7,250,381 | 36,524,821 | Virginea | 2.42x |
| E_phylacis_hap2_top11 | h2tg000011l | 17,323,908 | 6,663,693 | 35,797,366 | Virginea | 2.60x |
| E_phylacis_hap2_top11 | h2tg000008l | 5,139,230 | 14,287,990 | 30,255,362 | Decipiens | 2.78x |

[link to blob plot]

There are two explanations for this:

1. Frequent phase switching - seems unlikely because the ratio is so consistent, and there's always a clear dominant parent
2. Background noise - these are three individuals that are very different from one another, so it's feasible that the final probes (e.g. those uniquely in E. virginea that are also in the Meelup mallee) have a fair degree of noise, insofar as some proportion of them are also likely to hit E. decipiens chromosomes.

One way to distinguish these is to map the two final probe sets directly to the chromosomes, and look at the distribution. If 1 is correct, then we expect clear blocks that prefer one parent over the other (even with noise). If 2 is correct, then we expect a consistent signal along a chromosome where one parent is preferred compared to the other in a sliding window. The evidence from merqury failing suggests the window could be quite small, since it looks like there weren't runs of >100 unique kmers for any of the parents. So we could look at windows of ~1000 and if 2 is correct it should be a stable signal of ~2-3x preference for the dominant parent. 

### Mapping kmers to chromosomes

Let's map the kmers to the chromosomes with meryl, and do a sliding window to visualise it. I'll aim for two visualisations

First let's map the kmers with meryl:

```bash
# For Hap1
# Virginea k-mers on Hap1
meryl-lookup -bed \
  -sequence E_phylacis_hap1_top11.fasta \
  -mers E_virginea.final_probes.meryl \
  -output hap1_virginea.bed

# Decipiens k-mers on Hap1
meryl-lookup -bed \
  -sequence E_phylacis_hap1_top11.fasta \
  -mers E_decipiens.final_probes.meryl \
  -output hap1_decipiens.bed

# For Hap2
# Virginea k-mers on Hap2
meryl-lookup -bed \
  -sequence E_phylacis_hap2_top11.fasta \
  -mers E_virginea.final_probes.meryl \
  -output hap2_virginea.bed

# Decipiens k-mers on Hap2
meryl-lookup -bed \
  -sequence E_phylacis_hap2_top11.fasta \
  -mers E_decipiens.final_probes.meryl \
  -output hap2_decipiens.bed

```

Define boundaries:

```bash
# Create chromosome index files
samtools faidx E_phylacis_hap1_top11.fasta
samtools faidx E_phylacis_hap2_top11.fasta

# Extract contig names and lengths
cut -f 1,2 E_phylacis_hap1_top11.fasta.fai > hap1.genome
cut -f 1,2 E_phylacis_hap2_top11.fasta.fai > hap2.genome
```

make windows, I'll do 10K, and 100K to start with.

```bash
# Windows for Hap1
bedtools makewindows -g hap1.genome -w 10000 > hap1_10k.bed
bedtools makewindows -g hap1.genome -w 100000 > hap1_100k.bed

# Windows for Hap2
bedtools makewindows -g hap2.genome -w 10000 > hap2_10k.bed
bedtools makewindows -g hap2.genome -w 100000 > hap2_100k.bed
```

Now I'll get kmer coverage:

```bash
# Hap1 #
# 100k Windows
bedtools coverage -a hap1_100k.bed -b hap1_virginea.bed > h1_v_100k.txt
bedtools coverage -a hap1_100k.bed -b hap1_decipiens.bed > h1_d_100k.txt

# 10k Windows
bedtools coverage -a hap1_10k.bed -b hap1_virginea.bed > h1_v_10k.txt
bedtools coverage -a hap1_10k.bed -b hap1_decipiens.bed > h1_d_10k.txt

# Hap2 #
# 100k Windows
bedtools coverage -a hap2_100k.bed -b hap2_virginea.bed > h2_v_100k.txt
bedtools coverage -a hap2_100k.bed -b hap2_decipiens.bed > h2_d_100k.txt

# 10k Windows
bedtools coverage -a hap2_10k.bed -b hap2_virginea.bed > h2_v_10k.txt
bedtools coverage -a hap2_10k.bed -b hap2_decipiens.bed > h2_d_10k.txt
```

Now we plot those out along the genomes with an R script:

```bash
Rscript ../scripts/parental_kmer_plots.R . parental_kmer_plots
```

The plots below validate the assembly phasing and investigate parental inheritance patterns by mapping species-specific k-mer probes (derived from *E. virginea* and *E. decipiens* SRA data) across both haplotypes. 

The plots below compare raw k-mer density (counts per window) and the relative proportion of parental signal. Note that reciprocal crossovers in the proportion plots would indicate potential phase switches, while asymmetric "flickers" suggest complex ancestry or introgression. The former are very rare, suggesting that the phasing is good.

| Raw K-mer Density | Parental Proportions (Mirror Plots) |
| :---: | :---: |
| **100K Windows** | **100K Windows** |
| <img src="./04_parental_assignment/parental_kmer_plots/raw_density_100k.png" width="400"> | <img src="./04_parental_assignment/parental_kmer_plots/proportion_100k.png" width="400"> |
| **10K Windows** | **10K Windows** |
| <img src="./04_parental_assignment/parental_kmer_plots/raw_density_10k.png" width="400"> | <img src="./04_parental_assignment/parental_kmer_plots/proportion_10k.png" width="400"> |


Next, it's time to tidy everything up. 

```bash
cd ../04_parental_assignment
mkdir species_map
mv ../parental_spp_genomes/E_phylacis_* species_map/
mv ../parental_spp_genomes/completeness.stats species_map/
mv ../parental_spp_genomes/h* species_map/
mv ../parental_spp_genomes/parental_kmer_plots .
```

I added the species map file to the gitignore


## Parental assignment of the 11 major chromosomes

To summarise all of the above, we can now clearly assign the haplotypes into two parental groups.

This table summarizes the evidence used to assign the `h1` and `h2` contigs to their parental identities. Ratios in brackets indicate the strength of the signal (Dominant Parent / Secondary Parent).

What's amazing (and comforting) is that genome-alignment and unique parental species kmers agree perfectly. 

| Chr# | decipiens-original | virginea-original | h1 ID | h2 ID | h1 decip cov | h1 virg cov | h2 decip cov | h2 virg cov | h1 decip probes | h1 virg probes | h2 decip probes | h2 virg probes | h1 coverage parent | h1 probe parent | h2 coverage parent | h2 probe parent |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Chr1** | CM024609.1 | CM024514.1 | h1tg000007l | h2tg000010l | 29.8% | 61.8% | 52.0% | 32.1% | 7.7M | 21.2M | 23.1M | 6.1M | virginea (2.07x) | virginea (2.74x) | decipiens (1.62x) | decipiens (3.76x) |
| **Chr2** | CM024610.1 | CM024515.1 | h1tg000010l | h2tg000008l | 22.0% | 39.6% | 31.3% | 19.7% | 6.5M | 15.7M | 14.3M | 5.1M | virginea (1.80x) | virginea (2.42x) | decipiens (1.59x) | decipiens (2.78x) |
| **Chr3** | CM024611.1 | CM024516.1 | h1tg000006l | h2tg000009l | 21.0% | 13.4% | 24.2% | 44.0% | 13.5M | 4.2M | 11.8M | 23.8M | decipiens (1.57x) | decipiens (3.22x) | virginea (1.82x) | virginea (2.03x) |
| **Chr4** | CM024612.1 | CM024517.1 | h1tg000004l | h2tg000011l | 49.1% | 31.6% | 28.7% | 59.2% | 19.4M | 5.2M | 6.7M | 17.3M | decipiens (1.55x) | decipiens (3.72x) | virginea (2.07x) | virginea (2.60x) |
| **Chr5** | CM024613.1 | CM024518.1 | h1tg000003l | h2tg000307l | 24.6% | 47.0% | 42.0% | 26.2% | 12.0M | 24.2M | 29.0M | 8.4M | virginea (1.91x) | virginea (2.01x) | decipiens (1.60x) | decipiens (3.44x) |
| **Chr6** | CM024614.1 | CM024519.1 | h1tg000011l | h2tg000001l | 30.3% | 59.4% | 53.7% | 32.7% | 9.6M | 26.2M | 29.4M | 7.7M | virginea (1.96x) | virginea (2.73x) | decipiens (1.64x) | decipiens (3.82x) |
| **Chr7** | CM024615.1 | CM024520.1 | h1tg000001l | h2tg000002l | 40.7% | 24.9% | 22.5% | 40.7% | 25.8M | 7.2M | 10.0M | 18.6M | decipiens (1.63x) | decipiens (3.57x) | virginea (1.81x) | virginea (1.85x) |
| **Chr8** | CM024616.1 | CM024521.1 | h1tg000008l | h2tg000003l | 29.7% | 50.3% | 52.8% | 29.6% | 11.5M | 29.6M | 33.3M | 10.2M | virginea (1.69x) | virginea (2.57x) | decipiens (1.79x) | decipiens (3.27x) |
| **Chr9** | CM024617.1 | CM024522.1 | h1tg000009l | h2tg000006l | 49.9% | 35.0% | 28.8% | 62.6% | 20.4M | 5.7M | 7.3M | 17.5M | decipiens (1.43x) | decipiens (3.56x) | virginea (2.17x) | virginea (2.42x) |
| **Chr10** | CM024618.1 | CM024523.1 | h1tg000002l | h2tg000004l | 32.7% | 59.6% | 57.5% | 33.4% | 6.4M | 17.6M | 19.8M | 5.5M | virginea (1.82x) | virginea (2.75x) | decipiens (1.72x) | decipiens (3.58x) |
| **Chr11** | CM024619.1 | CM024524.1 | h1tg000005l | h2tg000007l | 31.4% | 54.9% | 57.3% | 29.5% | 7.5M | 20.2M | 23.3M | 6.0M | virginea (1.75x) | virginea (2.70x) | decipiens (1.94x) | decipiens (3.91x) |

Let's rearrange the 11 chromosomes of the two haplotypes accordingly:

```bash
# Define paths
ASSEMBLY_DIR="03_hifiasm_assembly"
HAP1="$ASSEMBLY_DIR/E_phylacis_hap1_top11.fa"
HAP2="$ASSEMBLY_DIR/E_phylacis_hap2_top11.fa"

# Define output files
VIRG_OUT="$ASSEMBLY_DIR/E_phylacis_virginea_parent_top11.fa"
DECIP_OUT="$ASSEMBLY_DIR/E_phylacis_decipiens_parent_top11.fa"

# Helper function to extract and rename
# Usage: rename_chr <input_file> <contig_id> <chr_number>
rename_chr() {
    seqkit grep -p "$2" "$1" | \
    seqkit replace -p "(.+)" -r "\$1 Eucalyptus phylacis chromosome $3"
}

echo "Sorting and renaming chromosomes..."

# --- 1. Assemble Virginea Parent (E. virginea) ---
# Chromosomes: 1(h1), 2(h1), 3(h2), 4(h2), 5(h1), 6(h1), 7(h2), 8(h1), 9(h2), 10(h1), 11(h1)
{
    rename_chr "$HAP1" "h1tg000007l" "1"
    rename_chr "$HAP1" "h1tg000010l" "2"
    rename_chr "$HAP2" "h2tg000009l" "3"
    rename_chr "$HAP2" "h2tg000011l" "4"
    rename_chr "$HAP1" "h1tg000003l" "5"
    rename_chr "$HAP1" "h1tg000011l" "6"
    rename_chr "$HAP2" "h2tg000002l" "7"
    rename_chr "$HAP1" "h1tg000008l" "8"
    rename_chr "$HAP2" "h2tg000006l" "9"
    rename_chr "$HAP1" "h1tg000002l" "10"
    rename_chr "$HAP1" "h1tg000005l" "11"
} > "$VIRG_OUT"

# --- 2. Assemble Decipiens Parent (E. decipiens) ---
# Chromosomes: 1(h2), 2(h2), 3(h1), 4(h1), 5(h2), 6(h2), 7(h1), 8(h2), 9(h1), 10(h2), 11(h2)
{
    rename_chr "$HAP2" "h2tg000010l" "1"
    rename_chr "$HAP2" "h2tg000008l" "2"
    rename_chr "$HAP1" "h1tg000006l" "3"
    rename_chr "$HAP1" "h1tg000004l" "4"
    rename_chr "$HAP2" "h2tg000307l" "5"
    rename_chr "$HAP2" "h2tg000001l" "6"
    rename_chr "$HAP1" "h1tg000001l" "7"
    rename_chr "$HAP2" "h2tg000003l" "8"
    rename_chr "$HAP1" "h1tg000009l" "9"
    rename_chr "$HAP2" "h2tg000004l" "10"
    rename_chr "$HAP2" "h2tg000007l" "11"
} > "$DECIP_OUT"

echo "Done."
echo "------------------------------------------"
echo "VIRGINEA FILE: $VIRG_OUT"
grep ">" "$VIRG_OUT"
echo "------------------------------------------"
echo "DECIPIENS FILE: $DECIP_OUT"
grep ">" "$DECIP_OUT"

mv 03_hifiasm_assembly/E_phylacis_virginea_parent_top11.fa 04_parental_assignment/
mv 03_hifiasm_assembly/E_phylacis_decipiens_parent_top11.fa 04_parental_assignment/
```

## Who's the mother?

Now let's use the cp genome to figure out who the mother is, then add that chromosome to the assembly.

First let's get the two cp genomes of the parents:

```bash
seqkit grep -p CM024632.1 parental_spp_genomes/E_decipiens.fa > 04_parental_assignment/E_decipiens_cp_ref.fa
seqkit grep -p CM024525.1 parental_spp_genomes/E_virginea.fa > 04_parental_assignment/E_virginea_cp_ref.fa
```

```bash
cd 04_parental_assignment

# Create blast database
makeblastdb -in 03_hifiasm_assembly/E_phylacis_asm.bp.p_ctg.fa -dbtype nucl -out e_phylacis_db

# Run BLAST
blastn -query E_virginea_cp_ref.fa -db e_phylacis_db \
       -outfmt "7 sseqid pident length qlen slen qcovs evalue" \
       -out cp_blast_results_virginea.txt

blastn -query E_decipiens_cp_ref.fa -db e_phylacis_db \
       -outfmt "7 sseqid pident length qlen slen qcovs evalue" \
       -out cp_blast_results_decipiens.txt

# Look at the top hit
head -n 25 cp_blast_results_decipiens.txt
head -n 25 cp_blast_results_virginea.txt
```

This clearly pulls out `ptg000022l` as the chloroplast contig:

```
(eucalypt_asm) rob@nova:~/MM_genome/04_parental_assignment$ head -n 35 cp_blast_results_decipiens.txt
# BLASTN 2.16.0+
# Query: CM024632.1 Eucalyptus decipiens isolate CCA4570 chloroplast, complete sequence, whole genome shotgun sequence
# Database: e_phylacis_db
# Fields: subject id, % identity, alignment length, query length, subject length, % query coverage per subject, evalue
# 4924 hits found
ptg000022l      99.357  49925   163268  198530  99      0.0
ptg000022l      99.674  35302   163268  198530  99      0.0
ptg000022l      99.095  34265   163268  198530  99      0.0
ptg000022l      99.281  31705   163268  198530  99      0.0
ptg000022l      99.641  29785   163268  198530  99      0.0
ptg000022l      99.788  26402   163268  198530  99      0.0
ptg000022l      99.739  26429   163268  198530  99      0.0
ptg000022l      99.790  26138   163268  198530  99      0.0
ptg000022l      99.346  7037    163268  198530  99      0.0
ptg000022l      99.033  6825    163268  198530  99      0.0
ptg000022l      98.860  4123    163268  198530  99      0.0
ptg000022l      99.384  3086    163268  198530  99      0.0
ptg000022l      98.182  55      163268  198530  99      5.08e-16
ptg000022l      100.000 39      163268  198530  99      8.56e-09
ptg000022l      97.619  42      163268  198530  99      8.56e-09
ptg000022l      97.500  40      163268  198530  99      1.11e-07
ptg000022l      97.500  40      163268  198530  99      1.11e-07
ptg000022l      97.500  40      163268  198530  99      1.11e-07
ptg000022l      97.500  40      163268  198530  99      1.11e-07
ptg000022l      97.500  40      163268  198530  99      1.11e-07
ptg000022l      79.121  91      163268  198530  99      2.40e-04
ptg000022l      100.000 30      163268  198530  99      8.62e-04
ptg000022l      100.000 30      163268  198530  99      8.62e-04
ptg000022l      100.000 30      163268  198530  99      8.62e-04
ptg000022l      100.000 28      163268  198530  99      0.011
ptg000033l      99.291  38084   163268  79721   53      0.0
ptg000033l      99.674  35566   163268  79721   53      0.0
ptg000033l      99.638  27101   163268  79721   53      0.0
ptg000033l      99.788  26402   163268  79721   53      0.0
ptg000033l      98.836  4124    163268  79721   53      0.0
(eucalypt_asm) rob@nova:~/MM_genome/04_parental_assignment$ head -n 35 cp_blast_results_virginea.txt
# BLASTN 2.16.0+
# Query: CM024525.1 Eucalyptus virginea isolate CCA4170 chloroplast, complete sequence, whole genome shotgun sequence
# Database: e_phylacis_db
# Fields: subject id, % identity, alignment length, query length, subject length, % query coverage per subject, evalue
# 4829 hits found
ptg000022l      99.794  80388   160251  198530  99      0.0
ptg000022l      99.699  79985   160251  198530  99      0.0
ptg000022l      99.892  35265   160251  198530  99      0.0
ptg000022l      99.883  26413   160251  198530  99      0.0
ptg000022l      99.788  26415   160251  198530  99      0.0
ptg000022l      99.786  26151   160251  198530  99      0.0
ptg000022l      99.935  3086    160251  198530  99      0.0
ptg000022l      100.000 39      160251  198530  99      8.40e-09
ptg000022l      97.619  42      160251  198530  99      8.40e-09
ptg000022l      97.500  40      160251  198530  99      1.09e-07
ptg000022l      97.500  40      160251  198530  99      1.09e-07
ptg000022l      97.500  40      160251  198530  99      1.09e-07
ptg000022l      97.500  40      160251  198530  99      1.09e-07
ptg000022l      97.500  40      160251  198530  99      1.09e-07
ptg000022l      79.121  91      160251  198530  99      2.35e-04
ptg000022l      100.000 30      160251  198530  99      8.46e-04
ptg000022l      100.000 30      160251  198530  99      8.46e-04
ptg000022l      100.000 30      160251  198530  99      8.46e-04
ptg000022l      100.000 28      160251  198530  99      0.011
ptg000021l      99.699  79985   160251  272542  77      0.0
ptg000021l      99.902  15309   160251  272542  77      0.0
ptg000021l      99.726  15312   160251  272542  77      0.0
ptg000021l      99.763  10966   160251  272542  77      0.0
ptg000021l      99.838  1851    160251  272542  77      0.0
ptg000021l      99.838  1849    160251  272542  77      0.0
ptg000021l      99.028  1851    160251  272542  77      0.0
ptg000021l      74.016  889     160251  272542  77      1.63e-80
ptg000021l      73.903  889     160251  272542  77      7.61e-79
ptg000021l      74.943  435     160251  272542  77      6.18e-45
ptg000021l      100.000 39      160251  272542  77      8.40e-09
```

Now let's align the three cp genomes and see what we see:

```bash
seqkit grep -p ptg000022l ../03_hifiasm_assembly/E_phylacis_asm.bp.p_ctg.fa > E_phylacis_cp_ref.fa

minimap2 -x asm5 E_decipiens_cp_ref.fa E_phylacis_cp_ref.fa > cp_phylacis_vs_decipiens.paf
minimap2 -x asm5 E_virginea_cp_ref.fa E_phylacis_cp_ref.fa > cp_phylacis_vs_virginea.paf
```

Comparing cols 10 (matching bases) and 11 (alignment length) makes it look a lot like virginea is the winner. It has 113749/115277 = 98.7% matching, while decipiens has 110447/115514 = 95.6% matching: 

```
(eucalypt_asm) rob@nova:~/MM_genome/04_parental_assignment$ more cp_phylacis_vs_decipiens.paf 
ptg000022l      198530  83434   198526  +       CM024632.1      163268  18      114838  110447  115514  60      tp:A:P  cm:i:10934      s1:i:110216     s2:i:25721      dv:f:0.0004     rl:i:0
ptg000022l      198530  12      83388   +       CM024632.1      163268  76475   163265  80252   86911   60      tp:A:P  cm:i:7928       s1:i:79660      s2:i:54296      dv:f:0.0010     rl:i:0
ptg000022l      198530  172420  198526  -       CM024632.1      163268  137154  163265  25733   26131   0       tp:A:S  cm:i:2560       s1:i:25721      dv:f:0.0009     rl:i:0

(eucalypt_asm) rob@nova:~/MM_genome/04_parental_assignment$ more cp_phylacis_vs_virginea.paf 
ptg000022l      198530  83434   198526  -       CM024525.1      160251  45056   160234  113749  115277  60      tp:A:P  cm:i:11365      s1:i:113681     s2:i:25669      dv:f:0.0002     rl:i:0
ptg000022l      198530  12      83415   -       CM024525.1      160251  2       83403   82339   83489   60      tp:A:P  cm:i:8222       s1:i:82305      s2:i:52099      dv:f:0.0007     rl:i:0
ptg000022l      198530  172393  198526  +       CM024525.1      160251  2       26106   25677   26144   0       tp:A:S  cm:i:2552       s1:i:25669      dv:f:0.0015     rl:i:0

```

Let's also do dotplots to visualise it:

```bash
# Dot plot vs. Virginea
Rscript ../scripts/pafCoordsDotPlotly.R \
    -i cp_phylacis_vs_virginea.paf \
    -o cp_phylacis_vs_virginea \
    -s -t -m 500 -q 100000

# Dot plot vs. Decipiens
Rscript ../scripts/pafCoordsDotPlotly.R \
    -i cp_phylacis_vs_decipiens.paf \
    -o cp_phylacis_vs_decipiens \
    -s -t -m 500 -q 100000
```

| Comparison: *E. phylacis* vs. *E. virginea* | Comparison: *E. phylacis* vs. *E. decipiens* |
| :---: | :---: |
| ![Phylacis vs Virginea](04_parental_assignment/cp_phylacis_vs_virginea.png) | ![Phylacis vs Decipiens](04_parental_assignment/cp_phylacis_vs_decipiens.png) |


The dotplots reveal nothing amiss. 

The mother is clearly E virginea.

# 05 Organellar genome assembly

## Read filtering

```bash
# Setup Directories
source_file="02_filtering/E_phylacis_filtered.fastq.gz"
temp_fq="05_organelle_genomes/temp_unzipped.fastq"
organelle_reads_dir="05_organelle_genomes/organelle_subset"
mkdir -p ${organelle_reads_dir}

# 1. Decompress once to save CPU/IO
echo "Decompressing source file..."
pigz -dc -p 64 $source_file > $temp_fq

# 2. Run the sweep
for Q in 15 20; do
    echo "--- Filtering for Q${Q} ---"
    
    # Process from the temp unzipped file
    chopper -q $Q < $temp_fq | bgzip -@ 64 > ${organelle_reads_dir}/E_phylacis_Q${Q}.fastq.gz
    
    # Report results
    echo "Stats for Q${Q}:"
    seqkit stats ${organelle_reads_dir}/E_phylacis_Q${Q}.fastq.gz
done

# shows that I have a lot of Q20 or greater data (see below) so now let's try read lengths:

# Sweep for Length while keeping Quality constant at Q20
for L in 20000 30000 40000 50000; do
    # Convert length to 'k' format for filename
    L_kb=$(($L / 1000))
    echo "--- Filtering for Q20 and Length > ${L_kb}kb ---"
    
    # Process: filter for Q20 AND length L
    chopper -q 20 -l $L < $temp_fq | bgzip -@ 64 > ${organelle_reads_dir}/E_phylacis_Q20_L${L_kb}k.fastq.gz
    
    # Report results
    echo "Stats for Q20_L${L_kb}k:"
    seqkit stats ${organelle_reads_dir}/E_phylacis_Q20_L${L_kb}k.fastq.gz
done

# 3. Cleanup
rm $temp_fq
echo "All subsets created in ${organelle_reads_dir} and temp file removed."

```



This gives:

```
--- Filtering for Q15 ---
Kept 1165231 reads out of 1345923 reads
Stats for Q15:
file                                                           format  type   num_seqs         sum_len  min_len   avg_len  max_len
05_organelle_genomes/organelle_subset/E_phylacis_Q15.fastq.gz  FASTQ   DNA   1,165,231  29,326,043,409   15,000  25,167.6  147,172
--- Filtering for Q20 ---
Kept 696853 reads out of 1345923 reads
Stats for Q20:
file                                                           format  type  num_seqs         sum_len  min_len   avg_len  max_len
05_organelle_genomes/organelle_subset/E_phylacis_Q20.fastq.gz  FASTQ   DNA    696,853  17,203,549,475   15,000  24,687.5  134,734


--- Filtering for Q20 and Length > 20kb ---
Kept 463140 reads out of 1345923 reads
Stats for Q20_L20k:
file                                                                format  type  num_seqs         sum_len  min_len   avg_len  max_len
05_organelle_genomes/organelle_subset/E_phylacis_Q20_L20k.fastq.gz  FASTQ   DNA    463,140  13,056,995,216   20,000  28,192.3  134,734
--- Filtering for Q20 and Length > 30kb ---
Kept 144088 reads out of 1345923 reads
Stats for Q20_L30k:
file                                                                format  type  num_seqs        sum_len  min_len   avg_len  max_len
05_organelle_genomes/organelle_subset/E_phylacis_Q20_L30k.fastq.gz  FASTQ   DNA    144,088  5,333,865,113   30,000  37,018.1  134,734
--- Filtering for Q20 and Length > 40kb ---
Kept 34958 reads out of 1345923 reads
Stats for Q20_L40k:
file                                                                format  type  num_seqs        sum_len  min_len   avg_len  max_len
05_organelle_genomes/organelle_subset/E_phylacis_Q20_L40k.fastq.gz  FASTQ   DNA     34,958  1,627,555,827   40,000  46,557.5  134,734
--- Filtering for Q20 and Length > 50kb ---
Kept 7428 reads out of 1345923 reads
Stats for Q20_L50k:
file                                                                format  type  num_seqs      sum_len  min_len   avg_len  max_len
05_organelle_genomes/organelle_subset/E_phylacis_Q20_L50k.fastq.gz  FASTQ   DNA      7,428  421,042,149   50,000  56,683.1  134,734

```

I'll try starting with the Q20 40KB subset: `E_phylacis_Q20_L40k.fastq.gz`

## oatk fails with raw reads

Let's try oatk...

First we get the databases:

```bash
db_dir="05_organelle_genomes/db"
mkdir -p ${db_dir}

# Base URL for the raw files
base_url="https://raw.githubusercontent.com/c-zhou/OatkDB/main/v20230921"

# List of files to download
files=(
    "magnoliopsida_mito.fam"
    "magnoliopsida_mito.fam.h3f"
    "magnoliopsida_mito.fam.h3i"
    "magnoliopsida_mito.fam.h3m"
    "magnoliopsida_mito.fam.h3p"
    "magnoliopsida_pltd.fam"
    "magnoliopsida_pltd.fam.h3f"
    "magnoliopsida_pltd.fam.h3i"
    "magnoliopsida_pltd.fam.h3m"
    "magnoliopsida_pltd.fam.h3p"
)

echo "Downloading Magnoliopsida HMM profiles..."

for file in "${files[@]}"; do
    wget -q "${base_url}/${file}" -O "${db_dir}/${file}"
    echo "  Downloaded: ${file}"
done

echo "Done. Databases are in ${db_dir}"
```

Some trial and error with oatk shows that raw nanopore reads don't work. So I'll need to error correct them with hifiasm first.

Various iterations of this approach with different k and c all failed. Just to much noise (homopolymer errors) in the reads I suppose.

```bash
# Define paths
PT_FAM="05_organelle_genomes/db/magnoliopsida_pltd.fam"
MT_FAM="05_organelle_genomes/db/magnoliopsida_mito.fam"
out_dir="05_organelle_genomes/oatk_results"
mkdir -p $out_dir

# Run OatK
# We use the Q20_L40k subset which has 1.6GB of high-quality long reads
oatk -k 1001 \
     -c 3 \
     -t 64 \
     -p $PT_FAM \
     -m $MT_FAM \
     -o ${out_dir}/e_phylacis \
     05_organelle_genomes/organelle_subset/E_phylacis_Q20_L30k.fastq.gz
```

## oatk with hifiasm corrected reads

Let's retry the above with hifiasm corrected reads

First we get the corrected reads

```bash
cd 03_hifiasm_assembly/
hifiasm -o E_phylacis_asm -t 124 --write-ec /dev/null
seqkit stats -T E_phylacis_asm.ec.fa
```

```
file    format  type    num_seqs        sum_len min_len avg_len max_len
E_phylacis_asm.ec.fa    FASTA   DNA     1345923 33991131314     14202   25254.9 147315
```

That's about 28x coverage in these reads (of the haploid genome).

Let's try oatk with c set to ~10x the nuclear coverage to start with (280).

If this fails (too much noise...) then I might need to drop c and/or k

I'll try a few and name the output files differently...

* c280, k1001: run completed, plastid assembled well, mt didn't assemble
* c140, k1001: run has zero coverage arcs, mt genome has 1 contig of 10kb (linear)
* c70, k1001: run has zero coverage arcs, mt genome is in 3 contigs totalling 450kb (one circular, two linear)
* c50, k1001:  run has zero coverage arcs, mt genome is in 3 contigs totalling 450kb (one circular, two linear) basically idential to c70, with an extra 450bp
* c40, k1001:  run has zero coverage arcs, mt genome is in 3 contigs totalling 450kb (one circular, two linear) basically idential to c70, with an extra 450bp
* c30, k1001: as above, but main contig of mt genome increases by 2.7kb
* c25, k1001: this is getting into nuclear genome territory..., but the assembly is identical to c30

Conclusion: c30 k1001 is the best we can probably do. I also tried a few lower values of k, but it made no difference initially, then down at k=601 the mt genome assembly broke, there were two ~identical 70kb linear fragments, and the nv for the main contig went from 11 to 13. 

So we stick with c30 k1001 as the final organelle assemblies.

```bash
# Define paths
PT_FAM="05_organelle_genomes/db/magnoliopsida_pltd.fam"
MT_FAM="05_organelle_genomes/db/magnoliopsida_mito.fam"
out_dir="05_organelle_genomes/oatk_results_c30k1001"
mkdir -p $out_dir

# Run OatK
oatk -k 1001 \
     -c 30 \
     -t 64 \
     -p $PT_FAM \
     -m $MT_FAM \
     -o ${out_dir}/e_phylacis \
     03_hifiasm_assembly/E_phylacis_asm.ec.fa

```

## Final Organelle assemblies

The mitochondrial genome assembles into 3 pieces. One is circular, the other two are linear, but the gfa shows that it's probably a complex mixture of major and minor circles, with potentially some linear pieces in there.

>ctg000001l     length=355146 wlength=41066840.0 nv=11 circular=false path=u151-,u152+,u150-,u148-,u146-,u147+,u141-,u12277+,u139+,u146-,u138-
>ctg000002l     length=23355 wlength=2312145.0 nv=1 circular=false path=u149+
>ctg000003c     length=79427 wlength=7942700.0 nv=1 circular=true path=u12148-

The cp genome assembles into a single circle, and is basically perfect as expected. The gfa shows the typical dumbell structure with LSC, SSC, and IR regions.

>ctg000001c     length=160190 wlength=150009126.0 nv=12 circular=true path=u140+,u12277+,u137-,u7+,u12278-,u143-,u142+,u143+,u12278+,u7-,u137+,u12277-

NB, with a higher c value of 140 or 280 which are the recommended values in the oatk docs of ~5-10x the nuclear coverage, the contig is identical but the graph structure is much simpler with an nv of 4 vs. 12 when c=30; this suggests that with more coverage or fresher samples we could probably do a better job of resolving the mt genome.

# 06 Final Haplotype Assembly

All of the above confirms that the genome is of high quality, that coverage from genome alignment and mapping parental species probes provide useful ways to sort contigs into parental groups, and that the organelle genomes assemble very well. We also know that the maternal haplotype is closer to virginea than decipiens.

Now we need to put this all together, to get the best assembly. A key limitation of most of the above is that I just looked at the biggest 11 contigs, but of course there were plenty of additional big contigs in the haplotype files. The main challenges are to address that, and then do the parental assignment. 

## 1. Re-run hifiasm

Ash Jones helpfully pointed out that hifiasm got a bit confused about the data because it's so heterozygous. Specifically  ```[M::ha_pt_gen] peak_hom: 29; peak_het: -1``` shows it didn't find the het peak, because it thought the hom peak was 29. Actually the homozygous peak is ~60x, and the heterozygous peak is 29x. The kmer distribution just looks odd because the heterozygosity is SO high. 

To see if I can improve the assembly, I'll run three new assemblies, all with `--hom-cov 60` and at all value of `-l [0, 1, 2, 3]`. I can then compare the assemblies, and look at the BUSCO completeness of the two haplotypes in each, as well as the cumulative distribution of contig sizes, and the BUSCO completeness of the major contigs (crossing my fingers that there are only 11 big contigs this time). 


```bash
mkdir -p 06_1_hifiasm_assemblies

# 1. Create a mount point
sudo mkdir -p /mnt/ramdisk

# 2. Mount 1500GB of the 2.2TB RAM as a disk
sudo mount -t tmpfs -o size=1500G tmpfs /mnt/ramdisk
sudo chown $USER /mnt/ramdisk

# 3. Copy your filtered data THERE
cp 02_filtering/E_phylacis_filtered.fastq.gz /mnt/ramdisk/

# 4. Run hifiasm inside the ramdisk
cd /mnt/ramdisk


mkdir l0
mkdir l1
mkdir l2
mkdir l3
mkdir h0 # what happens if we pretend it's haploid...

hifiasm -o l0/E_phylacis -t 160 --ont --hom-cov 60 --hg-size 520m -l 0 --telo-m AAACCCT --dual-scaf E_phylacis_filtered.fastq.gz 2>&1 | tee l0/hifiasm.log
hifiasm -o l1/E_phylacis -t 160 --ont --hom-cov 60 --hg-size 520m -l 1 --telo-m AAACCCT --dual-scaf E_phylacis_filtered.fastq.gz 2>&1 | tee l1/hifiasm.log
hifiasm -o l2/E_phylacis -t 160 --ont --hom-cov 60 --hg-size 520m -l 2 --telo-m AAACCCT --dual-scaf E_phylacis_filtered.fastq.gz 2>&1 | tee l2/hifiasm.log
hifiasm -o l3/E_phylacis -t 160 --ont --hom-cov 60 --hg-size 520m -l 3 --telo-m AAACCCT --dual-scaf E_phylacis_filtered.fastq.gz 2>&1 | tee l3/hifiasm.log
hifiasm -o h0/E_phylacis -t 160 --ont --hom-cov 30 --hg-size 1040m --telo-m AAACCCT E_phylacis_filtered.fastq.gz 2>&1 | tee h0/hifiasm.log



```

Now to pick a winner we can run compleasm and merqury on all of those assemblies. Note that I'll do this either for the primary assembly (for l0 and h0, which don't have haplotype files), or on the h1 and h2 files separately AND together (i.e. cat of the hap1 and hap2 files). This is because there's no guarantee that h1 and h2 will be balanced and sorted correctly. However, their union should have a genome size of ~1040MB, perfect kmer coverage, and 100% BUSCO duplication rate. If I can find a primary assembly or cat'ted h1 and h2 assembly with these characteristics, it's a good one. 

First we convert the gfa files to fasta, and build the meryl database, then run the QC.



```bash
# Build the read k-mer database once
meryl count k=21 threads=160 memory=100G \
    E_phylacis_filtered.fastq.gz output reads.meryl

THREADS=160
LINEAGE="eudicotyledons"
READS_MERYL="/mnt/ramdisk/reads.meryl"

## l0
# Make FASTA; compleasm; merqury
DIR="l0" # for this one, there's only the p_ctg.fa
mkdir -p ${DIR}/qc_results
gfatools gfa2fa ${DIR}/E_phylacis.bp.p_ctg.gfa > ${DIR}/qc_results/p_ctg.fasta
compleasm run -t $THREADS -l $LINEAGE -a ${DIR}/qc_results/p_ctg.fasta -o ${DIR}/qc_results/compleasm_p_ctg
mkdir -p ${DIR}/qc_results/merq_p_ctg && cd ${DIR}/qc_results/merq_p_ctg
merqury.sh $READS_MERYL ../p_ctg.fasta p_ctg_merq
cd ../../..

## l1, 2, and 3, and h0 change dir sequentially and re-run it
DIR="l3" # now we have two haplotypes, which we analyse separately and together

# 1. Prepare FASTAs
mkdir -p ${DIR}/qc_results
gfatools gfa2fa ${DIR}/E_phylacis.bp.hap1.p_ctg.gfa > ${DIR}/qc_results/hap1.fasta
gfatools gfa2fa ${DIR}/E_phylacis.bp.hap2.p_ctg.gfa > ${DIR}/qc_results/hap2.fasta
cat ${DIR}/qc_results/hap1.fasta ${DIR}/qc_results/hap2.fasta > ${DIR}/qc_results/union.fasta

# 2. Compleasm (Hap1, Hap2, then Union)
compleasm run -t $THREADS -l $LINEAGE -a ${DIR}/qc_results/hap1.fasta -o ${DIR}/qc_results/compleasm_hap1
compleasm run -t $THREADS -l $LINEAGE -a ${DIR}/qc_results/hap2.fasta -o ${DIR}/qc_results/compleasm_hap2
compleasm run -t $THREADS -l $LINEAGE -a ${DIR}/qc_results/union.fasta -o ${DIR}/qc_results/compleasm_union

# 3. Merqury (Hap1, Hap2, then Union)
mkdir -p ${DIR}/qc_results/merq_hap1 && cd ${DIR}/qc_results/merq_hap1
merqury.sh $READS_MERYL ../hap1.fasta hap1_merq
cd ../../../

mkdir -p ${DIR}/qc_results/merq_hap2 && cd ${DIR}/qc_results/merq_hap2
merqury.sh $READS_MERYL ../hap2.fasta hap2_merq
cd ../../../

mkdir -p ${DIR}/qc_results/merq_union && cd ${DIR}/qc_results/merq_union
merqury.sh $READS_MERYL ../union.fasta union_merq
cd ../../../

# Then repeat the above for h0, l2, and l3


#########
### Copy things you want back to the hard drive
### Ignore all the big files, but of course keep the assemblies.

# 5. VERY IMPORTANT: Move results back to /home before rebooting!
rsync -avh --progress --exclude='*.fastq.gz' /mnt/ramdisk/ ~/MM_genome/06_1_hifiasm_assemblies/

cd ~/MM_genome

# 6. Unmount (only after you've verified the files are safe!)
# cd ~
# sudo umount /mnt/ramdisk

# ignore this stuff for git, mostly
echo "06_1_hifiasm_assemblies/" >> .gitignore

```

Get assembly basic stats and organise them:

```bash
# get all the .fasta except those in the compleasm folders, and get lots of stats
find . -name "*.fasta" | grep -E "hap1|hap2|union|p_ctg" | grep -vE "compleasm|merq" | xargs seqkit stats -a -j 160 -N 50,90,95,99 | column -t
```

```
file                         format  type  num_seqs  sum_len        min_len  avg_len       max_len     Q1        Q2       Q3            sum_gap  N50         N50_num  Q20(%)  Q30(%)  AvgQual  GC(%)  sum_n   N50         N90         N95         N99
./l3/qc_results/union.fasta  FASTA   DNA   131       1,097,940,037  30,917   8,381,221.7   65,688,633  66,540.5  92,347   258,145       0        45,062,906  10       0       0       0        39.63  43,947  45,062,906  33,511,041  20,522,527  504,488
./l3/qc_results/hap2.fasta   FASTA   DNA   45        573,953,311    42,116   12,754,518    65,688,633  66,777    151,689  26,036,980    0        43,652,620  6        0       0       0        39.63  43,747  43,652,620  26,036,980  20,522,527  14,252,391
./l3/qc_results/hap1.fasta   FASTA   DNA   86        523,986,726    30,917   6,092,868.9   61,716,932  65,110    85,704   191,335       0        45,216,074  5        0       0       0        39.64  200     45,216,074  33,914,441  33,914,441  191,335
./l2/qc_results/union.fasta  FASTA   DNA   132       1,213,886,437  24,825   9,196,109.4   65,674,939  67,058    94,767   338,826.5     0        45,708,075  11       0       0       0        39.61  300     45,708,075  33,387,322  20,442,905  13,773,696
./l2/qc_results/hap2.fasta   FASTA   DNA   58        610,625,537    42,116   10,528,026.5  61,716,932  76,003    118,755  380,136       0        45,708,075  6        0       0       0        39.61  200     45,708,075  33,387,322  26,036,980  20,442,905
./l2/qc_results/hap1.fasta   FASTA   DNA   74        603,260,900    24,825   8,152,174.3   65,674,939  56,568    88,451   293,829       0        45,081,242  6        0       0       0        39.62  100     45,081,242  33,719,123  16,961,878  407,150
./l1/qc_results/union.fasta  FASTA   DNA   119       1,069,391,022  30,917   8,986,479.2   65,674,939  64,743    87,899   303,443.5     0        45,081,242  10       0       0       0        39.62  74,339  45,081,242  33,719,123  20,310,318  13,773,696
./l1/qc_results/hap2.fasta   FASTA   DNA   45        466,004,620    35,842   10,355,658.2  61,716,932  68,000    94,551   550,701       0        41,950,308  5        0       0       0        39.64  100     41,950,308  33,465,692  20,310,318  20,310,318
./l1/qc_results/hap1.fasta   FASTA   DNA   74        603,386,402    30,917   8,153,870.3   65,674,939  56,855    86,360   232,920       0        45,081,242  6        0       0       0        39.61  74,239  45,081,242  33,719,123  16,961,878  518,253
./l0/qc_results/p_ctg.fasta  FASTA   DNA   79        1,014,773,291  30,917   12,845,231.5  65,674,939  69,098.5  142,767  33,781,247.5  0        43,652,620  10       0       0       0        39.61  0       43,652,620  33,930,661  20,458,351  13,773,696
```



## L1 compleasm and kmers

For all of the following, run this code from within the relevant in the qc_results folder for each assembly
```
echo "| Metric | Hap1 | Hap2 | Union |"
echo "| :--- | :--- | :--- | :--- |"

# Create temporary files for each column to handle the side-by-side join
for mode in hap1 hap2 union; do
    grep -P "S:|D:|F:|M:|N:" compleasm_${mode}/summary.txt | awk -F':|,' '{print $2}' > ${mode}_col.tmp
done

# Add labels and join the columns
paste -d '|' <(echo -e "S\nD\nF\nM\nN") hap1_col.tmp hap2_col.tmp union_col.tmp | sed 's/^/| /' | sed 's/$/ |/'

# Clean up
rm *.tmp

echo "| Metric | Hap1 | Hap2 | Union |"
echo "| :--- | :--- | :--- | :--- |"

# Extract QV (Column 4 of the .qv file)
qv_h1=$(awk '{print $4}' merq_hap1/hap1_merq.qv)
qv_h2=$(awk '{print $4}' merq_hap2/hap2_merq.qv)
qv_un=$(awk '{print $4}' merq_union/union_merq.qv)

# Extract Completeness (Column 5 of the generic completeness.stats file)
comp_h1=$(awk '{print $5}' merq_hap1/completeness.stats)
comp_h2=$(awk '{print $5}' merq_hap2/completeness.stats)
comp_un=$(awk '{print $5}' merq_union/completeness.stats)

echo "| **QV** | $qv_h1 | $qv_h2 | $qv_un |"
echo "| **K-mer Completeness %** | $comp_h1 | $comp_h2 | $comp_un |"

```

The h0 and l0 assemblies don't separate haplotypes, so I'll only look at L1, L2 and L3. 

### L1

| Metric | Hap1 | Hap2 | Union |
| :--- | :--- | :--- | :--- |
| S|78.36%|78.89%|2.92% |
| D|17.15%|6.35%|96.86% |
| F|0.71%|1.00%|0.18% |
| M|3.78%|13.76%|0.04% |
| N|2805|2805|2805 |
| **QV** | 64.9563 | 66.1094 | 65.4217 |
| **K-mer Completeness %** | 66.5101 | 54.6202 | 99.3052 |


### L2
| Metric | Hap1 | Hap2 | Union |
| :--- | :--- | :--- | :--- |
| S|78.36%|73.51%|2.71% |
| D|17.15%|19.00%|97.08% |
| F|0.71%|0.46%|0.18% |
| M|3.78%|7.02%|0.04% |
| N|2805|2805|2805 |
| **QV** | 65.0701 | 65.76 | 65.4034 |
| **K-mer Completeness %** | 66.5072 | 67.0275 | 99.3053 |

### L3
| Metric | Hap1 | Hap2 | Union |
| :--- | :--- | :--- | :--- |
| S|96.83%|91.69%|2.85% |
| D|2.00%|7.95%|96.93% |
| F|0.29%|0.32%|0.18% |
| M|0.89%|0.04%|0.04% |
| N|2805|2805|2805 |
| **QV** | 65.1172 | 65.6202 | 65.3729 |
| **K-mer Completeness %** | 62.3954 | 66.0179 | 99.3053 |

### Quast

```
quast.py l1/qc_results/hap1.fasta l1/qc_results/hap2.fasta \
    l2/qc_results/hap1.fasta l2/qc_results/hap2.fasta \
    l3/qc_results/hap1.fasta l2/qc_results/hap2.fasta \
    -o quast_comparison \
    --labels L1-hap1,L1-hap2,L2-hap1,L2-hap2,L3-hap1,L3-hap2
```
![Screenshot 2026-02-06 at 6 50 35 AM](https://github.com/user-attachments/assets/e4c27269-182b-4f78-8a26-4f09d078c4b4)
![Screenshot 2026-02-06 at 6 53 00 AM](https://github.com/user-attachments/assets/6f2d5c4f-1ff1-43d4-ba0c-73f4533f90d5)
![Screenshot 2026-02-06 at 6 53 15 AM](https://github.com/user-attachments/assets/9f165230-6466-4238-a15e-c16163461aa9)

Most of the length of the L3 haplotypes is contained in the first 12 (hap1) or 14 (hap2) contigs. There are essentially no N's. The union of hap1 and hap2 is 1.135Gb which is ~568Mb haploid genome size. This is a little big, so it seems likely that there's some duplication, especially in hap2 whose total size is 611Mb. This is confirmed by the ~8% BUSCO duplication rate in hap2, and is something to look for when cleaning the assembly.

### Conclusion: use L3

L3 is *clearly* the best assembly. The two haplotypes have incredibly high compleasm scores, little duplication, almost nothing missing. I suspect hap2 has a a few contigs that actually belong in hap1, (it's longer than hap1, and it has more duplicates), but parental binning will solve that in the next step. The stats on the union are also excellent - 97% of BUSCOs are duplicated, as expected. The QV scores are uniformly high. The kmer completeness is very very high for the Union graph (bearing in mind these are nanopore reads) and is just about right for the haplotypes (we expect a lot of missing kmers if it's an F1 hybrid of divergent species). The N50 is about 45MB for all three assemblies (union, H1, H2). 

Now I'll tidy up and get what I need off the ramdisk. No point keeping big files here, especially since the assemblies take only an hour to re-do.

```bash
cd /mnt/ramdisk

#get rid of big duplicate files I don't need
find . -name "E_phylacis_filtered.fastq.gz" -delete
find . -name "*ec.bin" -delete
rm -rf mb_downloads/

# move to hard drive
rsync -avP . ~/MM_genome/06_1_hifiasm_assemblies/

# check sizes to ~verify
du -sh .
du -sh ~/MM_genome/06_1_hifiasm_assemblies/

# unmount and delete mount point
cd ~/MM_genome/
sudo umount /mnt/ramdisk
sudo rmdir /mnt/ramdisk
```

## L3 additional QC

## What does it look like

```bash
# 1. Get stats, 2. Sort by length, 3. Calculate cumulative sum and proportion
seqkit fx2tab -n -l -g 06_1_hifiasm_assemblies/l3/E_phylacis.hap1.fasta | \
sort -k2,2nr | \
awk 'BEGIN {OFS="\t"; total='$(seqkit stats -T 06_1_hifiasm_assemblies/l3/E_phylacis.hap1.fasta | awk 'NR==2 {print $5}')'} \
{sum+=$2; print $1, $2, $3, sum/total}'


seqkit fx2tab -n -l -g 06_1_hifiasm_assemblies/l3/E_phylacis.hap2.fasta | \
sort -k2,2nr | \
awk 'BEGIN {OFS="\t"; total='$(seqkit stats -T 06_1_hifiasm_assemblies/l3/E_phylacis.hap2.fasta | awk 'NR==2 {print $5}')'} \
{sum+=$2; print $1, $2, $3, sum/total}'
```

```
Contig          Size            GC      Cumulative proportion            
h1tg000010l     61716932        39.49   0.117783
h1tg000006l     60158556        39.36   0.232593
h1tg000002l     54083987        39.81   0.335809
h1tg000001l     51333920        39.26   0.433777
h1tg000014l     45216074        39.24   0.520069
h1tg000007l     41950308        39.66   0.600129
h1tg000008l     39810478        39.63   0.676105
h1tg000003l     38464560        39.82   0.749513
h1tg000005l     36525698        39.68   0.81922
h1tg000011l     35389336        39.37   0.886759
h1tg000009l     33914441        39.76   0.951483
h1tg000004l     16945854        39.58   0.983823
h1tg000021l     550977  53.13   0.984874
h1tg000019l     550701  38.47   0.985925
h1tg000027l     504488  37.83   0.986888
h1tg000018l     267638  44.29   0.987399
h1tg000031l     262582  40.03   0.9879
h1tg000017l     247779  52.98   0.988373
h1tg000044l     246068  42.27   0.988843
h1tg000015l     240840  53.76   0.989302
h1tg000045l     202662  54.80   0.989689
h1tg000023l     191335  38.27   0.990054
h1tg000012l     179571  54.71   0.990397
h1tg000047l     173338  54.33   0.990728
h1tg000040l     164941  54.08   0.991042
h1tg000066l     162610  39.77   0.991353
h1tg000046l     154476  44.83   0.991648
h1tg000060l     128970  42.15   0.991894
h1tg000059l     120315  39.46   0.992123
h1tg000071l     118804  44.88   0.99235
h1tg000051l     116421  55.11   0.992572
h1tg000024l     110453  55.37   0.992783
h1tg000058l     108016  51.82   0.992989
h1tg000035l     106259  52.63   0.993192
h1tg000025l     94874   54.08   0.993373
h1tg000043l     94628   38.58   0.993554
h1tg000085l     94551   39.56   0.993734
h1tg000016l     93941   45.54   0.993913
h1tg000028l     88856   35.42   0.994083
h1tg000063l     87206   41.05   0.994249
h1tg000038l     87121   39.31   0.994416
h1tg000075l     86414   55.36   0.99458
h1tg000041l     86306   55.35   0.994745
h1tg000056l     85102   54.47   0.994908
h1tg000029l     84892   53.00   0.99507
h1tg000064l     82462   38.37   0.995227
h1tg000050l     81571   38.80   0.995383
h1tg000057l     81096   53.93   0.995537
h1tg000076l     80866   51.88   0.995692
h1tg000034l     80746   41.93   0.995846
h1tg000055l     80201   39.46   0.995999
h1tg000022l     79850   37.07   0.996151
h1tg000036l     78549   37.35   0.996301
h1tg000048l     77500   37.65   0.996449
h1tg000084l     75001   38.03   0.996592
h1tg000061l     73615   38.31   0.996733
h1tg000033l     72644   54.23   0.996871
h1tg000070l     71332   39.67   0.997008
h1tg000073l     71031   39.25   0.997143
h1tg000049l     68892   40.66   0.997275
h1tg000039l     68546   39.76   0.997405
h1tg000042l     68000   37.92   0.997535
h1tg000079l     67339   41.26   0.997664
h1tg000032l     67171   41.50   0.997792
h1tg000072l     65110   38.42   0.997916
h1tg000080l     64287   40.53   0.998039
h1tg000020l     62283   55.10   0.998158
h1tg000065l     61333   54.26   0.998275
h1tg000078l     60163   41.59   0.99839
h1tg000052l     59545   40.02   0.998503
h1tg000068l     57909   38.35   0.998614
h1tg000086l     57449   41.38   0.998723
h1tg000013l     56568   55.01   0.998831
h1tg000054l     54413   55.10   0.998935
h1tg000037l     54248   53.39   0.999039
h1tg000030l     53716   41.06   0.999141
h1tg000077l     53318   55.54   0.999243
h1tg000083l     53008   38.12   0.999344
h1tg000053l     50303   30.76   0.99944
h1tg000067l     49495   40.33   0.999535
h1tg000062l     49067   53.80   0.999628
h1tg000026l     48552   36.93   0.999721
h1tg000069l     46567   54.44   0.99981
h1tg000074l     35174   52.33   0.999877
h1tg000081l     33610   55.48   0.999941
h1tg000082l     30917   55.26   1
```


```
Contig          Size            GC      Cumulative proportion            
h2tg000005l     65688633        39.49   0.114449
h2tg000001l     59206199        39.36   0.217605
h2tg000008l     56293706        39.33   0.315685
h2tg000011l     51100880        39.83   0.404718
h2tg000004l     45062906        39.29   0.483232
h2tg000012l     43652620        39.65   0.559288
h2tg000006l     42555022        39.62   0.633431
h2tg000009l     39808945        39.70   0.70279
h2tg000003l     37845087        39.39   0.768728
h2tg000007l     34610167        39.91   0.829029
h2tg000010l     33511041        39.69   0.887416
h2tg000002l     26036980        39.80   0.93278
h2tg000014l     20522527        40.12   0.968536
h2tg000013l     14252391        39.71   0.993368
h2tg000015l     500252  53.89   0.99424
h2tg000023l     280860  53.22   0.994729
h2tg000040l     253708  53.17   0.995171
h2tg000028l     227859  54.63   0.995568
h2tg000018l     207414  45.42   0.99593
h2tg000016l     193346  41.89   0.996267
h2tg000027l     174770  51.63   0.996571
h2tg000019l     164625  54.10   0.996858
h2tg000017l     151689  37.21   0.997122
h2tg000031l     140617  44.85   0.997367
h2tg000034l     118804  44.88   0.997574
h2tg000045l     95361   39.89   0.99774
h2tg000032l     93954   52.39   0.997904
h2tg000021l     92347   38.27   0.998065
h2tg000026l     89563   53.97   0.998221
h2tg000020l     87899   38.96   0.998374
h2tg000038l     84936   39.18   0.998522
h2tg000043l     75596   36.81   0.998654
h2tg000030l     68747   52.26   0.998774
h2tg000035l     66777   43.49   0.99889
h2tg000033l     66304   37.50   0.999006
h2tg000044l     66155   37.62   0.999121
h2tg000037l     63729   41.59   0.999232
h2tg000042l     63218   39.06   0.999342
h2tg000024l     63097   37.31   0.999452
h2tg000041l     62567   55.42   0.999561
h2tg000039l     57976   39.97   0.999662
h2tg000025l     56855   40.42   0.999761
h2tg000022l     49225   39.34   0.999847
h2tg000029l     45841   40.27   0.999927
h2tg000036l     42116   54.81   1
```

## Unresolved bubbles
One question I had is whether there were unresolved bubbles in the unitig graph, but there aren't any at all. (`gfatools bubble E_phylacis.bp.p_utg.gfa > bubbles_l3_utg.bed` produces nothing). This suggests that all the bubbles were big and obvious, as expected for accurate long reads with very divergent haplotypes.

## Kmer spectra
Next are kmer spectra for the union, h1, and h2 from the merqury runs above. These are also incredibly good. The union plot has very clear 1x and 2x peaks, as expected. The h1 and h2 plots show just over half of the 1x kmers in each sample (as expected - many unique kmers), and a smaller number of 2x kmers (as expected - some homozygous regions). 

<table>
  <tr>
    <td align="center"><b>Union (Total)</b></td>
    <td align="center"><b>Haplotype 1</b></td>
    <td align="center"><b>Haplotype 2</b></td>
  </tr>
  <tr>
    <td><img src="https://github.com/user-attachments/assets/4a75b606-59ef-4e56-9d7a-411907c19c9a" width="100%" /></td>
    <td><img src="https://github.com/user-attachments/assets/41fdd64e-c2c0-4b36-b458-7cf0452297b5" width="100%" /></td>
    <td><img src="https://github.com/user-attachments/assets/4357dad5-cc60-4f3c-b415-1f7dcfa8533e" width="100%" /></td>
  </tr>
</table>


## Dotplots and contig level chromosome homology

Let's align the hap1 and hap2 files to the virginea and decipiens references, and pull out the best contigs in those file that they align to. This will help us sort things out later.

```
# Define paths
HAP1="06_1_hifiasm_assemblies/l3/E_phylacis.hap1.fasta"
HAP2="06_1_hifiasm_assemblies/l3/E_phylacis.hap2.fasta"
VIR_REF="parental_spp_genomes/E_virginea.fa"
DEC_REF="parental_spp_genomes/E_decipiens.fa"

OUT_DIR="06_1_hifiasm_assemblies/l3/genome_alignment"
mkdir -p $OUT_DIR

# 1. Run Alignments
echo "Running Minimap2 alignments..."

# Hap1 vs Parents
minimap2 -x asm5 -N 1000 --secondary=no -t 128 $VIR_REF $HAP1 > $OUT_DIR/hap1_vs_virginea.paf
minimap2 -x asm5 -N 1000 --secondary=no -t 128 $DEC_REF $HAP1 > $OUT_DIR/hap1_vs_decipiens.paf

# Hap2 vs Parents
minimap2 -x asm5 -N 1000 --secondary=no -t 128 $VIR_REF $HAP2 > $OUT_DIR/hap2_vs_virginea.paf
minimap2 -x asm5 -N 1000 --secondary=no -t 128 $DEC_REF $HAP2 > $OUT_DIR/hap2_vs_decipiens.paf

# CREATE HEADER MAPS
# Regex updated: Captures all text after the CCA isolate ID up until the first comma.
echo "Creating header mapping files..."
grep ">" "$VIR_REF" | sed -E 's/^>([^ ]+) .*isolate CCA[0-9]+ ([^,]+).*/\1\t\2/' > "$OUT_DIR/vir_map.tsv"
grep ">" "$DEC_REF" | sed -E 's/^>([^ ]+) .*isolate CCA[0-9]+ ([^,]+).*/\1\t\2/' > "$OUT_DIR/dec_map.tsv"

# GENERATE THE TABLES
for HAP in "hap1" "hap2"; do
    echo "Processing $HAP..."
    TABLE="$OUT_DIR/${HAP}_homology_table.tsv"
    
    echo -e "Your_Contig\tYour_Contig_Len\tvir_Ref_Chr1\tvir_MatchBases_Chr1\tvir_FullRefLen_Chr1\tvir_Ref_Chr2\tvir_MatchBases_Chr2\tvir_FullRefLen_Chr2\tvir_Ref_Chr3\tvir_MatchBases_Chr3\tvir_FullRefLen_Chr3\tdec_Ref_Chr1\tdec_MatchBases_Chr1\tdec_FullRefLen_Chr1\tdec_Ref_Chr2\tdec_MatchBases_Chr2\tdec_FullRefLen_Chr2\tdec_Ref_Chr3\tdec_MatchBases_Chr3\tdec_FullRefLen_Chr3" > "$TABLE"

    CONTIGS=$(cat "$OUT_DIR/${HAP}_vs_virginea.paf" "$OUT_DIR/${HAP}_vs_decipiens.paf" | cut -f1 | sort -u)

    for CONTIG in $CONTIGS; do
        # Get contig length (col 2) correctly
        QLEN=$(awk -v c="$CONTIG" '$1==c {print $2; exit}' "$OUT_DIR/${HAP}_vs_virginea.paf" "$OUT_DIR/${HAP}_vs_decipiens.paf")
        
        LINE="$CONTIG\t$QLEN"
        
        for REF_TYPE in "virginea" "decipiens"; do
            PAF="$OUT_DIR/${HAP}_vs_${REF_TYPE}.paf"
            MAP_ID=$(echo $REF_TYPE | cut -c1-3)
            MAP="$OUT_DIR/${MAP_ID}_map.tsv"
            
            # Sum column 10, then join with map
            # We use | as a temp separator to keep the name and length together during sorting
            HITS=$(awk -F'\t' -v c="$CONTIG" '
                NR==FNR { d[$1]=$2; next }
                $1==c { m[$6]+=$10; r[$6]=$11 }
                END {
                    for (i in m) {
                        print m[i] "###" i ", " d[i] "###" r[i]
                    }
                }' "$MAP" "$PAF" | sort -t'#' -k1,1nr | head -n 3)

            HIT_COUNT=$(echo "$HITS" | grep -v '^$' | wc -l)
            
            while read -r HIT; do
                [ -z "$HIT" ] && continue
                # Reformat the "Sum###Name###Len" back to tabs
                VALS=$(echo "$HIT" | sed 's/###/\t/g')
                # The columns we want are Name (2), Sum (1), Len (3)
                LINE+="\t$(echo "$VALS" | cut -f2)\t$(echo "$VALS" | cut -f1)\t$(echo "$VALS" | cut -f3)"
            done <<< "$HITS"
            
            for ((i=HIT_COUNT; i<3; i++)); do
                LINE+="\tNA\t0\t0"
            done
        done
        echo -e "$LINE" >> "$TABLE"
    done
done

echo "Success! Final tables with descriptive names are in $OUT_DIR"

# now we make percentage tables

for HAP in "hap1" "hap2"; do
    INPUT="$OUT_DIR/${HAP}_homology_table.tsv"
    OUTPUT="$OUT_DIR/${HAP}_homology_percentages.tsv"
    
    echo "Refactoring $HAP into percentages..."

    # Create new header
    echo -e "Your_Contig\tYour_Contig_Len\tvir_Ref_Chr1\t%matchChr1\tvir_Ref_Chr2\t%matchChr2\tvir_Ref_Chr3\t%matchChr3\tdec_Ref_Chr1\t%matchChr1\tdec_Ref_Chr2\t%matchChr2\tdec_Ref_Chr3\t%matchChr3" > "$OUTPUT"

    # Skip header (NR>1) and process rows
    awk -F'\t' 'NR > 1 {
        qlen = $2
        
        # Function to format percentage
        # pct = (matches / qlen) * 100
        
        printf "%s\t%s", $1, $2
        
        # Virginea columns: Names are 3, 6, 9. Matches are 4, 7, 10.
        for (i=3; i<=9; i+=3) {
            pct = (qlen > 0) ? ($ (i+1) / qlen) * 100 : 0
            printf "\t%s\t%.1f%%", $i, pct
        }
        
        # Decipiens columns: Names are 12, 15, 18. Matches are 13, 16, 19.
        for (i=12; i<=18; i+=3) {
            pct = (qlen > 0) ? ($ (i+1) / qlen) * 100 : 0
            printf "\t%s\t%.1f%%", $i, pct
        }
        
        printf "\n"
    }' "$INPUT" >> "$OUTPUT"

done
```

This code just looks at the genome alignments, and asks what each contig aligns to best in virginea and decipiens, then summarises the output. 

### Homology of contigs vs. virginea and decipiens

Easiest to view it like this:

```bash
# For Hap1
(head -n 1 06_1_hifiasm_assemblies/l3/genome_alignment/hap1_homology_percentages.tsv && \
 tail -n +2 06_1_hifiasm_assemblies/l3/genome_alignment/hap1_homology_percentages.tsv | sort -t$'\t' -k2,2rn) | \
 column -t -s $'\t'

# For Hap2
(head -n 1 06_1_hifiasm_assemblies/l3/genome_alignment/hap2_homology_percentages.tsv && \
 tail -n +2 06_1_hifiasm_assemblies/l3/genome_alignment/hap2_homology_percentages.tsv | sort -t$'\t' -k2,2rn) | \
 column -t -s $'\t'
```

```
Your_Contig  Your_Contig_Len  vir_Ref_Chr1                    %matchChr1  vir_Ref_Chr2               %matchChr2  vir_Ref_Chr3               %matchChr3  dec_Ref_Chr1                    %matchChr1  dec_Ref_Chr2               %matchChr2  dec_Ref_Chr3               %matchChr3
h1tg000010l  61716932         CM024521.1, chromosome 8        52.9%       CM024520.1, chromosome 7   4.2%        CM024517.1, chromosome 4   1.7%        CM024616.1, chromosome 8        31.5%       CM024614.1, chromosome 6   1.0%        CM024617.1, chromosome 9   0.4%
h1tg000006l  60158556         CM024516.1, chromosome 3        24.5%       CM024523.1, chromosome 10  1.9%        CM024518.1, chromosome 5   0.7%        CM024611.1, chromosome 3        45.7%       CM024614.1, chromosome 6   0.4%        CM024613.1, chromosome 5   0.2%
h1tg000002l  54083987         CM024519.1, chromosome 6        34.0%       CM024524.1, chromosome 11  0.4%        CM024518.1, chromosome 5   0.1%        CM024614.1, chromosome 6        60.2%       CM024615.1, chromosome 7   1.4%        CM024619.1, chromosome 11  0.8%
h1tg000001l  51333920         CM024520.1, chromosome 7        27.8%       CM024515.1, chromosome 2   0.7%        CM024521.1, chromosome 8   0.5%        CM024615.1, chromosome 7        47.7%       CM024613.1, chromosome 5   1.2%        CM024617.1, chromosome 9   1.1%
h1tg000014l  45216074         CM024518.1, chromosome 5        25.8%       CM024519.1, chromosome 6   1.4%        CM024520.1, chromosome 7   0.3%        CM024613.1, chromosome 5        46.7%       CM024616.1, chromosome 8   0.6%        CM024614.1, chromosome 6   0.5%
h1tg000007l  41950308         CM024514.1, chromosome 1        60.7%       CM024516.1, chromosome 3   3.6%        CM024520.1, chromosome 7   0.1%        CM024609.1, chromosome 1        33.8%       CM024615.1, chromosome 7   0.1%        CM024611.1, chromosome 3   0.1%
h1tg000008l  39810478         CM024524.1, chromosome 11       65.0%       CM024520.1, chromosome 7   1.3%        CM024519.1, chromosome 6   0.1%        CM024619.1, chromosome 11       37.0%       CM024617.1, chromosome 9   0.9%        CM024611.1, chromosome 3   0.2%
h1tg000003l  38464560         CM024523.1, chromosome 10       32.0%       CM024520.1, chromosome 7   1.0%        CM024524.1, chromosome 11  0.7%        CM024618.1, chromosome 10       57.2%       CM024612.1, chromosome 4   0.9%        CM024614.1, chromosome 6   0.7%
h1tg000005l  36525698         CM024522.1, chromosome 9        56.4%       CM024524.1, chromosome 11  6.7%        CM024520.1, chromosome 7   0.1%        CM024617.1, chromosome 9        34.9%       CM024612.1, chromosome 4   0.5%        CM024611.1, chromosome 3   0.2%
h1tg000011l  35389336         CM024517.1, chromosome 4        54.2%       CM024521.1, chromosome 8   4.4%        CM024523.1, chromosome 10  2.2%        CM024612.1, chromosome 4        34.9%       CM024614.1, chromosome 6   0.2%        CM024615.1, chromosome 7   0.2%
h1tg000009l  33914441         CM024515.1, chromosome 2        32.7%       CM024516.1, chromosome 3   1.5%        CM024518.1, chromosome 5   0.1%        CM024610.1, chromosome 2        57.1%       CM024609.1, chromosome 1   0.5%        CM024613.1, chromosome 5   0.4%
h1tg000004l  16945854         CM024515.1, chromosome 2        55.5%       CM024518.1, chromosome 5   0.0%        CM024523.1, chromosome 10  0.0%        CM024610.1, chromosome 2        29.5%       CM024615.1, chromosome 7   0.0%        CM024616.1, chromosome 8   0.0%
h1tg000021l  550977           CM024515.1, chromosome 2        79.9%       CM024516.1, chromosome 3   1.4%        NA                         0.0%        CM024610.1, chromosome 2        64.3%       CM024614.1, chromosome 6   0.6%        CM024617.1, chromosome 9   0.0%
h1tg000019l  550701           CM024518.1, chromosome 5        24.6%       CM024516.1, chromosome 3   1.0%        CM024522.1, chromosome 9   0.8%        CM024611.1, chromosome 3        26.4%       CM024612.1, chromosome 4   0.3%        CM024614.1, chromosome 6   0.0%
h1tg000027l  504488           CM024519.1, chromosome 6        1.3%        CM024518.1, chromosome 5   0.7%        CM024521.1, chromosome 8   0.6%        CM024617.1, chromosome 9        2.4%        CM024618.1, chromosome 10  1.6%        CM024614.1, chromosome 6   1.0%
h1tg000018l  267638           CM024524.1, chromosome 11       13.5%       CM024525.1, chloroplast    11.8%       CM024517.1, chromosome 4   3.0%        CM024611.1, chromosome 3        13.2%       CM024614.1, chromosome 6   8.9%        CM024615.1, chromosome 7   4.6%
h1tg000031l  262582           CM024524.1, chromosome 11       10.5%       CM024521.1, chromosome 8   6.1%        CM024520.1, chromosome 7   4.9%        CM024619.1, chromosome 11       90.8%       NA                         0.0%        NA                         0.0%
h1tg000017l  247779           CM024515.1, chromosome 2        74.3%       CM024516.1, chromosome 3   8.0%        NA                         0.0%        CM024610.1, chromosome 2        52.8%       CM024614.1, chromosome 6   3.4%        NA                         0.0%
h1tg000044l  246068           CM024515.1, chromosome 2        1.1%        CM024522.1, chromosome 9   1.1%        CM024518.1, chromosome 5   0.7%        JABKBM010000014.1, tig00112428  3.4%        CM024613.1, chromosome 5   1.3%        CM024611.1, chromosome 3   1.0%
h1tg000015l  240840           CM024515.1, chromosome 2        74.1%       CM024516.1, chromosome 3   1.2%        NA                         0.0%        CM024610.1, chromosome 2        20.3%       CM024614.1, chromosome 6   1.8%        NA                         0.0%
h1tg000045l  202662           CM024515.1, chromosome 2        85.8%       CM024516.1, chromosome 3   1.3%        NA                         0.0%        CM024610.1, chromosome 2        99.8%       CM024613.1, chromosome 5   0.0%        NA                         0.0%
h1tg000023l  191335           CM024525.1, chloroplast         86.5%       CM024514.1, chromosome 1   2.7%        CM024521.1, chromosome 8   0.9%        CM024611.1, chromosome 3        83.8%       CM024617.1, chromosome 9   2.8%        CM024609.1, chromosome 1   2.7%
h1tg000012l  179571           CM024515.1, chromosome 2        82.2%       CM024516.1, chromosome 3   3.3%        NA                         0.0%        CM024610.1, chromosome 2        60.7%       NA                         0.0%        NA                         0.0%
h1tg000047l  173338           CM024515.1, chromosome 2        88.1%       NA                         0.0%        NA                         0.0%        CM024610.1, chromosome 2        84.5%       CM024614.1, chromosome 6   0.1%        CM024619.1, chromosome 11  0.0%
h1tg000040l  164941           CM024515.1, chromosome 2        59.6%       NA                         0.0%        NA                         0.0%        CM024610.1, chromosome 2        80.7%       CM024614.1, chromosome 6   0.2%        NA                         0.0%
h1tg000066l  162610           CM024518.1, chromosome 5        15.5%       CM024514.1, chromosome 1   3.0%        CM024519.1, chromosome 6   2.2%        JABKBM010000015.1, tig00006793  52.6%       CM024613.1, chromosome 5   33.1%       CM024619.1, chromosome 11  0.1%
h1tg000046l  154476           CM024516.1, chromosome 3        12.6%       CM024521.1, chromosome 8   5.0%        CM024519.1, chromosome 6   4.1%        CM024609.1, chromosome 1        12.6%       CM024613.1, chromosome 5   9.3%        CM024617.1, chromosome 9   8.9%
h1tg000060l  128970           CM024518.1, chromosome 5        0.4%        CM024517.1, chromosome 4   0.2%        CM024521.1, chromosome 8   0.2%        JABKBM010000014.1, tig00112428  3.9%        CM024618.1, chromosome 10  0.8%        CM024614.1, chromosome 6   0.1%
h1tg000059l  120315           CM024520.1, chromosome 7        9.8%        CM024518.1, chromosome 5   7.9%        CM024517.1, chromosome 4   6.2%        CM024615.1, chromosome 7        9.3%        CM024612.1, chromosome 4   8.0%        CM024613.1, chromosome 5   7.2%
h1tg000071l  118804           CM024517.1, chromosome 4        25.7%       CM024523.1, chromosome 10  13.2%       CM024516.1, chromosome 3   8.0%        CM024614.1, chromosome 6        49.7%       CM024613.1, chromosome 5   5.7%        CM024611.1, chromosome 3   2.7%
h1tg000024l  110453           CM024515.1, chromosome 2        82.3%       CM024516.1, chromosome 3   0.3%        NA                         0.0%        CM024610.1, chromosome 2        108.5%      NA                         0.0%        NA                         0.0%
h1tg000058l  108016           CM024515.1, chromosome 2        53.6%       CM024516.1, chromosome 3   1.5%        CM024520.1, chromosome 7   0.1%        CM024610.1, chromosome 2        57.3%       CM024614.1, chromosome 6   0.8%        NA                         0.0%
h1tg000035l  106259           CM024516.1, chromosome 3        63.4%       CM024515.1, chromosome 2   23.5%       CM024514.1, chromosome 1   0.7%        CM024610.1, chromosome 2        39.3%       CM024614.1, chromosome 6   2.3%        CM024618.1, chromosome 10  1.6%
h1tg000025l  94874            CM024515.1, chromosome 2        60.2%       CM024516.1, chromosome 3   8.0%        NA                         0.0%        CM024610.1, chromosome 2        87.9%       NA                         0.0%        NA                         0.0%
h1tg000043l  94628            CM024525.1, chloroplast         87.0%       CM024516.1, chromosome 3   4.8%        CM024517.1, chromosome 4   3.4%        CM024611.1, chromosome 3        89.1%       CM024613.1, chromosome 5   4.0%        CM024618.1, chromosome 10  0.5%
h1tg000085l  94551            CM024514.1, chromosome 1        15.9%       NA                         0.0%        NA                         0.0%        CM024609.1, chromosome 1        73.3%       NA                         0.0%        NA                         0.0%
h1tg000016l  93941            CM024514.1, chromosome 1        9.5%        CM024524.1, chromosome 11  5.9%        CM024516.1, chromosome 3   5.4%        CM024617.1, chromosome 9        18.8%       CM024614.1, chromosome 6   6.6%        CM024616.1, chromosome 8   4.3%
h1tg000028l  88856            CM024524.1, chromosome 11       40.6%       CM024516.1, chromosome 3   4.5%        CM024514.1, chromosome 1   1.5%        CM024609.1, chromosome 1        61.8%       NA                         0.0%        NA                         0.0%
h1tg000063l  87206            CM024514.1, chromosome 1        32.9%       CM024520.1, chromosome 7   7.9%        CM024516.1, chromosome 3   1.6%        CM024609.1, chromosome 1        63.5%       CM024619.1, chromosome 11  0.0%        NA                         0.0%
h1tg000038l  87121            CM024520.1, chromosome 7        9.8%        CM024518.1, chromosome 5   9.2%        CM024517.1, chromosome 4   3.5%        CM024615.1, chromosome 7        9.3%        CM024613.1, chromosome 5   8.0%        CM024612.1, chromosome 4   6.8%
h1tg000056l  85102            CM024515.1, chromosome 2        78.6%       CM024516.1, chromosome 3   15.0%       NA                         0.0%        CM024610.1, chromosome 2        106.1%      CM024614.1, chromosome 6   1.7%        NA                         0.0%
h1tg000029l  84892            CM024516.1, chromosome 3        43.5%       CM024515.1, chromosome 2   30.5%       NA                         0.0%        CM024610.1, chromosome 2        79.4%       CM024614.1, chromosome 6   1.3%        NA                         0.0%
h1tg000064l  82462            CM024521.1, chromosome 8        50.2%       NA                         0.0%        NA                         0.0%        CM024616.1, chromosome 8        64.5%       NA                         0.0%        NA                         0.0%
h1tg000050l  81571            CM024524.1, chromosome 11       30.1%       NA                         0.0%        NA                         0.0%        CM024619.1, chromosome 11       36.7%       NA                         0.0%        NA                         0.0%
h1tg000057l  81096            CM024515.1, chromosome 2        59.5%       CM024516.1, chromosome 3   12.6%       NA                         0.0%        CM024610.1, chromosome 2        78.4%       CM024614.1, chromosome 6   6.8%        NA                         0.0%
h1tg000076l  80866            CM024515.1, chromosome 2        36.1%       CM024516.1, chromosome 3   14.9%       NA                         0.0%        CM024610.1, chromosome 2        48.1%       NA                         0.0%        NA                         0.0%
h1tg000034l  80746            CM024521.1, chromosome 8        16.6%       CM024522.1, chromosome 9   4.6%        CM024518.1, chromosome 5   3.7%        CM024616.1, chromosome 8        19.6%       CM024613.1, chromosome 5   6.9%        CM024614.1, chromosome 6   0.6%
h1tg000055l  80201            CM024525.1, chloroplast         76.9%       CM024514.1, chromosome 1   6.5%        CM024521.1, chromosome 8   1.9%        CM024616.1, chromosome 8        74.0%       CM024609.1, chromosome 1   6.5%        CM024617.1, chromosome 9   6.3%
h1tg000022l  79850            CM024516.1, chromosome 3        71.5%       NA                         0.0%        NA                         0.0%        CM024611.1, chromosome 3        42.2%       CM024615.1, chromosome 7   0.4%        NA                         0.0%
h1tg000036l  78549            CM024524.1, chromosome 11       59.3%       CM024520.1, chromosome 7   1.3%        NA                         0.0%        CM024619.1, chromosome 11       72.0%       CM024611.1, chromosome 3   5.0%        NA                         0.0%
h1tg000048l  77500            CM024518.1, chromosome 5        30.4%       CM024519.1, chromosome 6   0.1%        CM024520.1, chromosome 7   0.1%        CM024613.1, chromosome 5        69.5%       CM024611.1, chromosome 3   0.2%        CM024615.1, chromosome 7   0.1%
h1tg000084l  75001            CM024515.1, chromosome 2        91.9%       NA                         0.0%        NA                         0.0%        CM024610.1, chromosome 2        21.1%       CM024613.1, chromosome 5   0.4%        NA                         0.0%
h1tg000061l  73615            CM024521.1, chromosome 8        44.6%       CM024518.1, chromosome 5   0.6%        NA                         0.0%        CM024616.1, chromosome 8        84.7%       NA                         0.0%        NA                         0.0%
h1tg000033l  72644            CM024515.1, chromosome 2        71.1%       CM024516.1, chromosome 3   17.6%       NA                         0.0%        CM024610.1, chromosome 2        91.7%       CM024614.1, chromosome 6   1.9%        NA                         0.0%
h1tg000070l  71332            CM024515.1, chromosome 2        36.6%       NA                         0.0%        NA                         0.0%        CM024610.1, chromosome 2        50.7%       NA                         0.0%        NA                         0.0%
h1tg000073l  71031            CM024515.1, chromosome 2        34.2%       CM024520.1, chromosome 7   1.9%        CM024516.1, chromosome 3   0.3%        CM024610.1, chromosome 2        43.3%       CM024616.1, chromosome 8   8.2%        CM024613.1, chromosome 5   1.5%
h1tg000049l  68892            CM024519.1, chromosome 6        15.9%       CM024518.1, chromosome 5   11.8%       CM024520.1, chromosome 7   6.4%        CM024614.1, chromosome 6        30.9%       CM024612.1, chromosome 4   6.9%        CM024618.1, chromosome 10  0.8%
h1tg000039l  68546            CM024521.1, chromosome 8        29.2%       NA                         0.0%        NA                         0.0%        CM024616.1, chromosome 8        69.4%       NA                         0.0%        NA                         0.0%
h1tg000042l  68000            CM024521.1, chromosome 8        7.9%        CM024516.1, chromosome 3   4.8%        CM024517.1, chromosome 4   1.6%        CM024611.1, chromosome 3        7.7%        CM024612.1, chromosome 4   3.3%        CM024619.1, chromosome 11  1.3%
h1tg000079l  67339            CM024514.1, chromosome 1        47.7%       NA                         0.0%        NA                         0.0%        CM024609.1, chromosome 1        67.6%       NA                         0.0%        NA                         0.0%
h1tg000032l  67171            CM024520.1, chromosome 7        72.2%       NA                         0.0%        NA                         0.0%        CM024615.1, chromosome 7        53.7%       CM024613.1, chromosome 5   0.2%        CM024610.1, chromosome 2   0.2%
h1tg000072l  65110            CM024523.1, chromosome 10       6.0%        CM024518.1, chromosome 5   3.8%        CM024515.1, chromosome 2   3.4%        CM024610.1, chromosome 2        35.3%       CM024616.1, chromosome 8   2.9%        CM024615.1, chromosome 7   1.9%
h1tg000080l  64287            CM024518.1, chromosome 5        11.0%       CM024524.1, chromosome 11  4.8%        CM024514.1, chromosome 1   4.5%        CM024619.1, chromosome 11       90.7%       NA                         0.0%        NA                         0.0%
h1tg000020l  62283            CM024515.1, chromosome 2        69.7%       CM024516.1, chromosome 3   11.4%       NA                         0.0%        CM024610.1, chromosome 2        103.5%      CM024614.1, chromosome 6   0.8%        NA                         0.0%
h1tg000065l  61333            CM024515.1, chromosome 2        42.1%       CM024516.1, chromosome 3   39.5%       NA                         0.0%        CM024610.1, chromosome 2        82.1%       NA                         0.0%        NA                         0.0%
h1tg000078l  60163            CM024520.1, chromosome 7        4.2%        CM024516.1, chromosome 3   1.2%        CM024524.1, chromosome 11  1.0%        CM024619.1, chromosome 11       19.2%       NA                         0.0%        NA                         0.0%
h1tg000052l  59545            CM024521.1, chromosome 8        8.4%        CM024520.1, chromosome 7   8.4%        CM024524.1, chromosome 11  5.7%        CM024619.1, chromosome 11       85.4%       NA                         0.0%        NA                         0.0%
h1tg000068l  57909            CM024516.1, chromosome 3        55.0%       NA                         0.0%        NA                         0.0%        CM024611.1, chromosome 3        29.3%       CM024618.1, chromosome 10  1.4%        NA                         0.0%
h1tg000086l  57449            CM024514.1, chromosome 1        20.1%       CM024517.1, chromosome 4   3.3%        NA                         0.0%        CM024617.1, chromosome 9        7.5%        CM024609.1, chromosome 1   5.0%        CM024612.1, chromosome 4   2.5%
h1tg000037l  54248            CM024515.1, chromosome 2        63.7%       NA                         0.0%        NA                         0.0%        CM024610.1, chromosome 2        60.1%       NA                         0.0%        NA                         0.0%
h1tg000030l  53716            CM024514.1, chromosome 1        9.5%        CM024515.1, chromosome 2   8.1%        CM024517.1, chromosome 4   4.5%        CM024610.1, chromosome 2        17.7%       CM024618.1, chromosome 10  4.9%        CM024612.1, chromosome 4   4.3%
h1tg000083l  53008            CM024515.1, chromosome 2        41.2%       NA                         0.0%        NA                         0.0%        CM024610.1, chromosome 2        56.1%       NA                         0.0%        NA                         0.0%
h1tg000053l  50303            JABKBC010000014.1, tig00007996  6.3%        NA                         0.0%        NA                         0.0%        NA                              0.0%        NA                         0.0%        NA                         0.0%
h1tg000067l  49495            CM024524.1, chromosome 11       7.4%        CM024521.1, chromosome 8   7.1%        CM024522.1, chromosome 9   2.7%        CM024619.1, chromosome 11       37.5%       CM024615.1, chromosome 7   4.1%        CM024618.1, chromosome 10  2.5%
h1tg000062l  49067            CM024515.1, chromosome 2        47.8%       NA                         0.0%        NA                         0.0%        CM024610.1, chromosome 2        65.1%       NA                         0.0%        NA                         0.0%
h1tg000026l  48552            CM024525.1, chloroplast         73.4%       CM024516.1, chromosome 3   1.5%        CM024521.1, chromosome 8   1.0%        CM024611.1, chromosome 3        69.9%       CM024617.1, chromosome 9   8.4%        CM024616.1, chromosome 8   2.9%
h1tg000069l  46567            CM024515.1, chromosome 2        52.8%       CM024516.1, chromosome 3   29.9%       CM024522.1, chromosome 9   2.8%        CM024610.1, chromosome 2        62.5%       NA                         0.0%        NA                         0.0%
h1tg000074l  35174            CM024515.1, chromosome 2        43.4%       NA                         0.0%        NA                         0.0%        CM024610.1, chromosome 2        45.8%       CM024614.1, chromosome 6   1.5%        NA                         0.0%
h1tg000081l  33610            CM024515.1, chromosome 2        65.3%       CM024516.1, chromosome 3   14.4%       NA                         0.0%        CM024610.1, chromosome 2        91.2%       NA                         0.0%        NA                         0.0%
```


```
Your_Contig  Your_Contig_Len  vir_Ref_Chr1               %matchChr1  vir_Ref_Chr2               %matchChr2  vir_Ref_Chr3               %matchChr3  dec_Ref_Chr1                    %matchChr1  dec_Ref_Chr2               %matchChr2  dec_Ref_Chr3               %matchChr3
h2tg000005l  65688633         CM024521.1, chromosome 8   29.2%       CM024520.1, chromosome 7   1.8%        CM024517.1, chromosome 4   0.9%        CM024616.1, chromosome 8        52.7%       CM024614.1, chromosome 6   1.3%        CM024617.1, chromosome 9   0.7%
h2tg000001l  59206199         CM024518.1, chromosome 5   45.2%       CM024519.1, chromosome 6   1.9%        CM024520.1, chromosome 7   0.6%        CM024613.1, chromosome 5        25.6%       CM024611.1, chromosome 3   0.6%        CM024614.1, chromosome 6   0.4%
h2tg000008l  56293706         CM024516.1, chromosome 3   46.7%       CM024523.1, chromosome 10  3.1%        CM024518.1, chromosome 5   1.1%        CM024611.1, chromosome 3        28.8%       CM024614.1, chromosome 6   0.4%        CM024612.1, chromosome 4   0.2%
h2tg000011l  51100880         CM024519.1, chromosome 6   65.2%       CM024524.1, chromosome 11  1.0%        CM024515.1, chromosome 2   0.0%        CM024614.1, chromosome 6        35.9%       CM024615.1, chromosome 7   0.8%        CM024619.1, chromosome 11  0.5%
h2tg000004l  45062906         CM024520.1, chromosome 7   51.8%       CM024521.1, chromosome 8   0.5%        CM024524.1, chromosome 11  0.4%        CM024615.1, chromosome 7        30.0%       CM024617.1, chromosome 9   0.8%        CM024614.1, chromosome 6   0.8%
h2tg000012l  43652620         CM024514.1, chromosome 1   30.3%       CM024516.1, chromosome 3   1.5%        CM024518.1, chromosome 5   0.2%        CM024609.1, chromosome 1        56.7%       CM024619.1, chromosome 11  0.1%        CM024610.1, chromosome 2   0.1%
h2tg000006l  42555022         CM024524.1, chromosome 11  32.7%       CM024520.1, chromosome 7   1.1%        CM024521.1, chromosome 8   0.4%        CM024619.1, chromosome 11       63.7%       CM024617.1, chromosome 9   1.4%        CM024616.1, chromosome 8   0.1%
h2tg000009l  39808945         CM024522.1, chromosome 9   29.0%       CM024524.1, chromosome 11  3.4%        CM024518.1, chromosome 5   0.3%        CM024617.1, chromosome 9        55.6%       CM024612.1, chromosome 4   0.7%        CM024611.1, chromosome 3   0.2%
h2tg000003l  37845087         CM024517.1, chromosome 4   27.0%       CM024521.1, chromosome 8   2.3%        CM024523.1, chromosome 10  1.2%        CM024612.1, chromosome 4        55.9%       CM024615.1, chromosome 7   0.2%        CM024613.1, chromosome 5   0.2%
h2tg000007l  34610167         CM024523.1, chromosome 10  63.4%       CM024524.1, chromosome 11  1.7%        CM024520.1, chromosome 7   1.6%        CM024618.1, chromosome 10       36.1%       CM024612.1, chromosome 4   0.7%        CM024614.1, chromosome 6   0.5%
h2tg000010l  33511041         CM024515.1, chromosome 2   55.8%       CM024516.1, chromosome 3   2.6%        CM024519.1, chromosome 6   0.1%        CM024610.1, chromosome 2        34.8%       CM024609.1, chromosome 1   0.3%        CM024613.1, chromosome 5   0.3%
h2tg000002l  26036980         CM024519.1, chromosome 6   30.9%       CM024524.1, chromosome 11  0.8%        CM024518.1, chromosome 5   0.2%        CM024614.1, chromosome 6        57.8%       CM024615.1, chromosome 7   1.7%        CM024616.1, chromosome 8   0.0%
h2tg000014l  20522527         CM024515.1, chromosome 2   29.6%       CM024518.1, chromosome 5   0.3%        CM024516.1, chromosome 3   0.2%        CM024610.1, chromosome 2        47.9%       CM024614.1, chromosome 6   1.1%        CM024613.1, chromosome 5   0.4%
h2tg000013l  14252391         CM024518.1, chromosome 5   22.9%       CM024521.1, chromosome 8   0.5%        CM024517.1, chromosome 4   0.3%        CM024613.1, chromosome 5        33.9%       CM024611.1, chromosome 3   3.8%        CM024615.1, chromosome 7   0.5%
h2tg000015l  500252           CM024515.1, chromosome 2   69.7%       CM024516.1, chromosome 3   13.7%       CM024514.1, chromosome 1   0.3%        CM024610.1, chromosome 2        67.3%       CM024614.1, chromosome 6   0.5%        CM024615.1, chromosome 7   0.3%
h2tg000023l  280860           CM024515.1, chromosome 2   76.8%       CM024516.1, chromosome 3   8.1%        NA                         0.0%        CM024610.1, chromosome 2        62.5%       CM024614.1, chromosome 6   0.6%        NA                         0.0%
h2tg000040l  253708           CM024515.1, chromosome 2   42.6%       CM024516.1, chromosome 3   22.5%       NA                         0.0%        CM024610.1, chromosome 2        51.0%       CM024614.1, chromosome 6   2.9%        NA                         0.0%
h2tg000028l  227859           CM024515.1, chromosome 2   87.2%       NA                         0.0%        NA                         0.0%        CM024610.1, chromosome 2        83.5%       CM024614.1, chromosome 6   0.2%        NA                         0.0%
h2tg000018l  207414           CM024524.1, chromosome 11  18.6%       CM024514.1, chromosome 1   5.3%        CM024516.1, chromosome 3   5.2%        CM024617.1, chromosome 9        16.3%       CM024615.1, chromosome 7   9.7%        CM024612.1, chromosome 4   7.2%
h2tg000016l  193346           CM024525.1, chloroplast    74.5%       NA                         0.0%        NA                         0.0%        CM024611.1, chromosome 3        44.5%       CM024632.1, chloroplast    39.7%       NA                         0.0%
h2tg000027l  174770           CM024516.1, chromosome 3   63.4%       CM024515.1, chromosome 2   15.0%       NA                         0.0%        CM024610.1, chromosome 2        49.8%       NA                         0.0%        NA                         0.0%
h2tg000019l  164625           CM024515.1, chromosome 2   91.2%       CM024516.1, chromosome 3   3.5%        CM024514.1, chromosome 1   1.7%        CM024610.1, chromosome 2        99.2%       CM024615.1, chromosome 7   1.7%        CM024614.1, chromosome 6   1.3%
h2tg000017l  151689           CM024525.1, chloroplast    98.7%       NA                         0.0%        NA                         0.0%        CM024611.1, chromosome 3        95.8%       NA                         0.0%        NA                         0.0%
h2tg000031l  140617           CM024516.1, chromosome 3   14.1%       CM024521.1, chromosome 8   5.1%        CM024519.1, chromosome 6   3.6%        CM024609.1, chromosome 1        12.0%       CM024613.1, chromosome 5   10.2%       CM024611.1, chromosome 3   5.7%
h2tg000034l  118804           CM024517.1, chromosome 4   25.7%       CM024523.1, chromosome 10  13.2%       CM024516.1, chromosome 3   8.0%        CM024614.1, chromosome 6        49.7%       CM024613.1, chromosome 5   5.7%        CM024611.1, chromosome 3   2.7%
h2tg000045l  95361            CM024518.1, chromosome 5   3.6%        CM024520.1, chromosome 7   3.3%        CM024516.1, chromosome 3   2.8%        CM024610.1, chromosome 2        2.1%        CM024612.1, chromosome 4   1.7%        CM024618.1, chromosome 10  1.3%
h2tg000032l  93954            CM024515.1, chromosome 2   28.3%       CM024516.1, chromosome 3   22.3%       CM024517.1, chromosome 4   1.4%        CM024610.1, chromosome 2        32.7%       CM024614.1, chromosome 6   14.6%       NA                         0.0%
h2tg000021l  92347            CM024525.1, chloroplast    86.7%       CM024516.1, chromosome 3   4.9%        CM024517.1, chromosome 4   3.5%        CM024616.1, chromosome 8        93.8%       CM024613.1, chromosome 5   0.9%        CM024609.1, chromosome 1   0.6%
h2tg000026l  89563            CM024515.1, chromosome 2   55.8%       CM024516.1, chromosome 3   12.3%       NA                         0.0%        CM024610.1, chromosome 2        69.6%       CM024614.1, chromosome 6   12.6%       NA                         0.0%
h2tg000020l  87899            CM024518.1, chromosome 5   51.6%       NA                         0.0%        NA                         0.0%        CM024613.1, chromosome 5        79.0%       NA                         0.0%        NA                         0.0%
h2tg000038l  84936            CM024516.1, chromosome 3   36.1%       CM024521.1, chromosome 8   5.4%        CM024518.1, chromosome 5   1.3%        CM024611.1, chromosome 3        77.1%       CM024617.1, chromosome 9   0.4%        NA                         0.0%
h2tg000043l  75596            CM024520.1, chromosome 7   69.7%       CM024517.1, chromosome 4   0.2%        NA                         0.0%        CM024614.1, chromosome 6        60.7%       CM024609.1, chromosome 1   0.2%        CM024610.1, chromosome 2   0.2%
h2tg000030l  68747            CM024515.1, chromosome 2   63.3%       CM024522.1, chromosome 9   1.9%        CM024517.1, chromosome 4   1.6%        CM024610.1, chromosome 2        45.8%       NA                         0.0%        NA                         0.0%
h2tg000035l  66777            CM024522.1, chromosome 9   60.8%       CM024518.1, chromosome 5   0.3%        NA                         0.0%        CM024617.1, chromosome 9        26.2%       CM024619.1, chromosome 11  5.5%        CM024616.1, chromosome 8   2.0%
h2tg000033l  66304            CM024516.1, chromosome 3   54.4%       CM024514.1, chromosome 1   0.1%        NA                         0.0%        CM024611.1, chromosome 3        92.3%       NA                         0.0%        NA                         0.0%
h2tg000044l  66155            CM024524.1, chromosome 11  93.1%       NA                         0.0%        NA                         0.0%        CM024611.1, chromosome 3        23.9%       CM024619.1, chromosome 11  11.4%       NA                         0.0%
h2tg000037l  63729            CM024516.1, chromosome 3   23.9%       CM024519.1, chromosome 6   7.0%        CM024514.1, chromosome 1   1.2%        CM024611.1, chromosome 3        22.4%       CM024616.1, chromosome 8   6.5%        CM024615.1, chromosome 7   1.0%
h2tg000042l  63218            CM024523.1, chromosome 10  35.7%       CM024514.1, chromosome 1   4.8%        CM024521.1, chromosome 8   0.8%        CM024618.1, chromosome 10       45.8%       CM024615.1, chromosome 7   5.9%        CM024610.1, chromosome 2   0.5%
h2tg000024l  63097            CM024520.1, chromosome 7   51.3%       NA                         0.0%        NA                         0.0%        CM024614.1, chromosome 6        52.1%       CM024611.1, chromosome 3   0.3%        NA                         0.0%
h2tg000041l  62567            CM024515.1, chromosome 2   39.5%       CM024516.1, chromosome 3   13.1%       CM024514.1, chromosome 1   4.0%        CM024610.1, chromosome 2        57.9%       CM024615.1, chromosome 7   4.1%        CM024614.1, chromosome 6   3.4%
h2tg000039l  57976            CM024517.1, chromosome 4   40.7%       CM024524.1, chromosome 11  0.3%        NA                         0.0%        CM024612.1, chromosome 4        63.0%       NA                         0.0%        NA                         0.0%
h2tg000025l  56855            CM024524.1, chromosome 11  92.7%       NA                         0.0%        NA                         0.0%        CM024619.1, chromosome 11       39.2%       CM024614.1, chromosome 6   2.9%        CM024615.1, chromosome 7   0.7%
h2tg000022l  49225            CM024518.1, chromosome 5   26.3%       CM024514.1, chromosome 1   4.9%        CM024519.1, chromosome 6   3.0%        JABKBM010000015.1, tig00006793  62.5%       CM024613.1, chromosome 5   62.3%       NA                         0.0%
h2tg000029l  45841            CM024521.1, chromosome 8   37.3%       CM024516.1, chromosome 3   0.1%        CM024518.1, chromosome 5   0.1%        CM024616.1, chromosome 8        39.5%       CM024615.1, chromosome 7   3.9%        CM024612.1, chromosome 4   0.2%
h2tg000036l  42116            CM024515.1, chromosome 2   58.4%       CM024516.1, chromosome 3   33.9%       CM024522.1, chromosome 9   3.0%        CM024610.1, chromosome 2        67.5%       NA                         0.0%        NA                         0.0%
```

We can summarise this succinctly by, for each contig, figuring out which parent has better coverage, and what the likely chromosome is.

Here's a script to do that.

```bash
BASE_DIR="06_1_hifiasm_assemblies/l3/genome_alignment"
OUTPUT="$BASE_DIR/final_assembly_summary.tsv"

echo -e "Contig\tContig_Len\tPrimary_Chr\tMinimap\tVirginea_Cov%\tDecipiens_Cov%" > "$OUTPUT"

for HAP_FILE in "$BASE_DIR/hap1_homology_percentages.tsv" "$BASE_DIR/hap2_homology_percentages.tsv"; do
    if [ ! -f "$HAP_FILE" ]; then continue; fi

    awk -F'\t' 'NR > 1 {
        qlen = $2
        vir_name = $3
        dec_name = $9
        
        # Clean percentages
        v_pct_raw = $4; sub(/%/, "", v_pct_raw); v_pct = (v_pct_raw == "NA" ? 0 : v_pct_raw)
        d_pct_raw = $10; sub(/%/, "", d_pct_raw); d_pct = (d_pct_raw == "NA" ? 0 : d_pct_raw)

        # 1. Determine Winner
        if (v_pct >= d_pct && v_pct > 0) {
            winner = "virginea"
            raw_full_name = vir_name
        } else if (d_pct > v_pct) {
            winner = "decipiens"
            raw_full_name = dec_name
        } else {
            winner = "unassigned"
            raw_full_name = "unmapped"
        }

        # 2. Extract Name using Regex (Look for chromosome, chloroplast, or tig)
        # This looks for the keywords and captures everything after them
        clean_name = raw_full_name
        if (raw_full_name ~ /chromosome/) {
            match(raw_full_name, /chromosome [0-9]+/)
            clean_name = substr(raw_full_name, RSTART, RLENGTH)
        } else if (raw_full_name ~ /chloroplast/) {
            clean_name = "chloroplast"
        } else if (raw_full_name ~ /tig/) {
            match(raw_full_name, /tig[0-9]+/)
            clean_name = substr(raw_full_name, RSTART, RLENGTH)
        } else if (raw_full_name ~ /, /) {
            # Fallback: if keywords fail but comma exists, take part after comma
            split(raw_full_name, parts, ", ")
            clean_name = parts[2]
        }

        printf "%s\t%s\t%s\t%s\t%.1f%%\t%.1f%%\n", $1, $2, clean_name, winner, v_pct, d_pct
    }' "$HAP_FILE" >> "$OUTPUT"
done

# Sort by length
(head -n 1 "$OUTPUT" && tail -n +2 "$OUTPUT" | sort -t$'\t' -k2,2rn) > "${OUTPUT}.tmp" && mv "${OUTPUT}.tmp" "$OUTPUT"

echo "Fixed! Names should now appear correctly for all rows."
```

That gives the following table:

```
awk -F'\t' '
    BEGIN { print "| " }
    { 
        printf "|"; 
        for(i=1; i<=NF; i++) printf " %s |", $i; 
        print ""; 
        if(NR==1) {
            printf "|"; 
            for(i=1; i<=NF; i++) printf "---|"; 
            print ""
        }
    }' 06_1_hifiasm_assemblies/l3/genome_alignment/final_assembly_summary.tsv
```


The key is that for all of the largest contigs, the coverage gives a very clear winner. A few chromosomes are broken, but these mostly add up to the right number.

For example, Chr5 in decipiens is in one piece from virginea `h2tg000001l` at 59.2Mbp. The decipiens copy comprises two major pieces `h1tg000014l` is 45.2Mbp, and `h2tg000013l` at 14.2Mbp, which sum to 59.4Mbp, very close to the length of the virginea copy. 

Chr2 is certainly the most fragmented, and there are certain contigs where it's essentially impossible to pick a winning parent. E.g. `h2tg000028l` is 227859bp and matches virginea by 87.2% and decipiens by 83.5.

Despite this, one could still bin them and scaffold them so as not to lose the information. 


| Contig | Contig_Len | Primary_Chr | Minimap | Virginea_Cov% | Decipiens_Cov% |                                                                                                                                                                                                          
|---|---|---|---|---|---|                                                                                                                                                                                                                                                                 
| h2tg000005l | 65688633 | chromosome 8 | decipiens | 29.2% | 52.7% |                                                                                                                                                                                                                     
| h1tg000010l | 61716932 | chromosome 8 | virginea | 52.9% | 31.5% |                                                                                                                                                                                                                      
| h1tg000006l | 60158556 | chromosome 3 | decipiens | 24.5% | 45.7% |                                                                                                                                                                                                                     
| h2tg000001l | 59206199 | chromosome 5 | virginea | 45.2% | 25.6% |                                                                                                                                                                                                                      
| h2tg000008l | 56293706 | chromosome 3 | virginea | 46.7% | 28.8% |                                                                                                                                                                                                                      
| h1tg000002l | 54083987 | chromosome 6 | decipiens | 34.0% | 60.2% |                                                                                                                                                                                                                     
| h1tg000001l | 51333920 | chromosome 7 | decipiens | 27.8% | 47.7% |                                                                                                                                                                                                                     
| h2tg000011l | 51100880 | chromosome 6 | virginea | 65.2% | 35.9% |                                                                                                                                                                                                                      
| h1tg000014l | 45216074 | chromosome 5 | decipiens | 25.8% | 46.7% |                                                                                                                                                                                                                     
| h2tg000004l | 45062906 | chromosome 7 | virginea | 51.8% | 30.0% |                                                                                                                                                                                                                      
| h2tg000012l | 43652620 | chromosome 1 | decipiens | 30.3% | 56.7% |                                                                                                                                                                                                                     
| h2tg000006l | 42555022 | chromosome 11 | decipiens | 32.7% | 63.7% |                                                                                                                                                                                                                    
| h1tg000007l | 41950308 | chromosome 1 | virginea | 60.7% | 33.8% |                                                                                                                                                                                                                      
| h1tg000008l | 39810478 | chromosome 11 | virginea | 65.0% | 37.0% |                                                                                                                                                                                                                     
| h2tg000009l | 39808945 | chromosome 9 | decipiens | 29.0% | 55.6% |                                                                                                                                                                                                                     
| h1tg000003l | 38464560 | chromosome 10 | decipiens | 32.0% | 57.2% |                                                                                                                                                                                                                    
| h2tg000003l | 37845087 | chromosome 4 | decipiens | 27.0% | 55.9% |                                                                                                                                                                                                                     
| h1tg000005l | 36525698 | chromosome 9 | virginea | 56.4% | 34.9% |                                                                                                                                                                                                                      
| h1tg000011l | 35389336 | chromosome 4 | virginea | 54.2% | 34.9% |                                                                                                                                                                                                                      
| h2tg000007l | 34610167 | chromosome 10 | virginea | 63.4% | 36.1% |                                                                                                                                                                                                                     
| h1tg000009l | 33914441 | chromosome 2 | decipiens | 32.7% | 57.1% |                                                                                                                                                                                                                     
| h2tg000010l | 33511041 | chromosome 2 | virginea | 55.8% | 34.8% |                                                                                                                                                                                                                      
| h2tg000002l | 26036980 | chromosome 6 | decipiens | 30.9% | 57.8% |                                                                                                                                                                                                                     
| h2tg000014l | 20522527 | chromosome 2 | decipiens | 29.6% | 47.9% |                                                                                                                                                                                                                     
| h1tg000004l | 16945854 | chromosome 2 | virginea | 55.5% | 29.5% |                                                                                                                                                                                                                      
| h2tg000013l | 14252391 | chromosome 5 | decipiens | 22.9% | 33.9% |                                                                                                                                                                                                                     
| h1tg000021l | 550977 | chromosome 2 | virginea | 79.9% | 64.3% |                                                                                                                                                                                                                        
| h1tg000019l | 550701 | chromosome 3 | decipiens | 24.6% | 26.4% |                                                                                                                                                                                                                       
| h1tg000027l | 504488 | chromosome 9 | decipiens | 1.3% | 2.4% |                                                                                                                                                                                                                         
| h2tg000015l | 500252 | chromosome 2 | virginea | 69.7% | 67.3% |                                                                                                                                                                                                                        
| h2tg000023l | 280860 | chromosome 2 | virginea | 76.8% | 62.5% |                                                                                                                                                                                                                        
| h1tg000018l | 267638 | chromosome 11 | virginea | 13.5% | 13.2% |                                                                                                                                                                                                                       
| h1tg000031l | 262582 | chromosome 11 | decipiens | 10.5% | 90.8% |                                                                                                                                                                                                                      
| h2tg000040l | 253708 | chromosome 2 | decipiens | 42.6% | 51.0% |                                                                                                                                                                                                                       
| h1tg000017l | 247779 | chromosome 2 | virginea | 74.3% | 52.8% |                                                                                                                                                                                                                        
| h1tg000044l | 246068 | tig00112428 | decipiens | 1.1% | 3.4% |                                                                                                                                                                                                                          
| h1tg000015l | 240840 | chromosome 2 | virginea | 74.1% | 20.3% |                                                                                                                                                                                                                        
| h2tg000028l | 227859 | chromosome 2 | virginea | 87.2% | 83.5% |                                                                                                                                                                                                                        
| h2tg000018l | 207414 | chromosome 11 | virginea | 18.6% | 16.3% |                                                                                                                                                                                                                       
| h1tg000045l | 202662 | chromosome 2 | decipiens | 85.8% | 99.8% |                                                                                                                                                                                                                       
| h2tg000016l | 193346 | chloroplast | virginea | 74.5% | 44.5% |                                                                                                                                                                                                                         
| h1tg000023l | 191335 | chloroplast | virginea | 86.5% | 83.8% |                                                                                                                                                                                                                         
| h1tg000012l | 179571 | chromosome 2 | virginea | 82.2% | 60.7% |
| h2tg000027l | 174770 | chromosome 3 | virginea | 63.4% | 49.8% |
| h1tg000047l | 173338 | chromosome 2 | virginea | 88.1% | 84.5% |
| h1tg000040l | 164941 | chromosome 2 | decipiens | 59.6% | 80.7% |
| h2tg000019l | 164625 | chromosome 2 | decipiens | 91.2% | 99.2% |
| h1tg000066l | 162610 | tig00006793 | decipiens | 15.5% | 52.6% |
| h1tg000046l | 154476 | chromosome 3 | virginea | 12.6% | 12.6% |
| h2tg000017l | 151689 | chloroplast | virginea | 98.7% | 95.8% |
| h2tg000031l | 140617 | chromosome 3 | virginea | 14.1% | 12.0% |
| h1tg000060l | 128970 | tig00112428 | decipiens | 0.4% | 3.9% |
| h1tg000059l | 120315 | chromosome 7 | virginea | 9.8% | 9.3% |
| h1tg000071l | 118804 | chromosome 6 | decipiens | 25.7% | 49.7% |
| h2tg000034l | 118804 | chromosome 6 | decipiens | 25.7% | 49.7% |
| h1tg000024l | 110453 | chromosome 2 | virginea | 82.3% | 108.5% |
| h1tg000058l | 108016 | chromosome 2 | decipiens | 53.6% | 57.3% |
| h1tg000035l | 106259 | chromosome 3 | virginea | 63.4% | 39.3% |
| h2tg000045l | 95361 | chromosome 5 | virginea | 3.6% | 2.1% |
| h1tg000025l | 94874 | chromosome 2 | decipiens | 60.2% | 87.9% |
| h1tg000043l | 94628 | chromosome 3 | decipiens | 87.0% | 89.1% |
| h1tg000085l | 94551 | chromosome 1 | decipiens | 15.9% | 73.3% |
| h2tg000032l | 93954 | chromosome 2 | decipiens | 28.3% | 32.7% |
| h1tg000016l | 93941 | chromosome 1 | virginea | 9.5% | 18.8% |
| h2tg000021l | 92347 | chromosome 8 | decipiens | 86.7% | 93.8% |
| h2tg000026l | 89563 | chromosome 2 | decipiens | 55.8% | 69.6% |
| h1tg000028l | 88856 | chromosome 1 | decipiens | 40.6% | 61.8% |
| h2tg000020l | 87899 | chromosome 5 | decipiens | 51.6% | 79.0% |
| h1tg000063l | 87206 | chromosome 1 | decipiens | 32.9% | 63.5% |
| h1tg000038l | 87121 | chromosome 7 | virginea | 9.8% | 9.3% |
| h1tg000056l | 85102 | chromosome 2 | virginea | 78.6% | 106.1% |
| h2tg000038l | 84936 | chromosome 3 | decipiens | 36.1% | 77.1% |
| h1tg000029l | 84892 | chromosome 2 | decipiens | 43.5% | 79.4% |
| h1tg000064l | 82462 | chromosome 8 | decipiens | 50.2% | 64.5% |
| h1tg000050l | 81571 | chromosome 11 | decipiens | 30.1% | 36.7% |
| h1tg000057l | 81096 | chromosome 2 | decipiens | 59.5% | 78.4% |
| h1tg000076l | 80866 | chromosome 2 | decipiens | 36.1% | 48.1% |
| h1tg000034l | 80746 | chromosome 8 | decipiens | 16.6% | 19.6% |
| h1tg000055l | 80201 | chloroplast | virginea | 76.9% | 74.0% |
| h1tg000022l | 79850 | chromosome 3 | virginea | 71.5% | 42.2% |
| h1tg000036l | 78549 | chromosome 11 | decipiens | 59.3% | 72.0% |
| h1tg000048l | 77500 | chromosome 5 | decipiens | 30.4% | 69.5% |
| h2tg000043l | 75596 | chromosome 7 | virginea | 69.7% | 60.7% |
| h1tg000084l | 75001 | chromosome 2 | virginea | 91.9% | 21.1% |
| h1tg000061l | 73615 | chromosome 8 | decipiens | 44.6% | 84.7% |
| h1tg000033l | 72644 | chromosome 2 | decipiens | 71.1% | 91.7% |
| h1tg000070l | 71332 | chromosome 2 | decipiens | 36.6% | 50.7% |
| h1tg000073l | 71031 | chromosome 2 | decipiens | 34.2% | 43.3% |
| h1tg000049l | 68892 | chromosome 6 | decipiens | 15.9% | 30.9% |
| h2tg000030l | 68747 | chromosome 2 | virginea | 63.3% | 45.8% |
| h1tg000039l | 68546 | chromosome 8 | decipiens | 29.2% | 69.4% |
| h1tg000042l | 68000 | chromosome 8 | virginea | 7.9% | 7.7% |
| h1tg000079l | 67339 | chromosome 1 | decipiens | 47.7% | 67.6% |
| h1tg000032l | 67171 | chromosome 7 | virginea | 72.2% | 53.7% |
| h2tg000035l | 66777 | chromosome 9 | virginea | 60.8% | 26.2% |
| h2tg000033l | 66304 | chromosome 3 | decipiens | 54.4% | 92.3% |
| h2tg000044l | 66155 | chromosome 11 | virginea | 93.1% | 23.9% |
| h1tg000072l | 65110 | chromosome 10 | virginea | 6.0% | 35.3% |
| h1tg000080l | 64287 | chromosome 11 | decipiens | 11.0% | 90.7% |
| h2tg000037l | 63729 | chromosome 3 | virginea | 23.9% | 22.4% |
| h2tg000042l | 63218 | chromosome 10 | decipiens | 35.7% | 45.8% |
| h2tg000024l | 63097 | chromosome 6 | decipiens | 51.3% | 52.1% |
| h2tg000041l | 62567 | chromosome 2 | decipiens | 39.5% | 57.9% |
| h1tg000020l | 62283 | chromosome 2 | virginea | 69.7% | 103.5% |
| h1tg000065l | 61333 | chromosome 2 | decipiens | 42.1% | 82.1% |
| h1tg000078l | 60163 | chromosome 7 | virginea | 4.2% | 19.2% |
| h1tg000052l | 59545 | chromosome 11 | decipiens | 8.4% | 85.4% |
| h2tg000039l | 57976 | chromosome 4 | decipiens | 40.7% | 63.0% |
| h1tg000068l | 57909 | chromosome 3 | virginea | 55.0% | 29.3% |
| h1tg000086l | 57449 | chromosome 9 | decipiens | 20.1% | 7.5% |
| h2tg000025l | 56855 | chromosome 11 | virginea | 92.7% | 39.2% |
| h1tg000037l | 54248 | chromosome 2 | virginea | 63.7% | 60.1% |
| h1tg000030l | 53716 | chromosome 1 | virginea | 9.5% | 17.7% |
| h1tg000083l | 53008 | chromosome 2 | decipiens | 41.2% | 56.1% |
| h1tg000053l | 50303 | tig00007996 | virginea | 6.3% | 0.0% |
| h1tg000067l | 49495 | chromosome 11 | virginea | 7.4% | 37.5% |
| h2tg000022l | 49225 | tig00006793 | decipiens | 26.3% | 62.5% |
| h1tg000062l | 49067 | chromosome 2 | decipiens | 47.8% | 65.1% |
| h1tg000026l | 48552 | chloroplast | virginea | 73.4% | 69.9% |
| h1tg000069l | 46567 | chromosome 2 | decipiens | 52.8% | 62.5% |
| h2tg000029l | 45841 | chromosome 8 | decipiens | 37.3% | 39.5% |
| h2tg000036l | 42116 | chromosome 2 | decipiens | 58.4% | 67.5% |
| h1tg000074l | 35174 | chromosome 2 | decipiens | 43.4% | 45.8% |
| h1tg000081l | 33610 | chromosome 2 | decipiens | 65.3% | 91.2% |

## 2 Parental kmer checks

Similar to what I did on the initial assembly, now I want to see how the unique kmers from the putative parental genomes map to the contigs. I already have hapmers from the parents, with k=31. this is described above, but I get the kmers, look for the noise | signal trough, filter out the noise, and then use subtraction to get the kmers unique to each parent (hapmers). Now I need to do the same for the phylacis reads, and then I can run merqury.


```bash
mkdir 06_2_parental_kmer_checks
cd 06_2_parental_kmer_checks

# make phylacis meryl db (not sure how this will go with nanopore, but let's see)
meryl count k=31 threads=128 memory=1200G \
    output E_phylacis.meryl \
    ../02_filtering/E_phylacis_filtered.fastq.gz

meryl histogram E_phylacis.meryl > E_phylacis.hist

# look for the trough
E_phylacis.hist | more

# Filter based on trough - the number is the trough I observed
meryl greater-than 10 E_phylacis.meryl output E_phylacis.filtered.meryl

# Run merqury twice - once with filtered and once with unfilitered read kmers
cd ~/MM_genome/06_1_hifiasm_assemblies/l3/
gfatools gfa2fa E_phylacis.bp.hap1.p_ctg.gfa > E_phylacis.hap1.fasta
gfatools gfa2fa E_phylacis.bp.hap2.p_ctg.gfa > E_phylacis.hap2.fasta

# 3. The two commands above lead merqury to fail. So next I'll try running it as I ran it initially on the top 11 contigs
mkdir trio
cd trio

# merqury is very finicky about having the actual files available, and fails a lot when it tries to make symlinks if they're not right there, so we copy them over for this analysis
cp -r  ../E_phylacis.meryl/ . # phylacis kmers from nanopore reads, no filtering for noise
cp -r ~/MM_genome/parental_spp_genomes/E_decipiens.final_probes.meryl/ .  # decipiens hapmers, noise filtered
cp -r ~/MM_genome/parental_spp_genomes/E_virginea.final_probes.meryl/ .   # virginea hapmers, noise filtered
cp ../../06_1_hifiasm_assemblies/l3/E_phylacis.hap1.fasta . # phylacis assembly hap1
cp ../../06_1_hifiasm_assemblies/l3/E_phylacis.hap2.fasta . # phylacis assembly hap2

merqury.sh \
    E_phylacis.meryl \
    E_virginea.final_probes.meryl \
    E_decipiens.final_probes.meryl \
    E_phylacis.hap1.fasta \
    E_phylacis.hap2.fasta \
    E_phylacis_trio

```

## 3 ROADIES

Get the genomes

```bash

mkdir 06_3_roadies
cd 06_3_roadies

# 1a. Download the packages (Corymbia is taxid 87658, needed because it's also a beetle genus)
datasets download genome taxon "Eucalyptus" --reference --assembly-level chromosome --filename eucs.zip
datasets download genome taxon "87658" --assembly-level chromosome,scaffold --filename corymbia.zip
datasets download genome taxon "Angophora" --assembly-level chromosome,scaffold --filename angophora.zip

# 1b. Create the reference directory
mkdir -p roadies_refs

# 1c. Unzip only the genomic fasta files (.fna) directly into the folder
# The -j (junk paths) flag is key here to avoid deep folder structures
unzip -j eucs.zip "ncbi_dataset/data/*/*.fna" -d roadies_refs/
unzip -j corymbia.zip "ncbi_dataset/data/*/*.fna" -d roadies_refs/
unzip -j angophora.zip "ncbi_dataset/data/*/*.fna" -d roadies_refs/

# Get metadata so we can rename files by species
cd roadies_refs
datasets summary genome taxon "Eucalyptus" --reference --assembly-level chromosome --as-json-lines > eucs_metadata.jsonl
datasets summary genome taxon "Angophora" --reference --assembly-level chromosome --as-json-lines > angophora_metadata.jsonl
datasets summary genome taxon "87658" --reference --assembly-level chromosome --as-json-lines > corymbia_metadata.jsonl

# finally, we rename the genome files like this:
python3 -c "
import os, json

# 1. Build a dictionary of Accession -> Species Name from all metadata files
mapping = {}
for meta_file in [f for f in os.listdir('.') if f.endswith('.jsonl')]:
    with open(meta_file) as f:
        for line in f:
            data = json.loads(line)
            # Access the nested accession and scientific name
            acc = data['accession']
            species = data['organism']['organism_name'].replace(' ', '_')
            mapping[acc] = species

# 2. Iterate through the genome files and rename them
for filename in os.listdir('.'):
    if filename.endswith('.fna'):
        # Check if the accession (GCA/GCF_...) is in the filename
        for acc, species in mapping.items():
            if acc in filename:
                new_name = f'{species}.fna'
                print(f'Renaming {filename} to {new_name}')
                os.rename(filename, new_name)
                break
"


```

That's our reference data sorted. Now we want to treat each of the contigs as its own 'species'. First let's see how long the contigs are...

```
seqkit stats ../../06_1_hifiasm_assemblies/l3/E_phylacis.hap*.fasta
```

```
file                                                    format  type  num_seqs      sum_len  min_len      avg_len     max_len
../../06_1_hifiasm_assemblies/l3/E_phylacis.hap1.fasta  FASTA   DNA         86  523,986,726   30,917  6,092,868.9  61,716,932
../../06_1_hifiasm_assemblies/l3/E_phylacis.hap2.fasta  FASTA   DNA         45  573,953,311   42,116   12,754,518  65,688,633
```

The plan now is to treat each of these as a species, and see where they land w.r.t. the other genomes.

After installing ROADIES from source, we'll start by making an input folder with the reference genomes from above, renamed for roadies


```
cd ~/MM_genome/ROADIES
mkdir -p euc_refs
cp ../06_3_roadies/roadies_refs/*.fna euc_refs/

cd euc_refs

for file in *.fna; do
    # Replace 'subsp._' with 'subsp_' to remove the extra dot
    new_name=$(echo "$file" | sed 's/subsp\._/subsp_/g')
    # Change extension from .fna to .fa
    new_name="${new_name%.fna}.fa"    
    mv "$file" "$new_name"
done

cd ..
```

OK, now let's make a file with all of our haplotype contigs as individual fasta files

```bash
mkdir -p query_contigs
cp ../06_3_roadies/roadies_input/h[12]tg*.fa query_contigs/
```

NExt, I modified the ROADIES code to allow me to sample from a focal species, instead of from all species. This means that I can specify a single one of hte input genomes, and all of the loci will be sampled from that genome. This is what I need to allow me to try and 'place' individual contigs. This is particulary the case when an individual contig is <<1 chromosome. In other words, what ROADIES will now do is sample N genes from the focal species, and look for homologues on the other genomes.

This invovles edits to just two files. To run this yourself, you'll need to install the latest version of ROADIES then update these two files on this repo:

```bash
git add -f ROADIES/config/config.yaml
git add -f ROADIES/workflow/rules/sampling.smk
```

Now we want to run roadies for each of the contigs that we have for our assembly...

Run the following from the ROADIES/ directory. My initial CONFIG is set up to:

* sample 1000 loci of 500bp each (if there's enough in the focal contig to do so)
* minimum 95% uppercase in the alignments
* minimum identity 85%
* coverage and continuity also have to be at last 85%
* max duplications 4
* 64 instances

The rest are just the defaults. This is set for a fairly conservative approach, looking for quite similar seuqences among the genomes.

In this part I just get all the contigs, and run ROADIES sequentially, once for each contig. Running ROADIES like this is very I/O heavy, so to keep things running quickly and to save my disks, I do it all on a big ramdisk, and just copy the output back to where it needs to go. That part is in the loop.

```bash
# --- PREP RAMDISK ---
sudo mkdir -p /mnt/ramdisk
sudo mount -t tmpfs -o size=1000G tmpfs /mnt/ramdisk
sudo chown rob:rob /mnt/ramdisk

# --- DEFINE ABSOLUTE PATHS (The "Safe" Way) ---
# Use $(pwd) or full paths so these work from anywhere
BASE_DIR="/home/rob/MM_genome/ROADIES"
REF_DIR="$BASE_DIR/euc_refs"
QUERY_DIR="$BASE_DIR/query_contigs"
FINAL_DEST="$BASE_DIR/euc_contigs_ramdisk"

mkdir -p "$FINAL_DEST"

RAM_DIR="/mnt/ramdisk/roadies_work"
RAM_ROADIES="/mnt/ramdisk/ROADIES_CODE"

# 1. Copy ROADIES code to RAM
echo "Migrating ROADIES code to RAM..."
cp -r "$BASE_DIR" "$RAM_ROADIES"

# 2. MOVE into the RAM directory to execute
cd "$RAM_ROADIES"


# 3. Get the list of all contigs
CONTIG_LIST=$(ls "$QUERY_DIR"/*.fa | xargs -n 1 basename)

# --- PILOT OVERWRITE ---
CONTIG_LIST="h1tg000010l.fa h2tg000020l.fa"
echo "Running PILOT mode with: $CONTIG_LIST"

# 4. The Loop
for CONTIG in $CONTIG_LIST; do
    ID="${CONTIG%.fa}"
    echo "Processing $ID..."
    
    # Clean and set up the RAM workspace
    rm -rf "$RAM_DIR"
    mkdir -p "$RAM_DIR/input" "$RAM_DIR/output"

    # Copy data TO ramdisk
    cp "$REF_DIR"/*.fa "$RAM_DIR/input/"
    cp "$QUERY_DIR/$CONTIG" "$RAM_DIR/input/"   

    #Calculate target loci to avoid overlap (1 per 100,000bp)
    # initial testing suggests that 1 500bp locus per roughly 100Kb is good, and mostly avoids overlap in the *accepted* loci well)
    CONTIG_PATH="$QUERY_DIR/$CONTIG"
    CONTIG_LEN=$(awk '/^>/ {next} {len += length($0)} END {print len}' "$CONTIG_PATH")
    N_LOCI=$(( CONTIG_LEN / 100000 ))

    # Apply constraints: Min 25, Max 1000
    if [ "$N_LOCI" -lt 25 ]; then
        N_LOCI=25
    elif [ "$N_LOCI" -gt 100 ]; then
        N_LOCI=100
    fi

    echo "------------------------------------------------"
    echo "Processing: $ID"
    echo "Contig length: $CONTIG_LEN bp"
    echo "Targeting 50% coverage -> N_LOCI: $N_LOCI"
    echo "------------------------------------------------"

    # 3. Update the config with the dynamic N_LOCI
    sed -e "s|^GENOMES:.*|GENOMES: \"$RAM_DIR/input\"|" \
        -e "s|^FOCAL_GENOME:.*|FOCAL_GENOME: \"$CONTIG\"|" \
        -e "s|^OUT_DIR:.*|OUT_DIR: \"$RAM_DIR/output\"|" \
        -e "s|^GENE_COUNT:.*|GENE_COUNT: $N_LOCI|" \
        config/config.yaml > "$RAM_DIR/temp_config.yaml"

    # Run ROADIES (Now running the RAM version of python/snakemake)
    cp "$RAM_DIR/temp_config.yaml" config/config.yaml
    python run_roadies.py --noconverge --cores 256 --config "config/config.yaml"

    # HARVEST results
    echo "Harvesting key results for $ID..."
    CONTIG_DEST="$FINAL_DEST/$ID"
    mkdir -p "$CONTIG_DEST"

    # 1. Copy the trees
    cp "$RAM_DIR/output/roadies.nwk" "$CONTIG_DEST/" 2>/dev/null
    cp "$RAM_DIR/output/roadies_stats.nwk" "$CONTIG_DEST/" 2>/dev/null

    # 2. Copy the plots and statistics folders
    cp -r "$RAM_DIR/output/plots" "$CONTIG_DEST/"
    cp -r "$RAM_DIR/output/statistics" "$CONTIG_DEST/"

    # 3. Copy the timestamp and config for record-keeping
    cp "$RAM_DIR/output/time_stamps.csv" "$CONTIG_DEST/"
    cp "$RAM_DIR/temp_config.yaml" "$CONTIG_DEST/"

    cp "$RAM_DIR/temp_config.yaml" "$CONTIG_DEST/"
    
    # Optional: Wipe RAM work dir but keep RAM_ROADIES for next iteration
    rm -rf "$RAM_DIR"
done
```



### ROADIES results

One way to summarise the results is to just look at where the haplotigs land in the species tree. 

Gemini wrote me a little script for that. If you point it at the folder that has the results summaries, you can run it like this:

(the input folder is hardcoded at the top of the script, so you need to change it to suit wherever the roadies output summaries are)

```bash
python summarise_haplotigs.py
```

The results are absolutely completely amazing, and beyond my wildest expectations here...

The haplotigs are mostly landing with pp of 1.0 next to *either* virginea or decipiens, but no other taxa.

Just as remarkable, Q1 is usually >70%, and Q2 and Q3 are very close to equal, suggesting that incongruence is roughly attributable to ILS only.

I'll update these as the results come in...

```
Haplotig        Sister_Group         Genes  BLen       pp1    f1       f2       f3       q1       q2       q3      
-------------------------------------------------------------------------------------------------------------------
h1tg000001l     Eucalyptus_decipiens 880    0.000065   1.000000 199.85   22.352   25.764   0.805954 0.090144 0.103902
h1tg000002l     Eucalyptus_decipiens 899    0.000514   1.000000 318.55   24.294   27.117   0.861038 0.065665 0.073297
h1tg000003l     Eucalyptus_decipiens 883    0.000034   1.000000 316.61   26.294   28.147   0.853282 0.070862 0.075856
h1tg000004l     Eucalyptus_virginea  908    0.000433   1.000000 224.88   36.617   36.117   0.755608 0.123036 0.121356
h1tg000005l     Eucalyptus_virginea  890    0.000371   1.000000 288.67   42.823   40.558   0.775889 0.115099 0.109012
h1tg000006l     Eucalyptus_decipiens 878    0.000261   1.000000 164.08   33.088   28.852   0.725960 0.146389 0.127651
h1tg000007l     Eucalyptus_virginea  895    0.000953   1.000000 287.23   38.235   34.882   0.797094 0.106105 0.096801
h1tg000008l     Eucalyptus_virginea  906    0.001965   1.000000 276.17   40.352   39.882   0.774880 0.113220 0.111900
h1tg000009l     Eucalyptus_decipiens 893    0.001254   1.000000 268.55   34.323   26.000   0.816580 0.104364 0.079056
h1tg000010l     Eucalyptus_virginea  909    0.000484   1.000000 247.23   38.323   38.705   0.762449 0.118186 0.119365
```



## Parental assignment table

Finally, we can put all of the above together to make a table which allows us to assign contigs to the father or mother with some confidence. 

Everything above suggests that decipiens and virginea are the likely parents, and that virginea is the mother (chloroplast matches very closely compared to decipiens)



# Still to do...

Here's the plan from here


0. Dot plots
L3 dot plots to check for duplicate regions in h1 and h2 - h1 vs h2 ,h1 against itself ,h2 against itself. remove potentially duplicated contigs.

2. Parental binning
1.1. Remove organelles
We're doing organelles separately, and we already have them perfect thanks to oatk. So the first job is to align the haplotypes to the two organelle genomes I now have, and remove all contigs that are ~perfectly covered by the organelle genomes. Then we're dealing with the nuclear genome only. (We might lose some small nuclear contigs that are copies of the cp and mt genomes, but that's OK). Remove any contig that is >=95% covered by organellar genomes is probably a good start, but the best thing is to look at the coverage and see.

this gives hap1_nuclear.fa and hap2_nuclear.fa

1.2. Alignment:
* align hap1 to hap2
* dotplot of all contigs from hap1 vs all from hap2
* get coverage table of all pairs, with H2 AND H1 coverage

This tells us which contigs are ~interchangeable and align to each other. Important for helping with parental assignment because we can expect homologues to be reciprocally assigned.

1.3. Kmer assignment
* Use the unique kmers we got for decipiens and virginea to get the dominant parent and dominant ratio for every contig in hap1 and hap2. 
* Decide on a cutoff ratio. With the top 11 contigs, the LOWEST dominant ratio was 1.85x, but the cutoff should be determined empirically by looking at the table.
* Contigs get one of three determinants: decipiens, virginea, unknown

* NB: Also use ROADIES with each contig as a 'species'. If we can see where they map vs. the ~35 existing euc genomes, we could assign haplotypes easily! Also each contig might just pop out cleanly on the species tree of other euc genomes. In an ideal world, homologous contigs would land on different parts of the tree, near their parent species.


Then combine kmer dominance with coverage table to see how far we can get in sorting the contigs into groups. For example, if two contigs are clearly homologues, and the hap1 contig is BELOW the ratio, but the hap2 homologue is ABOVE it (but the assignment is reciprocal) we can still sort it into parents. Each contig now ends up assigned to one of three groups: decipiens-parent, virginea-parent, unknown-parent

This can be done with a table like this:


contig1_ID
contig2_ID
contig1_length
contig2_length
length_ratio (1 / 2)
1v2_coverage
2v1_coverage
recip_best_hit (TRUE/FALSE)
kmer_ratio (highest / lowest coverage, noting which parent is highest, e.g. "2.5x decipiens")
parental_complementarity (TRUE if contig1 assigns to one parent, and contig2 assigns to the other)
contig1_final_assignment


1.4. Parental haplotypes
Create paternal (decipiens parent) and maternal (virginea parent; we know this from the organelles) haplotype assemblies.
Add organelle genome to maternal parent genome.

This gives us e_phylacis_paternal_contigs.fa and e_phylacis_maternal_contigs.fa

2. Scaffolding
Use ragtag to scaffold paternal against decipiens and maternal against virginea. This gives e_phylacis_paternal_scaffolds.fa and e_phylacis_maternal_scaffolds.fa

3. QC

For all four assemblies (paternal contigs, maternal contigs, paternal scaffolds, maternal scaffolds) QC including:
* compleasm
* Basic stats (N50 etc)
* map the reads back and check coverage
* ?

4. Annotation
The organelles are already assembled, but a basic annotation will be useful. 


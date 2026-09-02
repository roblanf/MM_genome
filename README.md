# Project: Assembly of the Meelup Mallee genome.

[give a brief intro to what this is, why, and the raw data]

# Environment

This repo uses Conda/Micromamba to help with reproducibility. Everything here is run on a Linux machine.

To set up and activate the environment:
```bash
micromamba create -f eucalypt_asm.yml
micromamba activate eucalypt_asm
```

# 01_QC: Quality control of raw data

First we establish the quality of the raw data. This can be done primarily with the GenomeQC tool, which generates graphics and descriptive summaries for genome assembly. The metrics computed in this section can be re-ran after all measures have been implemented, making for useful comparison and reflection on the quality control process.

## NanoPlot
First we use NanoPlot. This visualises and summarises long-read sequencing data. The main plot shows read-length and quality scores which you will see below:
Running NanoPlot:

```
bash /01_qc/scripts/run_nanoplot.sh
```
If this is your first time running NanoPlot, you will need to have the `pillow` installed in the environment. Add it here:
```
micromamba install pillow
```
Visualising the read length vs read quality results after running NanoPlot:

![Length vs Quality Scatter Plot](01_qc/results/nanoplot_test_output/LengthvsQualityScatterPlot_kde.png)

The displayed plot is for the tiny test dataset I've made. This is code to replicate my test data:

```
bash /01_qc/scripts/make_test_data.sh
```
*I will eventually replace this with the real outputs from the real raw_data once I'm done testing*

*I want to make a table of the important summary stats here*

NanoPlot will used again after the rest of the QC and filtering is complete to compare to these baseline results.


## K-mer counting and GenomeScope

GenomeQC allows for stats about the long reads' k-mers to be computed e.g. how many total k-mers and how many are unique. It is a good way to examine the size and complexity of the genome prior to assembly, and to measure the error rates in sequencing of the data we have.

The following code is used to compute basic KMC statistics and create a GenomeScope profile histogram useful for gleaning information about the Meelup Mallee genome:

```
bash /01_qc/scripts/run_kmc.sh
```
(output is currently for test data, placeholder to be replaced by full genome data soon)

KMC basics output:
Stats:
   No. of k-mers below min. threshold :            0
   No. of k-mers above max. threshold :            0
   No. of unique k-mers               :     49528710
   No. of unique counted k-mers       :     49528710
   Total no. of k-mers                :     61800000
   Total no. of reads                 :         5131
   Total no. of super-k-mers          :      8974769

GenomeScope profile histogram:

![GenomeScope K-mer Profile](01_qc/results/kmer_distribution/genomescope_results/linear_plot.png)



## Contamination check
The next step is to inspect the raw reads, determining the extent of contamination. One can check GC contamination easily, but a more thorough examination can be done with the Kraken2 tool.

** Roadblock ** Need to see if there is a rona server path for kraken2 already, otherwise need to download 8GB mini database.

Running a quicker contamination check on GC content for now:
```
bash /01_qc/scripts/run_gc_check.sh
```
(these are placeholder results from the test dataset)
9%      1
23%     2
24%     2
25%     3
26%     8
27%     4
28%     10
29%     14
30%     27
31%     24
32%     33
33%     77
34%     91
35%     178
36%     292
37%     429
38%     686
39%     790
40%     712
41%     536
42%     386
43%     215
44%     137
45%     101
46%     59
47%     25
48%     43
49%     42
50%     30
51%     27
52%     21
53%     20
54%     22
55%     28
56%     18
57%     5
58%     5
59%     11
60%     2
61%     2
62%     2
65%     1
66%     1
67%     1
68%     2
69%     1
71%     1
72%     2
73%     1
77%     1

There doesn't seem to be major evidence for contamination, as there is a unimodal distribution of GC content in the reads tested of around 38-40%. Spike at 59% could be an issue or just coincidence. The contigs can be better decontaminated once I implement the more complex Kraken2 tool.


# 02_filtering: Preparing reads for genome assembly by trimming and filtering

We have established the quality of the raw data. Next, the reads must be filtered and trimmed. This will involve the removal of adapters done by Porechop (trimming reads), and filtering reads based on thresholding quality score and other key metrics using Chopper (trims out those which don't pass the filters). Once this process has been conducted with parameters/thresholds determined suitable for the MM genome, the reads will be much better quality and well-prepared for genome assembly.

## Trimming adapters with Porechop
Porechop discovers unknown adapters in ONT reads such as those we have for the MM raw data and then trims them off. This is done using a k-mer based algorithm. It is a more specialised tool than Chopper, which will be used afterwards to more holistically filter out reads which are low in quality score or not suitable for genome assembly use based on a range of filters/metrics.

Running Porechop on the data:

If this is your first time running Porechop, you will need to install it to your environment like so:
```
micromamba install -c bioconda porechop -y 
```

```
bash /02_filtering/scripts/run_porechop.sh
```

Results:
Adapter trimming and chimeric read checking were run on `raw_data/tiny_test.fastq` (5,131 total reads).  (to be replaced by full dataset)
* **Start Adapters:** 3,347 / 5,131 reads trimmed (65.2%) — 130,752 bp removed
* **End Adapters:** 1,882 / 5,131 reads trimmed (36.7%) — 16,102 bp removed
* **Chimeric Reads Split:** 1 read (0.02%)
* **Total Adapters Removed:** 146,854 bp

There was a high adapter rate (~65% start, ~37% end in adapter sequences), so trimming out this noise in the reads is a highly effective and important measure prior to genome assembly.
Only one out of 5,131 reads had an adapter in the middle of the read, meaning this one had two distinct fragments joined together during library prep by mistake. The low rate suggests reads are in tact and lengths are genuine for almost all cases.

## Further trimming by filtering reads with Chopper



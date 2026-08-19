# Project: Assembly of the Meelup Mallee genome.

[give a brief intro to what this is, why, and the raw data]

# Environment

This repo uses conda to help with reproducibility. Everything here is run on a linux machine.

To set up a conda environment with everything you need, use the `environment.yml` file as follows:

```
conda env create -f environment.yml
```

# 01_QC: Quality control of raw data

First we establish the quality of the raw data. This... [short description]

First we use nanoplot to do x

```
bash /01_QC/scripts/run_nanoplot.sh
```

Then this

And that

## Results summary 




The first two sections, QC and filtering, are the pre-assembly modifications and checks that are run to improve the reads for assembly. For the genomeQC tools, I will use:

Porechop_ABI - discovers unknown adapters which must be trimmed, followed by...

Chopper - trims these adaptes and further trims reads via quality scores and other metrics.

Kraken2 - used for a formal contamination check.

NanoPlot - used before and after read cleaning/filtering to compare.


Then the next section is assembly, which will be done using hifiasm. To perform model selection, metrics such as BUSCO will be utilised.

Once a good and well thought out assembly has been chosen, the focus will be on parental assignment. The work done in this section will allow for a table to be created mapping all contigs to one of the Meelup Mallee parent species, *E. decipiens* or *E. virginea*.

Finally, the mitochondrial and chloroplast genomes can be considered and mapped to one of the parent species also.

## 01 QC

First, NanoPlot is run on the raw reads prior to cleaning and filtering. This is important for getting baseline information about read lengths, quality, etc. It will be useful to compare the metrics before and after read cleaning/filtering to assess the work done and improvement.
01

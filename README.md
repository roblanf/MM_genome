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

First we use NanoPlot. This visualises and summarises long-read sequencing data. The main plot shows read-length and quality scores which you will see below:
Running NanoPlot:

```
bash /01_QC/scripts/run_nanoplot.sh
```
Visualising the read length vs read quality results:

![Length vs Quality Scatter Plot](01_QC/01_NanoPlot_Raw/LengthvsQualityScatterPlot_kde.png)

Interpreting...

*other stats go here*

NanoPlot will used again after the rest of the QC and filtering is complete to compare to these baseline results.

Then more GenomeQC tools...

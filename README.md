# Project: Assembly of the Meelup Mallee genome.

The first two sections, QC and filtering, are the pre-assembly modifications and checks that are run to improve the reads for assembly. For the genomeQC tools, I will use:
Porechop_ABI - discovers unknown adapters which must be trimmed, followed by...
Chopper - trims these adaptes and further trims reads via quality scores and other metrics.
Kraken2 - used for a formal contamination check.
NanoPlot - used before and after read cleaning/filtering to compare.

Then the next sections are assembly, parental assignment and 


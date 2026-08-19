# Project: Assembly of the Meelup Mallee genome.

The first two sections, QC and filtering, are the pre-assembly modifications and checks that are run to improve the reads for assembly. For the genomeQC tools, I will use:

Porechop_ABI - discovers unknown adapters which must be trimmed, followed by...

Chopper - trims these adaptes and further trims reads via quality scores and other metrics.

Kraken2 - used for a formal contamination check.

NanoPlot - used before and after read cleaning/filtering to compare.


Then the next section is assembly, which will be done using hifiasm. To perform model selection, metrics such as BUSCO will be utilised.

Once a good and well thought out assembly has been chosen, the focus will be on parental assignment. The work done in this section will allow for a table to be created mapping all contigs to one of the Meelup Mallee parent species, *E. decipiens* or *E. virginea*.

Finally, the mitochondrial and chloroplast genomes can be considered and mapped to one of the parent species also.

## 01 QC

# Scaffolding

This README is purposed to store additional details required to conduct the `Scaffolding` step from the main <a href="../../README.md">`README.md`</a>.

## Table of Contents
- <a href="#06a">Scaffolding against eucalypt reference genomes</a>
- <a href="#06b">Checking telomere motif</a>

## <a id="06a">Scaffolding against eucalypt reference genomes</a>

```
# change working directory
cd 07_scaffold/

# virginea-like to virginea
ragtag.py scaffold -t 128 --mm2-params "-x asm20" -f 10000 -u -s 0.1 \
    -o virginea_like ../06_parental_assign/parents/roadies/Eucalyptus_virginea_GCA_014182375.1.fna \
    ../06_parental_assign/final/E_phylacis_virginea_like.fa

# decipiens-like to decipiens
ragtag.py scaffold -t 128 --mm2-params "-x asm20" -f 10000 -u -s 0.1 \
    -o decipiens_like ../06_parental_assign/parents/roadies/Eucalyptus_decipiens_GCA_014182575.1.fna \
    ../06_parental_assign/final/E_phylacis_decipiens_like.fa

# virginea-like to regnans
ragtag.py scaffold -t 128 --mm2-params "-x asm20" -f 10000 -u -s 0.1 \
    -o virginea_like_regnans ../06_parental_assign/parents/roadies/Eucalyptus_regnans_GCA_037349125.1.fna \
    ../06_parental_assign/final/E_phylacis_virginea_like.fa

# decipiens-like to regnans
ragtag.py scaffold -t 128 --mm2-params "-x asm20" -f 10000 -u -s 0.1 \
    -o decipiens_like_regnans ../06_parental_assign/parents/roadies/Eucalyptus_regnans_GCA_037349125.1.fna \
    ../06_parental_assign/final/E_phylacis_decipiens_like.fa

# change working directory
cd ..
```

| Metric               | virginea-like,<br>ref. *E. virginea* | virginea-like,<br>ref. *E. regnans* | decipiens-like,<br>ref. *E. decipiens* | decipiens-like,<br>ref. *E. regnans* |
| :------------------- | -----------------------------------: | ----------------------------------: | -------------------------------------: | -----------------------------------: |
| Input contigs        |                             277 |                            277 |                               245 |                             245 |
| Chromosome scaffolds |                              11 |                             11 |                                11 |                              11 |
| Contigs placed       |                             193 |                            159 |                               181 |                             153 |
| Placed (bp)          |                     534,221,985 |                    530,423,308 |                       565,214,048 |                     561,954,513 |
| Placed (%)           |                            98.2 |                           97.5 |                              98.6 |                            98.0 |
| Contigs unplaced     |                              84 |                            118 |                                64 |                              92 |
| `-s` < 0.1 (Mb)      |                       20 (5.19) |                      32 (5.96) |                         21 (2.43) |                       32 (5.99) |
| Unplaced (bp)        |                       9,686,931 |                     13,485,608 |                         8,269,866 |                      11,529,401 |
| Gaps (100 N)         |                             182 |                            147 |                               169 |                             142 |

## <a id="06b">Checking telomere motif</a>

```
# extract the 11 chromosomes
seqkit grep -r -p '^CM0245' 07_scaffold/virginea_like_a01/ragtag.scaffold.fasta  > 07_scaffold/virginea_like_a01/vir_chromosomes.fa
seqkit grep -r -p '^CM0246' 07_scaffold/decipiens_like_a01/ragtag.scaffold.fasta > 07_scaffold/decipiens_like_a01/dec_chromosomes.fa

# extract telomeres
tidk search -s AAACCCT -o vir_tidk -d 07_scaffold/ 07_scaffold/virginea_like_a01/vir_chromosomes.fa
tidk plot   -t 07_scaffold/vir_tidk_telomeric_repeat_windows.tsv -o vir_telomeres
tidk search -s AAACCCT -o dec_tidk -d 07_scaffold/ 07_scaffold/decipiens_like_a01/dec_chromosomes.fa
tidk plot   -t 07_scaffold/dec_tidk_telomeric_repeat_windows.tsv -o dec_telomeres
```

|                    decipiens-like                    |                   virginea-like                     |
| :--------------------------------------------------: | :-------------------------------------------------: |
| ![Decipiens-like telomeres](./figures/fig15_dec.svg) | ![Virginea-like telomeres](./figures/fig15_vir.svg) |

It seems like scaffolding does not really work that well, especially for the repetitive motifs on the virginea-like haplotype. Let's try another approach to scaffold the contigs...

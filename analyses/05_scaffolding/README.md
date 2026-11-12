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
tidk plot   -t 07_scaffold/vir_tidk_telomeric_repeat_windows.tsv -o 07_scaffold/vir_telomeres
tidk search -s AAACCCT -o dec_tidk -d 07_scaffold/ 07_scaffold/decipiens_like_a01/dec_chromosomes.fa
tidk plot   -t 07_scaffold/dec_tidk_telomeric_repeat_windows.tsv -o 07_scaffold/dec_telomeres
```

|                    decipiens-like                    |                   virginea-like                     |
| :--------------------------------------------------: | :-------------------------------------------------: |
| ![Decipiens-like telomeres](./figures/fig15_dec.svg) | ![Virginea-like telomeres](./figures/fig15_vir.svg) |

It seems like scaffolding does not really work that well, especially for the repetitive motifs on the virginea-like haplotype. I tried to do `ragtag.py correct` prior to scaffolding, but it made the scaffolds worse. Let's try another approach...

## <a id="06c">Binning contigs to each chromosome</a>
Let's try to bin the contigs prior to scaffolding, so every contig is placed on the right chromosome. Claude generated a <a href="./scripts/scaffold_per_chromosome.sh">script</a> to do this based on the `parent_of_contig.tsv` file from the `Assigning individual contigs to each putative parent` step (see the main <a href="../../README.md">`README.md`</a> file). We also need a file that stores the putative parent of each chromosome-scale contigs, which can be generated from `dv_scoring.tsv` file.

```
awk -F'\t' 'NR>1 {len=$3; gsub(",","",len)}
            NR>1 && (len+0 >= 10000000 || $1=="h2tg000036l") && ($15=="decipiens" || $15=="virginea") {print $1"\t"$15}' \
    06_parental_assign/alignment/dv_scoring.tsv > 06_parental_assign/alignment/parent_of_chromosome.tsv
```

### tidk

|                    decipiens-like                    |                   virginea-like                     |
| :--------------------------------------------------: | :-------------------------------------------------: |
| ![Decipiens-like telomeres](./figures/fig16_dec.svg) | ![Virginea-like telomeres](./figures/fig16_vir.svg) |

### decipiens-like
```
scaffold               anchor        left end                                          right end                                        internal-in-anchor
CM024616.1_RagTag      h1tg000003l   telomere 3416 copies; 16 comp. beyond (1.59 Mb)   telomere 2704 copies; 8 comp. beyond (1.29 Mb)   -
CM024613.1_RagTag      h1tg000295l   telomere 4092 copies; 10 comp. beyond (1.07 Mb)   no telomere                                      -
CM024611.1_RagTag      h1tg000007l   telomere 2989 copies                              telomere 3490 copies; 10 comp. beyond (1.21 Mb)  -
CM024610.1_RagTag      h1tg000014l   telomere 3261 copies; 1 comp. beyond (0.11 Mb)    telomere 3031 copies                             -
CM024614.1_RagTag      h1tg000004l   telomere 2241 copies; 2 comp. beyond (0.16 Mb)    telomere 5491 copies; 4 comp. beyond (0.25 Mb)   -
CM024615.1_RagTag      h1tg000001l   telomere 6223 copies; 17 comp. beyond (1.44 Mb)   telomere 5531 copies                             -
CM024609.1_RagTag      h1tg000009l   telomere 3625 copies; 11 comp. beyond (0.76 Mb)   telomere 9685 copies; 1 comp. beyond (0.06 Mb)   -
CM024619.1_RagTag      h2tg000006l   telomere 8149 copies; 1 comp. beyond (0.07 Mb)    no telomere                                      -
CM024617.1_RagTag      h2tg000011l   telomere 2900 copies; 4 comp. beyond (0.83 Mb)    telomere 4534 copies; 6 comp. beyond (0.48 Mb)   -
CM024618.1_RagTag      h2tg000003l   telomere 4288 copies; 9 comp. beyond (0.72 Mb)    telomere 5002 copies                             -
CM024612.1_RagTag      h2tg000001l   telomere 2522 copies; 6 comp. beyond (0.46 Mb)    telomere 3084 copies; 8 comp. beyond (0.62 Mb)   -

chromosome ends with a telomere: 20/22; components placed beyond a telomere: 114
```

### virginea-like
```
scaffold               anchor        left end                                          right end                                        internal-in-anchor
CM024521.1_RagTag      h2tg000008l   telomere 6775 copies; 1 comp. beyond (0.07 Mb)    telomere 3381 copies; 17 comp. beyond (1.52 Mb)  -
CM024518.1_RagTag      h2tg000004l   telomere 3019 copies                              telomere 1145 copies; 6 comp. beyond (0.69 Mb)   -
CM024516.1_RagTag      h2tg000010l   telomere 2636 copies; 7 comp. beyond (3.05 Mb)    telomere 4080 copies; 1 comp. beyond (0.09 Mb)   -
CM024515.1_RagTag      h2tg000005l   telomere 3133 copies; 2 comp. beyond (0.14 Mb)    telomere 4226 copies                             -
CM024519.1_RagTag      h2tg000012l   telomere 2537 copies                              telomere 4591 copies; 10 comp. beyond (1.01 Mb)  -
CM024520.1_RagTag      h2tg000002l   telomere 2831 copies; 11 comp. beyond (0.91 Mb)   telomere 4529 copies; 1 comp. beyond (0.09 Mb)   -
CM024514.1_RagTag      h2tg000007l   telomere 3814 copies                              telomere 8667 copies; 4 comp. beyond (0.30 Mb)   -
CM024524.1_RagTag      h1tg000006l   telomere 5317 copies                              telomere 2770 copies; 10 comp. beyond (1.24 Mb)  -
CM024522.1_RagTag      h1tg000005l   telomere 4125 copies; 7 comp. beyond (0.53 Mb)    telomere 4361 copies; 2 comp. beyond (0.30 Mb)   -
CM024517.1_RagTag      h1tg000012l   telomere 3057 copies; 14 comp. beyond (1.25 Mb)   telomere 4161 copies                             2701 copies at 1.82 Mb
CM024523.1_RagTag      h1tg000002l   telomere 5851 copies                              telomere 5680 copies; 10 comp. beyond (0.90 Mb)  -

chromosome ends with a telomere: 22/22; components placed beyond a telomere: 103
```

It seems like binning does help! Two chromosomes on the `decipiens-like` haplotype have missing telomeres, while all `virginea-like` chromosomes have both telomeres However, we still have two telomere peaks at the end of chromosomes 4 (`CM024517.1`) and 10 (`CM024523.1`) on the `virginea-like` haplotype. For chromosome 4, the two peaks are separated by a stretch of Ns, which was derived from Hifiasm. 
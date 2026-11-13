# Scaffolding

This README is purposed to store additional details required to conduct the `Scaffolding` step from the main <a href="../../README.md">`README.md`</a>.

## Table of Contents
- <a href="#06a">Scaffolding against *E. virginea* and *E. decipiens* reference genomes</a>
- <a href="#06b">Checking telomere motif</a>
- <a href="#06c">Scaffolding against *E. melliodora* reference genome</a>
- <a href="#06d">Trying different RagTag parameters</a>
    - <a href="#06d1">`-a`</a>
    - <a href="#06d2">`-f`</a>

## <a id="06a">Scaffolding against *E. virginea* and *E. decipiens* reference genomes</a>

```
# change working directory
cd 07_scaffold/

# virginea-like to virginea
ragtag.py scaffold -t 128 --mm2-params "-x asm20" -f 10000 -u -a 0.1 \
    -o virginea_like ../06_parental_assign/parents/roadies/Eucalyptus_virginea_GCA_014182375.1.fna \
    ../06_parental_assign/final/E_phylacis_virginea_like.dedup.fa

# decipiens-like to decipiens
ragtag.py scaffold -t 128 --mm2-params "-x asm20" -f 10000 -u -a 0.1 \
    -o decipiens_like ../06_parental_assign/parents/roadies/Eucalyptus_decipiens_GCA_014182575.1.fna \
    ../06_parental_assign/final/E_phylacis_decipiens_like.dedup.fa

# change working directory
cd ..
```

| Metric               | decipiens-like | virginea-like | 
| :------------------- | -------------: | ------------: |
| Chromosome scaffolds |             11 |            11 |
| Contigs placed       |            106 |           112 |
| Placed (bp)          |    557,947,248 |   523,001,838 |
| Placed (%)           |           98.9 |          98.7 |
| Contigs unplaced     |             42 |            45 |
| Unplaced (bp)        |      6,400,947 |     7,006,748 |
| Gaps (100 N)         |             94 |           101 |

## <a id="06b">Checking telomere motif</a>

```
# extract the 11 chromosomes
seqkit grep -r -p '^CM0245' 07_scaffold/virginea_like/ragtag.scaffold.fasta  > 07_scaffold/virginea_like/vir_chromosomes.fa
seqkit grep -r -p '^CM0246' 07_scaffold/decipiens_like/ragtag.scaffold.fasta > 07_scaffold/decipiens_like/dec_chromosomes.fa

# extract telomeres
tidk search -s AAACCCT -o vir_tidk -d 07_scaffold/ 07_scaffold/virginea_like/vir_chromosomes.fa
tidk plot   -t 07_scaffold/vir_tidk_telomeric_repeat_windows.tsv -o 07_scaffold/vir_telomeres
tidk search -s AAACCCT -o dec_tidk -d 07_scaffold/ 07_scaffold/decipiens_like/dec_chromosomes.fa
tidk plot   -t 07_scaffold/dec_tidk_telomeric_repeat_windows.tsv -o 07_scaffold/dec_telomeres
```

|                    decipiens-like                    |                   virginea-like                     |
| :--------------------------------------------------: | :-------------------------------------------------: |
| ![Decipiens-like telomeres](./figures/fig15_dec.png) | ![Virginea-like telomeres](./figures/fig15_vir.png) |

It seems like scaffolding works pretty well. Two chromosomes of the `decipiens-like` haplotype have missing telomeres. Moreover, chromosome 4 on `virginea-like` haplotype (`CM024517.1`) has two telomeres at the start. Let's extract the corresponding contig and map all contigs to it.

```
# extract the contig
grep "CM024517.1" 07_scaffold/virginea_like/ragtag.scaffold.agp | tail -1       # --> h1tg000012l

# extract the contig sequence
seqkit grep -p h1tg000012l 05_nuclear_qc/03_filtered/union.purged.fa > 05_nuclear_qc/03_filtered/h1tg000281l.fa

# map all contigs to the contig
minimap2 -cx asm20 05_nuclear_qc/03_filtered/union.purged.fa 05_nuclear_qc/03_filtered/h1tg000281l.fa \
    | awk '$11>=10000 {printf "%s\t%d-%d\tid=%.4f\tlen=%d\n", $6, $8, $9, $10/$11, $11}' | sort -t$'\t' -k4,4nr | head
```

This generated the following table, showing that there is a region on `h1tg000012l` that belongs to the `decipiens-like` haplotype, which is likely due to `hifiasm` wrongly stitched the two contigs together. Thus, we should remove the corresponding region (`35919968:35956317`), but retain the rest.
```
h1tg000012l	13475-35379892	    id=1.0000	len=35366417
h1tg000012l	35385041-35398699	id=1.0000	len=13658
h1tg000012l	35515957-35705868	id=1.0000	len=189911
h1tg000012l	35919968-35956317	id=1.0000	len=36349
h2tg000001l	37974011-38010248	id=0.9977	len=36305
```

Upon further inspection, we know that the `35398700:35515957` region consists of Ns. So, let's divide the contig into two to exclude the Ns.

```
# extract regions of interest
seqkit subseq --chr h1tg000012l -r 1:35398699        06_parental_assign/final/E_phylacis_virginea_like.dedup.fa | seqkit replace -p '.*' -r 'h1tg000012l'       > 06_parental_assign/final/h1tg000012l_main.fa
seqkit subseq --chr h1tg000012l -r 35515958:35705868 06_parental_assign/final/E_phylacis_virginea_like.dedup.fa | seqkit replace -p '.*' -r 'h1tg000012l_head'  > 06_parental_assign/final/h1tg000012l_head.fa

# update the haplotypes (remove `h1tg000012l`, add `h1tg000012l_main` and `h1tg000012l_head`)
seqkit grep -v -p h1tg000012l 06_parental_assign/final/E_phylacis_virginea_like.dedup.fa | \
    cat - 06_parental_assign/final/h1tg000012l_main.fa 06_parental_assign/final/h1tg000012l_head.fa \
    > 06_parental_assign/final/E_phylacis_virginea_like.dedup.split.fa
```

Now, let's redo the scaffolding for the `virginea-like` haplotype.

```
# change working directory
cd 07_scaffold/

# virginea-like to virginea
ragtag.py scaffold -t 128 --mm2-params "-x asm20" -f 10000 -u -a 0.1 \
    -o virginea_like_split ../06_parental_assign/parents/roadies/Eucalyptus_virginea_GCA_014182375.1.fna \
    ../06_parental_assign/final/E_phylacis_virginea_like.dedup.split.fa

# extract the 11 chromosomes
seqkit grep -r -p '^CM0245' virginea_like_split/ragtag.scaffold.fasta > virginea_like_split/vir_chromosomes.fa

# extract telomeres
tidk search -s AAACCCT -o vir_tidk_split -d ./ virginea_like_split/vir_chromosomes.fa
tidk plot   -t vir_tidk_split_telomeric_repeat_windows.tsv -o vir_telomeres_split

# change working directory
cd ..
```

|                virginea-like (before)                |                virginea-like (after)                |
| :--------------------------------------------------: | :-------------------------------------------------: |
| ![Virginea-like telomeres](./figures/fig15_vir.png)  | ![Virginea-like telomeres](./figures/fig16_vir.png) |

## <a id="06c">Scaffolding against *E. melliodora* reference genome</a>
I tried to scaffold the two haplotypes against a higher-quality *E. melliodora* reference genome. There might be a trade-off between the quality of the reference genome, as well as the genetic distance between the two species being compared.

```
# change working directory
cd 07_scaffold/

# virginea-like to melliodora
ragtag.py scaffold -t 64 --mm2-params "-x asm20" -f 10000 -u -a 0.1 \
    -o virginea_like_mel ../06_parental_assign/parents/roadies/Eucalyptus_melliodora_GCA_004368105.4.fna \
    ../06_parental_assign/final/E_phylacis_virginea_like.dedup.split.fa

# decipiens-like to melliodora
ragtag.py scaffold -t 64 --mm2-params "-x asm20" -f 10000 -u -a 0.1 \
    -o decipiens_like_mel ../06_parental_assign/parents/roadies/Eucalyptus_melliodora_GCA_004368105.4.fna \
    ../06_parental_assign/final/E_phylacis_decipiens_like.dedup.fa

# extract the 11 chromosomes
seqkit grep -r -p '^CM0234' virginea_like_mel/ragtag.scaffold.fasta  > virginea_like_mel/vir_chromosomes.fa
seqkit grep -r -p '^CM0234' decipiens_like_mel/ragtag.scaffold.fasta > decipiens_like_mel/dec_chromosomes.fa

# extract telomeres
tidk search -s AAACCCT -o vir_tidk_mel -d ./ virginea_like_mel/vir_chromosomes.fa
tidk plot   -t vir_tidk_mel_telomeric_repeat_windows.tsv -o vir_telomeres_mel
tidk search -s AAACCCT -o dec_tidk_mel -d ./ decipiens_like_mel/dec_chromosomes.fa
tidk plot   -t dec_tidk_mel_telomeric_repeat_windows.tsv -o dec_telomeres_mel

# change working directory
cd ..
```

| Metric               | decipiens-like | virginea-like | 
| :------------------- | -------------: | ------------: |
| Contigs placed       |             94 |           102 |
| Placed (bp)          |    553,425,857 |   522,408,916 |
| Contigs unplaced     |             54 |            56 |
| Unplaced (bp)        |     10,922,338 |     7,222,089 |
| Gaps (100 N)         |             82 |            91 |

|                    decipiens-like                    |                   virginea-like                     |
| :--------------------------------------------------: | :-------------------------------------------------: |
| ![Decipiens-like telomeres](./figures/fig17_dec.png) | ![Virginea-like telomeres](./figures/fig17_vir.png) |

It seems that the results are comparable with no clear winner. Let's try another approach by first mapping *all* contigs to *E. melliodora* reference genome, and then use the collapsed haplotype as the reference for each haplotype.

```
# change working directory
cd 07_scaffold/

# combine the two contigs
cat ../06_parental_assign/final/E_phylacis_virginea_like.dedup.split.fa \
    ../06_parental_assign/final/E_phylacis_decipiens_like.dedup.fa \
    > ../06_parental_assign/final/E_phylacis_union.dedup.split.fa

# union to melliodora
ragtag.py scaffold -t 64 --mm2-params "-x asm20" -f 10000 -u -a 0.1 \
    -o union_mel ../06_parental_assign/parents/roadies/Eucalyptus_melliodora_GCA_004368105.4.fna \
    ../06_parental_assign/final/E_phylacis_union.dedup.split.fa

# extract the 11 chromosomes
seqkit grep -r -p '^CM0234' union_mel/ragtag.scaffold.fasta  > union_mel/union_chromosomes.fa

# virginea-like to union
ragtag.py scaffold -t 64 --mm2-params "-x asm20" -f 10000 -u -a 0.1 \
    -o virginea_like_union union_mel/union_chromosomes.fa \
    ../06_parental_assign/final/E_phylacis_virginea_like.dedup.split.fa

# decipiens-like to union
ragtag.py scaffold -t 64 --mm2-params "-x asm20" -f 10000 -u -a 0.1 \
    -o decipiens_like_union union_mel/union_chromosomes.fa \
    ../06_parental_assign/final/E_phylacis_decipiens_like.dedup.fa

# extract the 11 chromosomes
seqkit grep -r -p '^CM0234' virginea_like_union/ragtag.scaffold.fasta  > virginea_like_union/vir_chromosomes.fa
seqkit grep -r -p '^CM0234' decipiens_like_union/ragtag.scaffold.fasta > decipiens_like_union/dec_chromosomes.fa

# extract telomeres
tidk search -s AAACCCT -o vir_tidk_union -d ./ virginea_like_union/vir_chromosomes.fa
tidk plot   -t vir_tidk_union_telomeric_repeat_windows.tsv -o vir_telomeres_union
tidk search -s AAACCCT -o dec_tidk_union -d ./ decipiens_like_union/dec_chromosomes.fa
tidk plot   -t dec_tidk_union_telomeric_repeat_windows.tsv -o dec_telomeres_union

# change working directory
cd ..
```

| Metric               | decipiens-like | virginea-like | 
| :------------------- | -------------: | ------------: |
| Contigs placed       |             71 |            84 |
| Placed (bp)          |    554,327,408 |   521,762,818 |
| Contigs unplaced     |             77 |            74 |
| Unplaced (bp)        |     10,020,787 |     7,868,187 |
| Gaps (100 N)         |             60 |            73 |

|                    decipiens-like                    |                   virginea-like                     |
| :--------------------------------------------------: | :-------------------------------------------------: |
| ![Decipiens-like telomeres](./figures/fig18_dec.png) | ![Virginea-like telomeres](./figures/fig18_vir.png) |

I think it looks much cleaner for most chromosomes, with telomere motifs found at both ends of the chromosome! Let's see if we can improve this further by changing the scaffolding parameters.

## <a id="06d">Trying different RagTag parameters</a>
As *E. melliodora* is not a direct close relative to sect. *Bisectae*, we need to be more careful on the initial scaffolding.

### <a id="06d1">`-a`</a>

```
# change working directory
cd 07_scaffold_a/

# iterate over different values of -a
for a in 0.0 0.1 0.2 0.3 0.4 0.5; do

    # union to melliodora
    ragtag.py scaffold -t 64 --mm2-params "-x asm20" -f 10000 -u -a ${a} \
        -o union_mel_a${a} ../06_parental_assign/parents/roadies/Eucalyptus_melliodora_GCA_004368105.4.fna \
        ../06_parental_assign/final/E_phylacis_union.dedup.split.fa

    # extract the 11 chromosomes
    seqkit grep -r -p '^CM0234' union_mel_a${a}/ragtag.scaffold.fasta \
        > union_mel_a${a}/union_chromosomes.fa

    # virginea-like to union
    ragtag.py scaffold -t 64 --mm2-params "-x asm20" -f 10000 -u -a ${a} \
        -o virginea_like_union_a${a} union_mel_a${a}/union_chromosomes.fa \
        ../06_parental_assign/final/E_phylacis_virginea_like.dedup.split.fa

    # decipiens-like to union
    ragtag.py scaffold -t 64 --mm2-params "-x asm20" -f 10000 -u -a ${a} \
        -o decipiens_like_union_a${a} union_mel_a${a}/union_chromosomes.fa \
        ../06_parental_assign/final/E_phylacis_decipiens_like.dedup.fa

    # extract the 11 chromosomes
    seqkit grep -r -p '^CM0234' virginea_like_union_a${a}/ragtag.scaffold.fasta \
        > virginea_like_union_a${a}/vir_chromosomes.fa
    seqkit grep -r -p '^CM0234' decipiens_like_union_a${a}/ragtag.scaffold.fasta \
        > decipiens_like_union_a${a}/dec_chromosomes.fa

    # extract telomeres
    tidk search -s AAACCCT -o vir_tidk_union_a${a} -d ./ virginea_like_union_a${a}/vir_chromosomes.fa
    tidk plot   -t vir_tidk_union_a${a}_telomeric_repeat_windows.tsv -o vir_telomeres_union_a${a}
    tidk search -s AAACCCT -o dec_tidk_union_a${a} -d ./ decipiens_like_union_a${a}/dec_chromosomes.fa
    tidk plot   -t dec_tidk_union_a${a}_telomeric_repeat_windows.tsv -o dec_telomeres_union_a${a}

done

# change working directory
cd ..
```

### Union

| Metric             | 0.0        | 0.1        | 0.2        | 0.3        | 0.4        | 0.5        |
| :----------------- | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| Contigs placed     | 201        | 174        | 165        | 161        | 159        | 157        |
| Placed (bp)        | 853369829  | 847845139  | 847077575  | 846003747  | 845696338  | 845490834  |
| Contigs unplaced   | 105        | 132        | 141        | 145        | 147        | 149        |
| Unplaced (bp)      | 240609371  | 246134061  | 246901625  | 247975453  | 248282862  | 248488366  |
| Gaps (100 N)       | 189        | 162        | 153        | 149        | 147        | 145        |

### Decipiens-like

| Metric             | 0.0        | 0.1        | 0.2        | 0.3        | 0.4        | 0.5        |
| :----------------- | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| Contigs placed     | 63         | 71         | 75         | 76         | 75         | 77         |
| Placed (bp)        | 555160364  | 554327408  | 554735695  | 554790709  | 554552711  | 554758215  |
| Contigs unplaced   | 85         | 77         | 73         | 72         | 73         | 71         |
| Unplaced (bp)      | 9187831    | 10020787   | 9612500    | 9557486    | 9795484    | 9589980    |
| Gaps (100 N)       | 52         | 60         | 64         | 65         | 64         | 66         |

| 0.0        | 0.1        | 0.2        | 0.3        | 0.4        | 0.5        |
| ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| ![](./figures/fig19_dec_a00.png) | ![](./figures/fig19_dec_a01.png) | ![](./figures/fig19_dec_a02.png) | ![](./figures/fig19_dec_a03.png) | ![](./figures/fig19_dec_a04.png) | ![](./figures/fig19_dec_a05.png) |

### Virginea-like

| Metric             | 0.0        | 0.1        | 0.2        | 0.3        | 0.4        | 0.5        |
| :----------------- | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| Contigs placed     | 86         | 84         | 85         | 86         | 86         | 86         |
| Placed (bp)        | 521223776  | 521762818  | 521823859  | 522714367  | 522714367  | 522714367  |
| Contigs unplaced   | 72         | 74         | 73         | 72         | 72         | 72         |
| Unplaced (bp)      | 8407229    | 7868187    | 7807146    | 6916638    | 6916638    | 6916638    |
| Gaps (100 N)       | 75         | 73         | 74         | 75         | 75         | 75         |

| 0.0        | 0.1        | 0.2        | 0.3        | 0.4        | 0.5        |
| ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| ![](./figures/fig19_vir_a00.png) | ![](./figures/fig19_vir_a01.png) | ![](./figures/fig19_vir_a02.png) | ![](./figures/fig19_vir_a03.png) | ![](./figures/fig19_vir_a04.png) | ![](./figures/fig19_vir_a05.png) |

It seems that the results are comparable when `-a` > 0.0. However, using `-a` ≥ 0.3 seems to move the telomere of `CM023406.3` closer to the middle. Let's use `-a 0.2` then!

### <a id="06d2">`-f`</a>

```
# change working directory
cd 07_scaffold_a02_f/

# iterate over -f values
for f_label in 5k 10k 25k 50k 75k 100k 250k 500k; do

    # convert f_label to bp
    f=$(( ${f_label%k} * 1000 ))

    # union to melliodora
    ragtag.py scaffold -t 64 --mm2-params "-x asm20" -f ${f} -u -a 0.2 \
        -o union_mel_f${f_label} ../06_parental_assign/parents/roadies/Eucalyptus_melliodora_GCA_004368105.4.fna \
        ../06_parental_assign/final/E_phylacis_union.dedup.split.fa

    # extract the 11 chromosomes
    seqkit grep -r -p '^CM0234' union_mel_f${f_label}/ragtag.scaffold.fasta \
        > union_mel_f${f_label}/union_chromosomes.fa

    # virginea-like to union
    ragtag.py scaffold -t 64 --mm2-params "-x asm20" -f ${f} -u -a 0.1 \
        -o virginea_like_union_f${f_label} union_mel_f${f_label}/union_chromosomes.fa \
        ../06_parental_assign/final/E_phylacis_virginea_like.dedup.split.fa

    # decipiens-like to union
    ragtag.py scaffold -t 64 --mm2-params "-x asm20" -f ${f} -u -a 0.1 \
        -o decipiens_like_union_f${f_label} union_mel_f${f_label}/union_chromosomes.fa \
        ../06_parental_assign/final/E_phylacis_decipiens_like.dedup.fa

    # extract the 11 chromosomes
    seqkit grep -r -p '^CM0234' virginea_like_union_f${f_label}/ragtag.scaffold.fasta \
        > virginea_like_union_f${f_label}/vir_chromosomes.fa
    seqkit grep -r -p '^CM0234' decipiens_like_union_f${f_label}/ragtag.scaffold.fasta \
        > decipiens_like_union_f${f_label}/dec_chromosomes.fa

    # extract telomeres
    tidk search -s AAACCCT -o vir_tidk_union_f${f_label} -d ./ virginea_like_union_f${f_label}/vir_chromosomes.fa
    tidk plot   -t vir_tidk_union_f${f_label}_telomeric_repeat_windows.tsv -o vir_telomeres_union_f${f_label}
    tidk search -s AAACCCT -o dec_tidk_union_f${f_label} -d ./ decipiens_like_union_f${f_label}/dec_chromosomes.fa
    tidk plot   -t dec_tidk_union_f${f_label}_telomeric_repeat_windows.tsv -o dec_telomeres_union_f${f_label}

done

# change working directory
cd ..
```

### Union

| Metric             | 5k         | 10k        | 25k        | 50k        | 75k        | 100k       |
| :----------------- | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| Contigs placed     | 170        | 165        | 157        | 147        | 143        | 142        |
| Placed (bp)        | 847748482  | 847077575  | 846561745  | 846726970  | 846184561  | 846092539  |
| Contigs unplaced   | 136        | 141        | 149        | 159        | 163        | 164        |
| Unplaced (bp)      | 246230718  | 246901625  | 247417455  | 247252230  | 247794639  | 247886661  |
| Gaps (100 N)       | 157        | 153        | 145        | 136        | 132        | 131        |

### Decipiens-like

| Metric             | 5k         | 10k        | 25k        | 50k        | 75k        | 100k       |
| :----------------- | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| Contigs placed     | 76         | 75         | 77         | 83         | 83         | 82         |
| Placed (bp)        | 554737707  | 554735695  | 554759526  | 557068269  | 557063298  | 556971276  |
| Contigs unplaced   | 72         | 73         | 71         | 65         | 65         | 66         |
| Unplaced (bp)      | 9610488    | 9612500    | 9588669    | 7279926    | 7284897    | 7376919    |
| Gaps (100 N)       | 65         | 64         | 66         | 72         | 72         | 71         |

| 5k         | 10k        | 25k        | 50k        | 75k        | 100k       |
| ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| ![](./figures/fig20_dec_f5k.png) | ![](./figures/fig20_dec_f10k.png) | ![](./figures/fig20_dec_f25k.png) | ![](./figures/fig20_dec_f50k.png) | ![](./figures/fig20_dec_f75k.png) | ![](./figures/fig20_dec_f100k.png) |

### Virginea-like

| Metric             | 5k         | 10k        | 25k        | 50k        | 75k        | 100k       |
| :----------------- | ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| Contigs placed     | 82         | 85         | 86         | 90         | 93         | 92         |
| Placed (bp)        | 521653739  | 521823859  | 521933233  | 522859729  | 523313177  | 523222714  |
| Contigs unplaced   | 76         | 73         | 72         | 68         | 65         | 66         |
| Unplaced (bp)      | 7977266    | 7807146    | 7697772    | 6771276    | 6317828    | 6408291    |
| Gaps (100 N)       | 71         | 74         | 75         | 79         | 82         | 81         |

| 5k         | 10k        | 25k        | 50k        | 75k        | 100k       |
| ---------: | ---------: | ---------: | ---------: | ---------: | ---------: |
| ![](./figures/fig20_vir_f5k.png) | ![](./figures/fig20_vir_f10k.png) | ![](./figures/fig20_vir_f25k.png) | ![](./figures/fig20_vir_f50k.png) | ![](./figures/fig20_vir_f75k.png) | ![](./figures/fig20_vir_f100k.png) |

Similarly, it seems that increasing the value of `-f` does not have a big impact, except on one chromosome on the `decipiens-like` haplotype. Let's stick with `-f 25000`.


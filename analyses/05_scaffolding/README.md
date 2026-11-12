# Scaffolding

This README is purposed to store additional details required to conduct the `Scaffolding` step from the main <a href="../../README.md">`README.md`</a>.

## Table of Contents
- <a href="#06a">Scaffolding against eucalypt reference genomes</a>
- <a href="#06b">Checking telomere motif</a>
- <a href="#06c">Excluding contigs beyond telomeres</a>

## <a id="06a">Scaffolding against eucalypt reference genomes</a>

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
| Input contigs        |            277 |           277 |
| Chromosome scaffolds |             11 |            11 |
| Contigs placed       |            193 |           159 |
| Placed (bp)          |    534,221,985 |   530,423,308 |
| Placed (%)           |           98.2 |          97.5 |
| Contigs unplaced     |             84 |           118 |
| Unplaced (bp)        |      9,686,931 |    13,485,608 |
| Gaps (100 N)         |            182 |           147 |

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

![Virginea-like telomeres](./figures/fig16_vir.png)

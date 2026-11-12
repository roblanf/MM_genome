#!/usr/bin/env bash
# Per-chromosome scaffolding of one parental genome of E. phylacis.
#
#   1. (no contig correction - see step 3)
#   2. chromosome bins: each contig assigned to a chromosome from the parental-assignment tables
#      (chromosome contigs by themselves; fragments via their homologue in the other haplotype)
#   3. one RagTag scaffold run per chromosome, with only that chromosome's contigs and reference sequence
#   4. final FASTA = 11 chromosome scaffolds + everything unplaced; combined AGP
#   5. telomere check (tidk + telomere_check.py) -> list of fragments placed beyond a telomere
#
# Usage:  bash scaffold_per_chromosome.sh vir     (or dec)
# Edit the INPUTS block once; everything else is derived.
set -euo pipefail

PARENT=${1:?usage: $0 vir|dec}

# ---------------------------------------------------------------- INPUTS ----------------------
BASE=/home/jeremias/meelupmallee/long_reads
PARENT_OF_CONTIG=$BASE/06_parental_assign/alignment/parent_of_contig.tsv
PARENT_OF_CHROM=$BASE/06_parental_assign/alignment/parent_of_chromosome.tsv
TELOMERE_CHECK=$BASE/scripts/telomere_check.py
THREADS=64
case $PARENT in
  vir) LABEL=virginea;  ASM=$BASE/06_parental_assign/final/E_phylacis_virginea_like.fa
       REF=$BASE/06_parental_assign/roadies/euc_refs/Eucalyptus_virginea_GCA_014182375.1.fa
       FIRST_AGP=$BASE/07_scaffold/virginea_like_a01/ragtag.scaffold.agp ;;
  dec) LABEL=decipiens; ASM=$BASE/06_parental_assign/final/E_phylacis_decipiens_like.fa
       REF=$BASE/06_parental_assign/roadies/euc_refs/Eucalyptus_decipiens_GCA_014182575.1.fa
       FIRST_AGP=$BASE/07_scaffold/decipiens_like_a01/ragtag.scaffold.agp ;;
  *)   echo "parent must be vir or dec"; exit 1 ;;
esac
WORK=$BASE/07c_scaffold_correct_bin/$PARENT
# ------------------------------------------------------------------------------------------------

for f in "$PARENT_OF_CONTIG" "$PARENT_OF_CHROM" "$TELOMERE_CHECK" "$ASM" "$REF" "$FIRST_AGP"; do
    [ -s "$f" ] || { echo "missing input: $f"; exit 1; }
done
mkdir -p "$WORK"/{bins,scaffolds,tidk} && cd "$WORK"
log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ---- 0. homologue pairs (hap1 contig <tab> hap2 contig), incl. the chromosome-2 pieces -----
cat > pairs.tsv << 'EOF'
h1tg000003l	h2tg000008l
h1tg000295l	h2tg000004l
h1tg000007l	h2tg000010l
h1tg000004l	h2tg000012l
h1tg000001l	h2tg000002l
h1tg000009l	h2tg000007l
h1tg000006l	h2tg000006l
h1tg000005l	h2tg000011l
h1tg000012l	h2tg000001l
h1tg000002l	h2tg000003l
h1tg000008l	h2tg000009l
h1tg000014l	h2tg000005l
h1tg000008l	h2tg000036l
EOF

# ---- 1. chromosome contig -> reference chromosome (from the first, whole-genome RagTag run) ---
log "anchor -> reference chromosome"
awk 'NR==FNR {a[$1]; next} $5=="W" && ($6 in a) {sub(/_RagTag$/, "", $1); print $6"\t"$1}' \
    "$PARENT_OF_CHROM" "$FIRST_AGP" > anchor_to_ref.tsv

# ---- 2. bins on original contig names: ref chromosome <tab> contig -------------------------
log "chromosome bins"
awk -F'\t' -v P="$LABEL" '
    NR==FNR {h1of[$2]=$1; h2of[$1]=$2; next}
    FNR>1 && $5==P {
        if ($4=="chromosome")   own=$1;
        else if ($4=="group")   own = ($1 ~ /^h1/) ? h1of[$6] : h2of[$6];
        else                    next;
        print own"\t"$1 }' pairs.tsv "$PARENT_OF_CONTIG" \
  | awk -F'\t' 'NR==FNR {ref[$1]=$2; next} ($1 in ref) {print ref[$1]"\t"$2}' anchor_to_ref.tsv - \
  | sort -u > chr_bins.tsv
cut -f1 chr_bins.tsv | sort | uniq -c | awk '{printf "    %s: %d contigs\n", $2, $1}'

# ---- 3. contig correction: OFF (tested and rejected: against a reference of another species,
#         ragtag correct breaks contigs at every real rearrangement; the chromosome contigs are
#         read-supported and intact). The 'correct/' outputs are written as identity so the rest
#         of the workflow is unchanged. Edit the parental FASTA by hand for any contig that has to
#         be split (e.g. h1tg000012l at its N run) before running this script.
mkdir -p correct
ln -sf "$(realpath "$ASM")" correct/ragtag.correct.fasta
[ -s "$ASM.fai" ] || samtools faidx "$ASM"
awk 'BEGIN{OFS="\t"} {print $1, 1, $2, 1, "W", $1, 1, $2, "+"}' "$ASM.fai" > correct/ragtag.correct.agp
awk '$5=="W" {print $6"\t"$1}' correct/ragtag.correct.agp > piece_to_contig.tsv       # identity: contig -> contig

# ---- 4. bins in corrected names; contigs outside every bin stay unplaced ---------------------
awk -F'\t' 'NR==FNR {bin[$2]=$1; next} ($2 in bin) {print bin[$2]"\t"$1}' chr_bins.tsv piece_to_contig.tsv > chr_bins_corrected.tsv
cut -f2 chr_bins_corrected.tsv | sort > binned.txt
seqkit seq -n -i correct/ragtag.correct.fasta | sort > all_corrected.txt
comm -23 all_corrected.txt binned.txt > unbinned.txt
log "$(wc -l < binned.txt) contig pieces binned, $(wc -l < unbinned.txt) unbinned (will be unplaced)"

# ---- 5. one RagTag run per chromosome -------------------------------------------------------
> chromosomes.fa; > unplaced.fa; > combined.agp
for chr in $(cut -f1 chr_bins_corrected.tsv | sort -u); do
    awk -v c=$chr '$1==c {print $2}' chr_bins_corrected.tsv > bins/$chr.txt
    seqkit grep -f bins/$chr.txt correct/ragtag.correct.fasta > bins/$chr.fa
    seqkit grep -p "$chr"          "$REF"                       > bins/$chr.ref.fa
    nq=$(grep -c '>' bins/$chr.fa); nr=$(grep -c '>' bins/$chr.ref.fa)
    log "$chr: $nq query sequences, $nr reference sequence"
    if [ "$nq" -eq 0 ] || [ "$nr" -eq 0 ]; then echo "    empty input for $chr, skipped"; continue; fi
    ragtag.py scaffold -t $THREADS --mm2-params "-x asm20" -f 10000 -u -a 0.1 \
        ${EXCLUDE:+-e $EXCLUDE} -o scaffolds/$chr bins/$chr.ref.fa bins/$chr.fa > scaffolds/$chr.log 2>&1
    # the chromosome scaffold is named after the reference sequence; everything else in the output is unplaced
    seqkit grep    -p "${chr}_RagTag" scaffolds/$chr/ragtag.scaffold.fasta >> chromosomes.fa
    seqkit grep -v -p "${chr}_RagTag" scaffolds/$chr/ragtag.scaffold.fasta >> unplaced.fa
    grep -v '^#' scaffolds/$chr/ragtag.scaffold.agp >> combined.agp
done
seqkit grep -f unbinned.txt correct/ragtag.correct.fasta >> unplaced.fa
cat chromosomes.fa unplaced.fa > E_phylacis_${LABEL}_like.scaffolds.fa
log "final assembly:"; seqkit stats chromosomes.fa unplaced.fa E_phylacis_${LABEL}_like.scaffolds.fa

# ---- 6. telomere check ------------------------------------------------------------------------
log "telomeres"
tidk search -s AAACCCT -o ${PARENT} -d tidk chromosomes.fa > /dev/null 2>&1
python3 "$TELOMERE_CHECK" --windows tidk/${PARENT}_telomeric_repeat_windows.tsv --agp combined.agp \
    --out telomeres.tsv --exclude-out beyond_telomere.txt
if [ -s beyond_telomere.txt ]; then
    log "$(wc -l < beyond_telomere.txt) fragments lie beyond a telomere (beyond_telomere.txt)."
    log "to leave them unplaced: EXCLUDE=$WORK/beyond_telomere.txt bash $0 $PARENT   (rerun; correct is reused)"
fi

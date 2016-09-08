#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
	echo -e ""
	echo -e "Usage: $(basename $0) <peak_summits.bed> <slop_radius> <fimo_threshold>"
	echo -e ""
	exit
fi

summits=$1
radius=$2
threshold=$3

# Annotations
LENGTHS="$HOME/data/bed/genome_3D7_v24.lengths"
REF="$HOME/data/refs/PlasmoDB-24_Pfalciparum3D7_Genome.fasta"
DREME="$HOME/projects/ap2I-chipseq/data/dreme/real_peaks/dreme.xml"

# output files
PREFIX="${summits%.*}"
SUMMITS_SLOP="${PREFIX}_${radius}bp_slop.bed"
SUMMITS_SLOP_FASTA="${SUMMITS_SLOP%.*}.fasta"
BACKGROUND="${SUMMITS_SLOP_FASTA%.*}.background"
FIMO_OUT="./motif1"
PEAK_OVERLAPS="${FIMO_OUT}/peak_overlaps.bed"
MAX_SCORES="${FIMO_OUT}/max_peak_scores.tsv"
DISTANCES="${FIMO_OUT}/distances.tsv"
TOP_HITS="${FIMO_OUT}/top_hits.bed"
TOP_HITS_SLOP="${TOP_HITS%.*}_30bps.bed"
TOP_HITS_SLOP_FASTA="${TOP_HITS%.*}_30bps.fasta"

# prepare peak summits by grabbing base pairs around them, retrieving the sequence, and calculating a background
bedtools slop -b "$radius" -i "$summits" -g "$LENGTHS" > "$SUMMITS_SLOP"
bedtools getfasta -fi "$REF" -bed "$SUMMITS_SLOP" -fo "$SUMMITS_SLOP_FASTA"
cat "$SUMMITS_SLOP_FASTA" | fasta-get-markov -m 2 > "$BACKGROUND"

# run fimo on the slopped peak summits file from above
fimo --bgfile "$BACKGROUND" --motif DTGCAC --oc "$FIMO_OUT" --thresh "$threshold" --parse-genomic-coord "$DREME" "$SUMMITS_SLOP_FASTA"

# annotate motif hits with peaks, keep only max scoring motifs for each peak, then keep only the closest motifs of the duplicates
bedtools intersect -a "$SUMMITS_SLOP" -b "$FIMO_OUT/fimo.gff" -wo | awk -v var="$radius" 'BEGIN {FS=OFS="\t"} {if($12 == "+") {print $6,$9,$10,$4,$11,$12,$9 - ($3-var)} else if($12 == "-") {print $6,$9,$10,$4,$11,$12,$10 - ($3-var)}}' > "$PEAK_OVERLAPS"
bedtools groupby -i "$PEAK_OVERLAPS" -g 4 -c 5 -o max > "$MAX_SCORES"
awk 'BEGIN {FS=OFS="\t"} FNR == NR {n[$1] = $2; next} n[$4] == $5 {print $0}' "$MAX_SCORES" "$PEAK_OVERLAPS" > /tmp/tmp.file
awk 'BEGIN {FS=OFS="\t"} function abs(v) {return v < 0 ? -v : v}; FNR == NR {if(!n[$4]) {n[$4] = $7} else {if(abs(n[$4]) > abs($7)) {n[$4] = $7}} next} {if(n[$4] == $7) {print $0}}' /tmp/tmp.file /tmp/tmp.file | uniq > "$TOP_HITS"
bedtools slop -l 12 -r 13 -g "$LENGTHS" -i "$TOP_HITS" | sort -k5 -gr > "$TOP_HITS_SLOP"
bedtools getfasta -fi "$REF" -bed "$TOP_HITS_SLOP" -fo "$TOP_HITS_SLOP_FASTA"

# run Rscript to generate heatmap
run_rscript plot_motif_heatmap.R -p "$TOP_HITS" -f motif1/top_hits_30bps.fasta -m motif1/motif_heatmap.eps -l motif1/logo.eps -d motif1/distances.eps

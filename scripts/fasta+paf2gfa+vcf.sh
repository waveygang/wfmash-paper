#!/bin/bash

# Inputs
FASTA=$1
PAF=$2
REFERENCE_PREFIX=$3
REFERENCE_SEPARATOR=$4
OUTPUT_PREFIX=$5
THREADS=$6

# Paths
SEQWISH=/home/guarracino/tools/seqwish/bin/seqwish
VG=/home/guarracino/tools/vg

echo "Graph induction"
GFA=${OUTPUT_PREFIX}.gfa
# -k 0` is for getting a lossless representation of the pairwise alignment
\time -v "$SEQWISH" -s "$FASTA" -p "$PAF" -g "$GFA" -k 0 -B 10M -t "$THREADS" -P

echo "Variant identification"
VCF=${OUTPUT_PREFIX}.vcf.gz
\time -v "$VG" deconstruct -P "$REFERENCE_PREFIX" -H "$REFERENCE_SEPARATOR" -e -a -t "$THREADS" "$GFA" | bgzip -c -@ "$THREADS" -l 9 > "$VCF" && tabix "$VCF"

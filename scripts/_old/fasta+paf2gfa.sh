#!/bin/bash

PATH_SEQWISH=/home/guarracino/tools/seqwish/bin/seqwish-ccfefb016fcfc9937817ce61dc06bbcf382be75e
PATH_VG=/home/guarracino/tools/vg

# Inputs
FASTA=$1
SEGMENT_LEN=$2
BLOCK_LEN=$3
IDENTITY=$4
NUM_HAPLOTYPES=$5
PATH_WFMASH=$6
PREFIX=$7
THREADS=$8

# Paths
PAF=${PREFIX}.s${SEGMENT_LEN}.l${BLOCK_LEN}.p$IDENTITY.n${NUM_HAPLOTYPES}.paf
GFA=${PREFIX}.s${SEGMENT_LEN}.l${BLOCK_LEN}.p$IDENTITY.n${NUM_HAPLOTYPES}.k0.gfa

echo "All-vs-all alignment"
\time -v $PATH_WFMASH "$FASTA" "$FASTA" -X -s "$SEGMENT_LEN" -l "$BLOCK_LEN" -p "$IDENTITY" -n "$NUM_HAPLOTYPES" -t "$THREADS" > "$PAF"

echo "Graph induction"
# -k 0` is for getting a lossless representation of the pairwise alignment
\time -v ${PATH_SEQWISH} -s "$FASTA" -p "$PAF" -g "$GFA" -k 0 -B 10M -t "$THREADS" -P

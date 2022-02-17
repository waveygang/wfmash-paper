#!/bin/bash

# Example:
# bash scripts/fasta+paf2vcf.sh data/LPA.subset.fa.gz 100k 300k 98 14 chm13 '_' ~/tools/wfmash/build/bin/wfmash-master master

# Inputs
FASTA=$1
SEGMENT_LEN=$2
BLOCK_LEN=$3
IDENTITY=$4
NUM_HAPLOTYPES=$5
REFERENCE_PREFIX=$6
SEPARATOR=$7
PATH_WFMASH=$8
OUTPUT_PREFIX=$9

FASTA_PREFIX=$(basename $FASTA .fa.gz)

PATH_SEQWISH=~/tools/seqwish/bin/seqwish-ccfefb016fcfc9937817ce61dc06bbcf382be75e
PATH_VG=~/tools/vg

# Paths
PAF=${OUTPUT_PREFIX}.$FASTA_PREFIX.s${SEGMENT_LEN}.l${BLOCK_LEN}.p$IDENTITY.n${NUM_HAPLOTYPES}.k16.paf
GFA=$PAF.k0.B10M.gfa
VCF=$GFA.vcf.gz

echo "All-vs-all alignment"
\time -v $PATH_WFMASH $FASTA $FASTA -X -s ${SEGMENT_LEN} -l ${BLOCK_LEN} -p $IDENTITY -n ${NUM_HAPLOTYPES} -k 16 -t 48 > $PAF

echo "Graph induction"
# -k 0` is for getting a lossless representation of the pairwise alignment
\time -v ${PATH_SEQWISH} -s $FASTA -p $PAF -g $GFA -k 0 -B 10M -t 48 -P

echo "Identify variants"
\time -v ${PATH_VG} deconstruct -e -a -P $REFERENCE_PREFIX -H $SEPARATOR $GFA -t 48 | bgzip -c > $VCF && tabix $VCF

#!/bin/bash

# Inputs
SAMPLE=$1
TRUTH_VCF=$2
QUERY_VCF=$3
EASY_REGIONS_BED=$4
HARD_REGIONS_BED=$5
OUTPUT_DIR=$6
PATH_VCF_PREPROCESS_SH=$7
REF_SDF=$8
REF="/lizardfs/erikg/HPRC/year1v2genbank/evaluation/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
echo "VCF processing"
bash "$PATH_VCF_PREPROCESS_SH" "$QUERY_VCF".vcf.gz "$SAMPLE" 50
rm "$QUERY_VCF".renamed.vcf.gz

FNAME=$(dirname "$QUERY_VCF")/$(basename "$QUERY_VCF".vcf.gz)
PREFIX="${FNAME%.vcf.gz}"
NORMALIZED_VCF=${PREFIX}.norm.max50.vcf.gz

echo "VCF evaluation"

rtg vcfeval \
    -t "$REF_SDF" \
    -b "$TRUTH_VCF" \
    -e "$EASY_REGIONS_BED" \
    -c "$NORMALIZED_VCF" \
    -T 12 \
    -o "$OUTPUT_DIR"/easy

rtg vcfeval \
    -t "$REF_SDF" \
    -b "$TRUTH_VCF" \
    -e "$HARD_REGIONS_BED" \
    -c "$NORMALIZED_VCF" \
    -T 12 \
    -o "$OUTPUT_DIR"/hard

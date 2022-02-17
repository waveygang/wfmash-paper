#!/bin/bash

# Variables
REF=/lizardfs/erikg/HPRC/year1v2genbank/evaluation/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
REF_SDF=/lizardfs/erikg/HPRC/year1v2genbank/evaluation/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.sdf

# Inputs
SAMPLE=$1
TRUTH_VCF=$2
QUERY_VCF=$3
EASY_REGIONS_BED=$4
HARD_REGIONS_BED=$5
OUTPUT_DIR=$6
PATH_VCF_PREPROCESS_SH=$7

echo "VCF renaming"
zcat $QUERY_VCF | sed 's/^grch38#//g' | bgzip -c -@ 48 > "$QUERY_VCF".renamed.vcf.gz && tabix "$QUERY_VCF".renamed.vcf.gz

echo "VCF processing"
bash ${PATH_VCF_PREPROCESS_SH} "$QUERY_VCF".renamed.vcf.gz "$SAMPLE" 50
rm "$QUERY_VCF".renamed.vcf.gz

FNAME=$(dirname "$QUERY_VCF")/$(basename "$QUERY_VCF".renamed.vcf.gz)
PREFIX="${FNAME%.vcf.gz}"
NORMALIZED_VCF=${PREFIX}.max50.chr1-22.vcf.gz

echo "VCF evaluation"

/gnu/store/3vmp4dw8y0r49h0hbjbgv3bckgvz4k0m-rtg-tools-3.11/rtg vcfeval \
    -t "$REF_SDF" \
    -b "$TRUTH_VCF" \
    -e "$EASY_REGIONS_BED" \
    -c "$NORMALIZED_VCF" \
    -T 48 \
    -o "$OUTPUT_DIR"/easy

/gnu/store/3vmp4dw8y0r49h0hbjbgv3bckgvz4k0m-rtg-tools-3.11/rtg vcfeval \
    -t "$REF_SDF" \
    -b "$TRUTH_VCF" \
    -e "$HARD_REGIONS_BED" \
    -c "$NORMALIZED_VCF" \
    -T 48 \
    -o "$OUTPUT_DIR"/hard

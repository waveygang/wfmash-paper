#!/bin/bash

# From https://github.com/wwliao/pangenome-utils/blob/main/preprocess_vcf.sh
# Usage: preprocess_vcf.sh <VCF file> <sample name> <max variant size>

VCF=$1
SAMPLE=$2
REF=$3
MAXSIZE=$4

FNAME=$(basename $VCF)
PREFIX=$(dirname $VCF)/"${FNAME%.vcf.gz}"
CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
MEM="10G"

bcftools view -a -s ${SAMPLE} -Ou ${VCF} \
    | bcftools norm -f ${REF} -c s -m - -Ou \
    | bcftools view -e 'GT="ref" | GT~"\."' -f 'PASS,.' -Ou \
    | bcftools sort -m ${MEM} -T bcftools-sort.XXXXXX -Ou \
    | bcftools norm -d exact -Oz -o ${PREFIX}.norm.vcf.gz \
    && bcftools index -t ${PREFIX}.norm.vcf.gz \
    && bcftools view -e "STRLEN(REF)>${MAXSIZE} | STRLEN(ALT)>${MAXSIZE}" \
                 -r ${CHROMS} -Oz -o ${PREFIX}.max${MAXSIZE}.chr1-22.vcf.gz \
                 ${PREFIX}.norm.vcf.gz \
    && bcftools index -t ${PREFIX}.max${MAXSIZE}.chr1-22.vcf.gz \
    && rm ${PREFIX}.norm.vcf.gz*

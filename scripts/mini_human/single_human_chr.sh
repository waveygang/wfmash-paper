#!/bin/bash

# Input
FASTA=$1
PAF=$2
CHR=$3 # e.g. chr1
PREFIX=$4
DIR_OUTPUT=$5

# Paths
PATH_REF_FASTA=/lizardfs/erikg/HPRC/year1v2genbank/evaluation/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
PATH_REF_SDF=/lizardfs/erikg/HPRC/year1v2genbank/evaluation/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.sdf
DIR_TRUTH_VCF=/lizardfs/guarracino/wfmash-paper/human/truth/
DIR_GENOMES=/lizardfs/guarracino/wfmash-paper/human/genomes

# Clean up
rm -rf /scratch/human-$PREFIX-$CHR
mkdir -p /scratch/human-$PREFIX-$CHR
cd /scratch/human-$PREFIX-$CHR || exit

bash /lizardfs/guarracino/wfmash-paper/scripts/fasta+paf2gfa+vcf.sh \
  $FASTA $PAF \
  grch38 '#' \
  $CHR 48

PATH_VCF_GZ=$CHR.vcf.gz

# Fix CHROM (to match the CHROM in the truth variant set)
mv $PATH_VCF_GZ xxx.vcf.gz
zcat xxx.vcf.gz | sed 's/grch38#1#//g' | bgzip -c -@ 48 > $PATH_VCF_GZ
rm xxx.vcf.gz

# Compare query/truth
for SAMPLE in HG00438 HG00621 HG00673 HG00733 HG00735 HG00741; do
  PATH_TRUTH_VCF_GZ=${DIR_TRUTH_VCF}/$SAMPLE.GRCh38_no_alt.deepvariant.vcf.gz

  bash /lizardfs/guarracino/wfmash-paper/scripts/vcf_preprocess.sh \
    $PATH_TRUTH_VCF_GZ \
    $SAMPLE \
    $PATH_REF_FASTA \
    $CHR \
    50
  TRUTH_VCF_GZ=$SAMPLE.GRCh38_no_alt.deepvariant.norm.max50.vcf.gz

  bash /lizardfs/guarracino/wfmash-paper/scripts/vcf_preprocess.sh \
    $PATH_VCF_GZ \
    $SAMPLE \
    $PATH_REF_FASTA \
    $CHR \
    50
  QUERY_VCF_GZ=$CHR.norm.max50.vcf.gz

  rtg vcfeval \
    -t $PATH_REF_SDF \
    -b $TRUTH_VCF_GZ \
    -c $QUERY_VCF_GZ \
    -T 48 \
    -o vcfeval/$SAMPLE/vg
done

DIR_OUTPUT_PREFIX=$DIR_OUTPUT/$PREFIX
mkdir -p $DIR_OUTPUT_PREFIX

cd /scratch/
mv /scratch/human-$PREFIX-$CHR $DIR_OUTPUT_PREFIX/$CHR

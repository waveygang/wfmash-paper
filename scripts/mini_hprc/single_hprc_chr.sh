#!/bin/bash

i=$1
FASTA=$2
PREFIX=$3
PATH_WFMASH=$4
OUTPUT_DIR=$5
PATH_EASY_REGIONS=$6
PATH_HARD_REGIONS=$7
PATH_FASTA_PAF_TO_VCF=$8
DIR_TRUTH_VCF=$9
PATH_VCF_PREPROCESS_SH=${10}
PATH_VCF_EVALUATION_SH=${11}

REF=/lizardfs/erikg/HPRC/year1v2genbank/evaluation/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
REF_SDF=/lizardfs/erikg/HPRC/year1v2genbank/evaluation/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.sdf

mkdir -p /scratch/$PREFIX-chr$i
cd /scratch/$PREFIX-chr$i

bash ${PATH_FASTA_PAF_TO_VCF} $FASTA 10k 50k 98 14 grch38 '#' ${PATH_WFMASH} ${PREFIX};

FASTA_PREFIX=$(basename $FASTA .fa.gz)
VCF=$PREFIX.${FASTA_PREFIX}.s10k.l50k.p98.n14.k16.paf.k0.B10M.gfa.vcf.gz

for SAMPLE in HG00438 HG00621 HG00673 HG00733 HG00735 HG00741; do
  TRUTH_VCF=${DIR_TRUTH_VCF}/$SAMPLE.GRCh38_no_alt.deepvariant.vcf.gz

  echo "VCF renaming"
  zcat $VCF | sed 's/^grch38#//g' | bgzip -c -@ 12 > "$VCF".renamed.vcf.gz && tabix "$VCF".renamed.vcf.gz

  bash ${PATH_VCF_EVALUATION_SH} $SAMPLE $TRUTH_VCF "$VCF".renamed.vcf.gz ${PATH_EASY_REGIONS} ${PATH_HARD_REGIONS} $SAMPLE ${PATH_VCF_PREPROCESS_SH} $REF_SDF
done

mv /scratch/$PREFIX-chr$i/* ${OUTPUT_DIR}
rm /scratch/$PREFIX-chr$i/ -rf

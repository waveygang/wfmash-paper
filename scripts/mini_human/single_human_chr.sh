#!/bin/bash

# Input
s=$1
l=$2
p=$3
n=$4
i=$5 # Num. chromosome
PATH_WFMASH=$6
PREFIX=$7
DIR_OUTPUT=$8

# Paths
vcfwave=/gnu/store/hkbkw85gjvsyqvx9vv9bw0ynmad989ag-vcflib-1.0.3+a36dbe9-11/bin/vcfwave # Jul 10, 2022
PATH_REF_FASTA=/lizardfs/erikg/HPRC/year1v2genbank/evaluation/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
PATH_REF_SDF=/lizardfs/erikg/HPRC/year1v2genbank/evaluation/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.sdf
DIR_TRUTH_VCF=/lizardfs/guarracino/wfmash-paper/human/truth/
DIR_GENOMES=/lizardfs/guarracino/wfmash-paper/human/genomes


CHR=chr$i

rm -rf /scratch/$CHR
mkdir -p /scratch/$CHR
cd /scratch/$CHR

PATH_FASTA_GZ=${DIR_GENOMES}/chr${i}hg00.masked.fa.gz

bash /lizardfs/guarracino/wfmash-paper/scripts/fasta+paf2gfa+vcf.sh \
  $PATH_FASTA_GZ \
  $s $l $p $n \
  grch38 '#' \
  $PATH_WFMASH \
  $CHR \
  48

PATH_VCF_GZ=$CHR.s$s.l$l.p$p.n$n.k0.vcf.gz

# Fix CHROM (to match the CHROM in the truth variant set)
mv $PATH_VCF_GZ xxx.vcf.gz
zcat xxx.vcf.gz | sed 's/grch38#1#//g' | bgzip -c -@ 48 > $PATH_VCF_GZ
rm xxx.vcf.gz

#echo "Realign REF/ALT alleles"
#PATH_VCF_WAVED_GZ=$CHR.s$s.l$l.p$p.n$n.k0.waved.vcf.gz
#vcfbub -l 0 -a 100000 --input $PATH_VCF_GZ | $vcfwave -I 1000 -t 48 | bgzip -c -@ 48 > $PATH_VCF_WAVED_GZ

# Compare query/truth
for SAMPLE in HG00438 HG00621 HG00673 HG00733 HG00735 HG00741; do
  PATH_TRUTH_VCF_GZ=${DIR_TRUTH_VCF}/$SAMPLE.GRCh38_no_alt.deepvariant.vcf.gz

  bash /lizardfs/guarracino/wfmash-paper/scripts/vcf_preprocess.sh \
    $PATH_TRUTH_VCF_GZ \
    $SAMPLE \
    $PATH_REF_FASTA \
    $i \
    50
  TRUTH_VCF_GZ=$SAMPLE.GRCh38_no_alt.deepvariant.norm.max50.vcf.gz

  bash /lizardfs/guarracino/wfmash-paper/scripts/vcf_preprocess.sh \
    $PATH_VCF_GZ \
    $SAMPLE \
    $PATH_REF_FASTA \
    $i \
    50
  QUERY_VCF_GZ=$CHR.s$s.l$l.p$p.n$n.k0.norm.max50.vcf.gz
  rtg vcfeval \
    -t $PATH_REF_SDF \
    -b $TRUTH_VCF_GZ \
    -c $QUERY_VCF_GZ \
    -T 48 \
    -o vcfeval/$SAMPLE/vg

#  bash /lizardfs/guarracino/wfmash-paper/scripts/vcf_preprocess.sh \
#    $PATH_VCF_WAVED_GZ \
#    $SAMPLE \
#    $PATH_REF_FASTA \
#    $i \
#    50
#  QUERY_WAVED_VCF_GZ=$CHR.s$s.l$l.p$p.n$n.k0.waved.norm.max50.vcf.gz
#  rtg vcfeval \
#    -t $PATH_REF_SDF \
#    -b $TRUTH_VCF_GZ \
#    -c $QUERY_WAVED_VCF_GZ \
#    -T 48 \
#    -o vcfeval/$SAMPLE/vcfwave
done

DIR_OUTPUT_PREFIX=$DIR_OUTPUT/$PREFIX
mkdir -p $DIR_OUTPUT_PREFIX

cd /scratch/
mv /scratch/$CHR $DIR_OUTPUT_PREFIX

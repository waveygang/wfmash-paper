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
DIR_GENOMES=/lizardfs/guarracino/wfmash-paper/tomato/genomes


CHR=chr$i

rm -rf /scratch/tomato-nucmer-$PREFIX-$CHR
mkdir -p /scratch/tomato-nucmer-$PREFIX-$CHR
cd /scratch/tomato-nucmer-$PREFIX-$CHR || exit

PATH_FASTA_GZ=${DIR_GENOMES}/chr$i.fa.gz

bash /lizardfs/guarracino/wfmash-paper/scripts/fasta+paf2gfa+vcf.sh \
  $PATH_FASTA_GZ \
  $s $l $p $n \
  SL5 '#' \
  "$PATH_WFMASH" \
  $CHR \
  48

PATH_GFA=$CHR.s$s.l$l.p$p.n$n.k0.gfa

# Compare query/truth
bash /lizardfs/guarracino/pggb-paper/scripts/gfa2evaluation.sh $PATH_GFA SL5 evaluation 48

DIR_OUTPUT_PREFIX=$DIR_OUTPUT/$PREFIX
mkdir -p $DIR_OUTPUT_PREFIX

cd /scratch/
mv /scratch/tomato-nucmer-$PREFIX-$CHR $DIR_OUTPUT_PREFIX/$CHR

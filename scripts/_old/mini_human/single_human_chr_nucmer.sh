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
DIR_GENOMES=/lizardfs/guarracino/wfmash-paper/human/genomes


CHR=chr$i

rm -rf /scratch/human-nucmer-$PREFIX-$CHR
mkdir -p /scratch/human-nucmer-$PREFIX-$CHR
cd /scratch/human-nucmer-$PREFIX-$CHR || exit

PATH_FASTA_GZ=${DIR_GENOMES}/chr${i}hg00.masked.fa.gz

bash /lizardfs/guarracino/wfmash-paper/scripts/fasta+paf2gfa.sh \
  $PATH_FASTA_GZ \
  $s $l $p $n \
  "$PATH_WFMASH" \
  $CHR \
  48

PATH_GFA=$CHR.s$s.l$l.p$p.n$n.k0.gfa

# Compare query/truth
bash /lizardfs/guarracino/pggb-paper/scripts/gfa2evaluation.sh $PATH_GFA chm13 evaluation 48

DIR_OUTPUT_PREFIX=$DIR_OUTPUT/$PREFIX
mkdir -p $DIR_OUTPUT_PREFIX

cd /scratch/
mv /scratch/human-nucmer-$PREFIX-$CHR $DIR_OUTPUT_PREFIX/$CHR

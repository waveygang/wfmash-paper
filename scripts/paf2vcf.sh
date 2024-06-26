#!/bin/bash

PAF=$1
NAME=$2
FASTA1=$3
FASTA2=$4
PANDEPTH=$5
DIR_OUTPUT=$6

mkdir /scratch/$NAME
cd /scratch/$NAME
wgatools filter -f paf -a 50000 $PAF | sort -k1,1V -k3,3n > $NAME.filter.50kb.paf

$PANDEPTH -i $NAME.filter.50kb.paf -a -o $NAME.filter.50kb.depth
zcat $NAME.filter.50kb.depth.SiteDepth.gz | awk '$3<=1' | awk '{print $1"\t"$2"\t"$2+1}'| bedtools merge -d 1 > $NAME.depth1.bed

# Capitalize to avoid wgatools to call variants where, for example, reference/alternate alleles are C/c
# Remove non-standard bases like R/Y that make wgatools call unhappy
REFNAME=$(basename $FASTA2 .fasta)
QRYNAME=$(basename $FASTA1 .fasta)
cat $FASTA2 | sed -E  '/^[^>]/ s/[^ATGCNatgcn]/N/g' | tr 'atgcn' 'ATGCN' > $REFNAME.fix.fa && samtools faidx $REFNAME.fix.fa
cat $FASTA1 | sed -E  '/^[^>]/ s/[^ATGCNatgcn]/N/g' | tr 'atgcn' 'ATGCN' > $QRYNAME.fix.fa && samtools faidx $QRYNAME.fix.fa

wgatools paf2maf -g $REFNAME.fix.fa -q $QRYNAME.fix.fa $NAME.filter.50kb.paf > $NAME.filter.50kb.maf

wgatools rename --prefixs "${REFNAME}#,${QRYNAME}#1#" $NAME.filter.50kb.maf > $NAME.filter.50kb.rename.maf
wgatools maf-index $NAME.filter.50kb.rename.maf
wgatools call --snp --svlen 0 $NAME.filter.50kb.rename.maf | sed "s/${REFNAME}#1#//g" | bgzip -@ 48 -c > $NAME.vcf.gz

bcftools view -T $NAME.depth1.bed $NAME.vcf.gz | bcftools view -e 'SVLEN>=50' -O z -o $NAME.depth1.short.vcf.gz
bcftools view -T $NAME.depth1.bed $NAME.vcf.gz | bcftools view -i 'SVLEN>=50' -O z -o $NAME.depth1.svs.vcf.gz

rm $NAME.filter.50kb*maf $NAME.filter.50kb.rename.maf.index $NAME.filter.50kb.depth.chr.stat.gz $NAME.filter.50kb.depth.SiteDepth.gz $NAME.filter.50kb.paf $REFNAME.fix.fa* $QRYNAME.fix.fa*

mkdir -p $DIR_OUTPUT
mv * $DIR_OUTPUT

cd /scratch
rm /scratch/$NAME -rf

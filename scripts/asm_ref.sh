#!/bin/bash
sample=$1
ref=$2
fmt="%C\\n%Us user %Ss system %P cpu %es total %MKb max memory"
time="/usr/bin/time"

sed "s/>/>${sample}#1#/g" ${sample}.fa |bgzip -@ 12 -c > ${sample}.panSN.fa.gz
samtools faidx ${sample}.panSN.fa.gz

sed "s/>/>${ref}#1#/g" ${ref}.fa |bgzip -@ 12 -c > ${ref}.panSN.fa.gz
samtools faidx ${ref}.panSN.fa.gz

## wfmash (have to follow panSN)
for p in 80 90
do
    #$time -f "$fmt" wfmash -s 10000 -p $p -n 1 -k 19 -H 0.001 -Y '#' -t 36 --hg-filter-ani-diff 30 ${ref}.fa ${sample}.fa > ${sample}.wfmash_p${p}s10kn1k19.paf
    $time -f "$fmt" wfmash -a -s 10000 -p $p -n 1 -k 19 -H 0.001 -Y '#' -t 36 --hg-filter-ani-diff 30 ${ref}.panSN.fa.gz ${sample}.panSN.fa.gz > ${sample}.wfmash_p${p}s10kn1k19.sam
done


## minimap2
for p in asm5 asm20
do
#    #$time -f "$fmt" minimap2 -x $p -c --eqx -t 36 ${ref}.fa ${sample}.fa > ${sample}.mm2_${p}.paf
    $time -f "$fmt" minimap2 -ax $p -c --eqx -t 36 ${ref}.fa ${sample}.fa > ${sample}.mm2_${p}.sam
done

## winnowmap2
#$time -f "$fmt" meryl count k=19 output merylDB asm1.fa
#$time -f "$fmt" meryl print greater-than distinct=0.9998 merylDB > repetitive_k19.txt
#$time -f "$fmt" winnowmap -W repetitive_k19.txt -ax asm20 --eqx -t 20 ./${ref}.fa ${sample}.fa > ${sample}.winowmap2_asm20.sam


## AnchorWave

$time -f "$fmt" anchorwave gff2seq -i ${ref}.gene.gff -r ${ref}.fa -o ${ref}.filter.cds.fa
$time -f "$fmt" minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 ${ref}.fa ${ref}.filter.cds.fa > ${ref}.cds.sam

$time -f "$fmt" minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 ${sample}.fa ${ref}.filter.cds.fa > ${sample}.cds.sam

## within-species with inverison

$time -f "$fmt" anchorwave genoAli -i ${ref}.gene.gff -as ${ref}.filter.cds.fa -r ${ref}.fa -a ${sample}.cds.sam -ar ${ref}.cds.sam -s ${sample}.fa -n ${sample}.aw.anchors -o ${sample}.aw.maf -f ${sample}.aw.m.maf -IV
rm ${sample}.aw.m.maf

## https://gitlab.com/mcfrith/last/-/blob/main/bin/maf-convert

$time -f "$fmt" python maf-convert sam ${sample}.aw.maf > ${sample}.aw.sam

# mummmer
# $time -f "$fmt" nucmer --sam-long ${sample}.nucmer.sam --maxmatch -c 500 -b 500 -l 100 -t 8 ${ref}.fa ${sample}.fa

## Variant calling
### SVs

grep "^@" ${sample}.mm2_asm20.sam > header

for i in `ls *.sam|grep -v "cds.sam"|sed -e "s/.sam//g" -e "s/${sample}.//g"`;
do
    echo "$i"
    if [[ $i == *aw* ]]; then
        grep -v "^@" ${sample}.${i}.sam >> combine.sam
        cat header <(grep -v "^@" ${sample}.${i}.sam)|samtools sort -@ 12 -O bam -o ${sample}.${i}.sorted.bam
    elif [[ $i == *wfmash* ]]; then
        sed -e "s/${sample}#1#/${sample}#${i}#1#/g" -e "s/${ref}#1#//g" ${sample}.${i}.sam |grep -v "^@" >> combine.sam
        sed "s/${ref}#1#//g" ${sample}.${i}.sam|samtools sort -@ 12 -O bam -o ${sample}.${i}.sorted.bam -
    elif [[ $i == *mm2* ]];then
        grep -v "^@" ${sample}.${i}.sam|sed -e "s/^Chr/${sample}#${i}#1#Chr/g" >> combine.sam
        samtools sort -@ 12 -O bam -o ${sample}.${i}.sorted.bam ${sample}.${i}.sam
    fi
done

cat header combine.sam|samtools sort -@ 12 -O bam -o ${sample}.aw_wfmash_mm2.sorted.bam
rm combine.sam
samtools index -@ 12 ${sample}.aw_wfmash_mm2.sorted.bam

# cutesv
cuteSV ${sample}.aw_wfmash_mm2.sorted.bam ${ref}.fa ${sample}.aw_wfmash_mm2.cutesv.vcf ./ \
    -s 1 --genotype --report_readid -p -1 -mi 500 -md 500 --min_mapq 0 \
    --max_cluster_bias_INS 1000 \
    --diff_ratio_merging_INS 0.9 \
    --max_cluster_bias_DEL 1000 \
    --diff_ratio_merging_DEL 0.5 \
    --min_size 50

bgzip -@ 12 ${sample}.aw_wfmash_mm2.cutesv.vcf

# filter wfmash-only
bcftools query -f "%CHROM\t%POS\t%INFO/RE\t%INFO/RNAMES\n" ${sample}.aw_wfmash_mm2.cutesv.vcf.gz|awk '$3 == 1'|grep "wfmash"|awk '{print $1"\t"$2"\t"$2+1"\t"$5}' > ${sample}.wfmash_only.bed

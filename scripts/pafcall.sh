#!/bin/bash
sample=$1
name=$2

~/software/wgatools/wgatools filter -f paf -a 50000 ${sample}_${name}.paf > ${sample}_${name}.filter.50kb.paf
~/software/PanDepth/bin/pandepth -i ./${sample}_${name}.filter.50kb.paf -a -o ${sample}_${name}
zcat ${sample}_${name}.SiteDepth.gz|awk '$3>1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > ${sample}_${name}_depth_lt1.bed
zcat ${sample}_${name}.SiteDepth.gz|awk '$3==1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > ${sample}_${name}_depth_eq1.bed
zcat ${sample}_${name}.SiteDepth.gz|awk '$3<1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > ${sample}_${name}_depth_ls1.bed
~/software/wgatools/wgatools paf2maf -g ref/Col-CC.fa -q fasta/${sample}.fa ${sample}_${name}.filter.50kb.paf > ${sample}_${name}.filter.50kb.maf
~/software/wgatools/wgatools rename --prefixs "Col-CC#1#,${sample}#1#" ${sample}_${name}.filter.50kb.maf > ${sample}_${name}.filter.50kb.rename.maf 
~/software/wgatools/wgatools mi ${sample}_${name}.filter.50kb.rename.maf
~/software/wgatools/wgatools call -s -l 1 -n ${sample}_${name} ${sample}_${name}.filter.50kb.rename.maf|sed "s/Col-CC#1#//g"|bgzip -@ 12 -c > ${sample}_${name}.vcf.gz

if [ -s ${sample}_${name}_depth_lt1.bed ]; then
    bcftools view -T ^${sample}_${name}_depth_lt1.bed -O z -o ${sample}_${name}.filter.vcf.gz ${sample}_${name}.vcf.gz
else
    mv ${sample}_${name}.vcf.gz ${sample}_${name}.filter.vcf.gz
fi

rm ${sample}_${name}.SiteDepth.gz
rm ${sample}_${name}.filter.50kb.rename.maf ${sample}_${name}.filter.50kb.maf
bgzip -@ 12 ${sample}_${name}.filter.50kb.paf

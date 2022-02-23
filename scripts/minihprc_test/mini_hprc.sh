#!/bin/bash

PREFIX=$1
PATH_WFMASH=$2

DIR_FASTA=/lizardfs/guarracino/HPRC/mini_dataset
DIR_REGIONS=/lizardfs/guarracino/HPRC/mini_dataset/union

# From https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/B581EBA7-8BDE-4C7C-9DEA-78B99A051155--Yale_HPP_Year1_Variant_Calls/
DIR_TRUTH_VCF=/lizardfs/guarracino/HPRC/mini_dataset/truth/

OUTPUT_DIR=/lizardfs/guarracino/HPRC/mini_dataset
PATH_FASTA_PAF_TO_VCF=/lizardfs/guarracino/HPRC/mini_dataset/fasta+paf2vcf.mod.sh
PATH_VCF_PREPROCESS_SH=/lizardfs/guarracino/HPRC/mini_dataset/vcf_preprocess.sh
PATH_VCF_EVALUATION_SH=/lizardfs/guarracino/HPRC/mini_dataset/vcf_evaluation.sh


( seq 22 -1 1; echo X; echo Y; echo M ) | while read i; do
        #echo $i

        if [[ ! -s ${DIR_REGIONS}/GRCh38_notinalldifficultregions.chr$i.bed.gz ]]; then
                zgrep -P "^chr$i\t" ${DIR_REGIONS}/GRCh38_notinalldifficultregions.bed.gz | bgzip > ${DIR_REGIONS}/GRCh38_notinalldifficultregions.chr$i.bed.gz;
        fi
        if [[ ! -s ${DIR_REGIONS}/GRCh38_alldifficultregions.chr$i.bed.gz ]]; then
                zgrep -P "^chr$i\t" ${DIR_REGIONS}/GRCh38_alldifficultregions.bed.gz | bgzip > ${DIR_REGIONS}/GRCh38_alldifficultregions.chr$i.bed.gz;
        fi

        OUTPUT_DIR_CHR=${OUTPUT_DIR}/$PREFIX/chr${i}

        if [ ! -d ${OUTPUT_DIR_CHR} ]; then
                #echo $PREFIX

                mkdir -p $OUTPUT_DIR_CHR

                FASTA=${DIR_FASTA}/chr${i}hg00.fa.gz
                PATH_EASY_REGIONS=${DIR_REGIONS}/GRCh38_notinalldifficultregions.chr$i.bed.gz
                PATH_HARD_REGIONS=${DIR_REGIONS}/GRCh38_alldifficultregions.chr$i.bed.gz

                sbatch -p workers -c 48 --job-name minihprc --wrap 'bash single_hprc_chr.sh '$FASTA' '$PREFIX' '${PATH_WFMASH}' '${OUTPUT_DIR_CHR}' '${PATH_EASY_REGIONS}' '${PATH_HARD_REGIONS}' '${PATH_FASTA_PAF_TO_VCF}' '${DIR_TRUTH_VCF}' '${PATH_VCF_PREPROCESS_SH}' '${PATH_VCF_EVALUATION_SH}';'
        fi
done


# Dataset generation
#( seq 22; echo X; echo Y; echo M ) | while read i; do
#        echo $i
#        samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/parts/chr${i}.pan.fa $(grep 'chr\|HG00' /lizardfs/erikg/HPRC/year1v2genbank/parts/chr${i}.pan.fa.fai | cut -f 1) | bgzip -@ 48 > chr${i}hg00.fa.gz
#        samtools faidx chr${i}hg00.fa.gz
#done

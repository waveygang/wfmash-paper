# usage
create a ref fasta and sdi

```shell
#-------------------------------------------------------------------
# Input
s=10k
l=50k
p=95
n=7

PATH_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-fixed-pen
PREFIX=fixed

PATH_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-poly-pen
PREFIX=poly
#-------------------------------------------------------------------

# Paths
vcfwave=/gnu/store/hkbkw85gjvsyqvx9vv9bw0ynmad989ag-vcflib-1.0.3+a36dbe9-11/bin/vcfwave # Jul 10, 2022

PATH_REF_FASTA=/lizardfs/guarracino/wfmash-paper/tomato/genomes/SL5.fa
PATH_REF_SDF=/lizardfs/guarracino/wfmash-paper/tomato/genomes/SL5.sdf
DIR_TRUTH_VCF=/lizardfs/guarracino/wfmash-paper/tomato/variants
DIR_GENOMES=/lizardfs/guarracino/wfmash-paper/tomato/genomes/

mkdir -p /scratch/$PREFIX
cd /scratch/$PREFIX

seq 11 11 | while read i; do
    CHR=chr$i
    
    mkdir -p $CHR
    cd $CHR
    
    PATH_FASTA_GZ=${DIR_GENOMES}/chr$i.fa.gz
    
#    bash /lizardfs/guarracino/wfmash-paper/scripts/fasta+paf2gfa+vcf.sh \
#      $PATH_FASTA_GZ \
#      $s $l $p $n \
#      SL5 '#' \
#      $PATH_WFMASH \
#      $CHR \
#      48
    
    PATH_VCF_GZ=$CHR.s$s.l$l.p$p.n$n.k0.vcf.gz
    
    # Fix CHROM (w.r.t. the truth variant set)
#    mv $PATH_VCF_GZ xxx.vcf.gz
#    zcat xxx.vcf.gz | sed 's/SL5#1#ch//g' | bgzip -c -@ 48 > $PATH_VCF_GZ
#    rm xxx.vcf.gz
  
    PATH_VCF_WAVED_GZ=$CHR.s$s.l$l.p$p.n$n.k0.waved.vcf.gz
    vcfbub -l 0 -a 100000 --input $PATH_VCF_GZ | vcfwave -I 1000 -t 48 | bgzip -c -@ 48 > $PATH_VCF_WAVED_GZ
    
    # Compare query/truth
    for SAMPLE in PP TS204 TS281 TS413 TS629 TS96; do
      PATH_TRUTH_VCF_GZ=${DIR_TRUTH_VCF}/$SAMPLE.hifi.vcf.gz
    
      # Preprocessing
      bash /lizardfs/guarracino/wfmash-paper/scripts/vcf_preprocess.sh \
        $PATH_VCF_GZ \
        $SAMPLE \
        $PATH_REF_FASTA \
        $i \
        50
      bash /lizardfs/guarracino/wfmash-paper/scripts/vcf_preprocess.sh \
        $PATH_TRUTH_VCF_GZ \
        $SAMPLE \
        $PATH_REF_FASTA \
        $i \
        50

      TRUTH_VCF_GZ=$SAMPLE.hifi.norm.max50.vcf.gz        
        
      QUERY_VCF_GZ=$CHR.s$s.l$l.p$p.n$n.k0.norm.max50.vcf.gz       
      rtg vcfeval \
        -t $PATH_REF_SDF \
        -b $TRUTH_VCF_GZ \
        -c $QUERY_VCF_GZ \
        -T 48 \
        -o vcfeval/$SAMPLE
        
      QUERY_WAVED_VCF_GZ=$CHR.s$s.l$l.p$p.n$n.k0.waved.norm.max50.vcf.gz  
      rtg vcfeval \
        -t $PATH_REF_SDF \
        -b $TRUTH_VCF_GZ \
        -c $QUERY_WAVED_VCF_GZ \
        -T 48 \
        -o vcfeval/$SAMPLE.waved
    done
    
    cd /scratch/$PREFIX
done

cd /scratch/
mv /scratch/$PREFIX /lizardfs/guarracino/wfmash-paper/tomato/tests
```

# Mini human

Run:

```shell
s=10k
l=50k
p=80
n=7
DIR_OUTPUT=/lizardfs/guarracino/wfmash-paper/human/tests


PATH_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-fixed-low-pen
PREFIX=fixed-0-4-6-1

PATH_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-fixed-med-pen
PREFIX=fixed-0-7-11-1

PATH_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-fixed-high-pen
PREFIX=fixed-0-11-17-1

PATH_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-poly-pen
PREFIX=poly

(seq 1 22; echo X; echo Y) | while read i; do
  sbatch -p workers -c 48 --job-name tomato --wrap "hostname; \time -v bash /lizardfs/guarracino/wfmash-paper/scripts/mini_human/single_human_chr.sh $s $l $p $n $i $PATH_WFMASH $PREFIX $DIR_OUTPUT"
done
```

Statistics:

```shell
cd /lizardfs/guarracino/wfmash-paper/human/tests/

(echo branch chromosome sample vcf tp.baseline tp.call fp fn precision recall f1.score | tr ' ' '\t';
  grep None */*/vcfeval/*/*/summary.txt | cut -f 1,2,4,5,6 -d '/' | sed 's/summary.txt://g' |tr -s ' ' | sed 's,/ *, ,g' | cut -f 5 --complement -d ' ' | tr ' ' '\t') \
  > mini_tomato.stats.tsv
```

Utils:

```shell
# Dataset generation
( echo 1 ) | while read i; do
  echo $i

  samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/parts/chr${i}.pan.fa $(grep 'chr\|HG00' /lizardfs/erikg/HPRC/year1v2genbank/parts/chr${i}.pan.fa.fai | cut -f 1) | \
    bgzip -@ 48 > chr${i}hg00.fa.gz
  samtools faidx chr${i}hg00.fa.gz
done

# Truth from https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/B581EBA7-8BDE-4C7C-9DEA-78B99A051155--Yale_HPP_Year1_Variant_Calls/
```

# Mini human

Run:

```shell
FASTA=/lizardfs/guarracino/wfmash-paper/human/genomes/chr20hg00.masked.fa.gz
PAF=/lizardfs/guarracino/wfmash-paper/human/alignments/chr20hg00.masked.paf
CHR=chr20
PREFIX=wfmash
DIR_OUTPUT=/lizardfs/guarracino/wfmash-paper/human/tests

sbatch -p workers -c 48 --job-name human --wrap "hostname; \time -v bash /lizardfs/guarracino/wfmash-paper/scripts/mini_human/single_human_chr.sh $FASTA $PAF $CHR $PREFIX $DIR_OUTPUT"
```

Statistics:

```shell
cd /lizardfs/guarracino/wfmash-paper/human/tests/

(echo branch chromosome sample vcf tp.baseline tp.call fp fn precision recall f1.score | tr ' ' '\t';
  grep None */*/vcfeval/*/*/summary.txt | cut -f 1,2,4,5,6 -d '/' | sed 's/summary.txt://g' |tr -s ' ' | sed 's,/ *, ,g' | cut -f 5 --complement -d ' ' | tr ' ' '\t') \
  > mini_human.stats.tsv
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

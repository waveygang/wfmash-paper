# Mini tomato

Run:

```shell
PATH_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-fixed-low-pen
PREFIX=fixed-0-4-6-1

PATH_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-fixed-med-pen
PREFIX=fixed-0-7-11-1

PATH_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-fixed-high-pen
PREFIX=fixed-0-11-17-1

PATH_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-poly-pen
PREFIX=poly

s=10k
l=50k
p=80
n=7
DIR_OUTPUT=/lizardfs/guarracino/wfmash-paper/tomato/tests/

seq 1 12 | while read i; do
  sbatch -p workers -c 48 --job-name tomato --wrap "hostname; \time -v bash /lizardfs/guarracino/wfmash-paper/scripts/mini_tomato/single_tomato_chr.sh 10k 50k 95 7 $i $PATH_WFMASH $PREFIX $DIR_OUTPUT"
done
```

Statistics:

```shell
cd /lizardfs/guarracino/wfmash-paper/tomato/tests/

ls | while read PREFIX; do
  echo $PREFIX
  
  seq 1 12 | while read i; do
    grep None vcfeval/*/summary
  done
done

(echo branch chromosome sample vcf tp.baseline tp.call fp fn precision recall f1.score | tr ' ' '\t';
  grep None */*/vcfeval/*/*/summary.txt | cut -f 1,2,4,5,6 -d '/' | sed 's/summary.txt://g' |tr -s ' ' | sed 's,/ *, ,g' | cut -f 5 --complement -d ' ' | tr ' ' '\t') \
  > mini_tomato.stats.tsv
```

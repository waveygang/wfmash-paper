# Mini tomato

Run:

```shell
s=20k
l=100k
p=80
n=6 # 7 - 1
DIR_OUTPUT=/lizardfs/guarracino/wfmash-paper/tomato/tests/


PATH_WFMASH="/home/guarracino/tools/wfmash/build/bin/wfmash-fixed-0-4-6-1 -w 512"
PREFIX=fixed-0-4-6-1

PATH_WFMASH="/home/guarracino/tools/wfmash/build/bin/wfmash-fixed-0-7-11-1 -w 512"
PREFIX=fixed-0-7-11-1

PATH_WFMASH="/home/guarracino/tools/wfmash/build/bin/wfmash-fixed-0-11-17-1 -w 512"

PATH_WFMASH="/home/guarracino/tools/wfmash/build/bin/wfmash-poly -w 512"
PREFIX=poly

(seq 1 22; echo X; echo Y) | while read i; do
  sbatch -p workers -c 48 --job-name tomato --wrap "hostname; \time -v bash /lizardfs/guarracino/wfmash-paper/scripts/mini_tomato/single_tomato_chr.sh $s $l $p $n $i '$PATH_WFMASH' $PREFIX $DIR_OUTPUT"
done
```

Statistics:

```shell
cd /lizardfs/guarracino/wfmash-paper/tomato/tests/

(echo branch chromosome sample vcf tp.baseline tp.call fp fn precision recall f1.score | tr ' ' '\t';
  grep None */*/vcfeval/*/*/summary.txt | cut -f 1,2,4,5,6 -d '/' | sed 's/summary.txt://g' |tr -s ' ' | sed 's,/ *, ,g' | cut -f 5 --complement -d ' ' | tr ' ' '\t') \
  > mini_tomato.stats.tsv
```

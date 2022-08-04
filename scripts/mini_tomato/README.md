# usage
create a ref fasta and sdi

```shell
PATH_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-fixed-pen
PREFIX=fixed
#PATH_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-poly-pen
#PREFIX=poly

s=10k
l=50k
p=95
n=7
DIR_OUTPUT=/lizardfs/guarracino/wfmash-paper/tomato/tests/

seq 12 12 | while read i; do
  sbatch -p workers -c 48 --job-name tomato --wrap "hostname; \time -v bash /lizardfs/guarracino/wfmash-paper/scripts/mini_tomato/single_tomato_chr.sh 10k 50k 95 7 $i $PATH_WFMASH $PREFIX $DIR_OUTPUT"
done
```

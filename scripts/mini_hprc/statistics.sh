prefix1=2f9e6a889
prefix2=wflign-div-bffbe87


# Runtimes

# Memory
grep 'Command\|Maximum' HG00438_*/sl* | grep $prefix1 | grep wfmash -A 1 | grep Maximum | cut -f 1,3 -d ':' > tmp1
grep 'Command\|Maximum' HG00438_*/sl* | grep $prefix2 | grep wfmash -A 1 | grep Maximum | cut -f 1,3 -d ':' > tmp2
#paste tmp1 tmp2 | awk '{print($1,($4-$2)/1024/1024)}' | column -t # GB
paste tmp1 tmp2 | awk '{print($1,($2/$4))}' | column -t
rm tmp1 tmp2

# Num. matching bases
ls HG00438_*/*/*paf | grep $prefix1 | while read f; do nmatches=$(< $f awk '{ sum += $10 } END { print sum }'); echo "$f: $nmatches" >> tmp1; done
ls HG00438_*/*/*paf | grep $prefix2 | while read f; do nmatches=$(< $f awk '{ sum += $10 } END { print sum }'); echo "$f: $nmatches" >> tmp2; done
paste tmp1 tmp2 | awk '{print($1,$2-$4)}' | cut -f 3 -d '/' | column -t
rm tmp1 tmp2

# F1-scores
paste \
	<(grep None */*/*/*txt | grep $prefix1 | tr -s ' ' | cut -f 1,9 -d ' ') \
	<(grep None */*/*/*txt | grep $prefix2 | tr -s ' ' | cut -f 1,9 -d ' ') |\
	awk '{print($1, $2-$4)}'




# Num. matching bases (wfa2-integration - master-2f9e6a8)
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr1hg00.s100k.l300k.p98.n14.k16.paf:   1681514
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr10hg00.s100k.l300k.p98.n14.k16.paf:  949830
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr11hg00.s100k.l300k.p98.n14.k16.paf:  287760
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr12hg00.s100k.l300k.p98.n14.k16.paf:  439593
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr13hg00.s100k.l300k.p98.n14.k16.paf:  913379
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr14hg00.s100k.l300k.p98.n14.k16.paf:  706493
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr15hg00.s100k.l300k.p98.n14.k16.paf:  1447868
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr16hg00.s100k.l300k.p98.n14.k16.paf:  2100110
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr17hg00.s100k.l300k.p98.n14.k16.paf:  1008072
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr18hg00.s100k.l300k.p98.n14.k16.paf:  634935
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr19hg00.s100k.l300k.p98.n14.k16.paf:  1071772
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr2hg00.s100k.l300k.p98.n14.k16.paf:   866856
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr20hg00.s100k.l300k.p98.n14.k16.paf:  716981
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr21hg00.s100k.l300k.p98.n14.k16.paf:  523111
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr22hg00.s100k.l300k.p98.n14.k16.paf:  1192435
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr3hg00.s100k.l300k.p98.n14.k16.paf:   1316787
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr4hg00.s100k.l300k.p98.n14.k16.paf:   1510068
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr5hg00.s100k.l300k.p98.n14.k16.paf:   729678
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr6hg00.s100k.l300k.p98.n14.k16.paf:   356216
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr7hg00.s100k.l300k.p98.n14.k16.paf:   1709086
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr8hg00.s100k.l300k.p98.n14.k16.paf:   1411058
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chr9hg00.s100k.l300k.p98.n14.k16.paf:   1486277
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chrMhg00.s100k.l300k.p98.n14.k16.paf:   0
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chrXhg00.s100k.l300k.p98.n14.k16.paf:   496055
2f9e6a889809b6cb90aae5c06e896db0c59e23c5.chrYhg00.s100k.l300k.p98.n14.k16.paf:   9384

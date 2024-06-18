NAME1=$1
NAME2=$2
THREADS=$3
FASTA1=$4
FASTA2=$5

nucmer --maxmatch -c 500 -b 500 -l 100 -p $NAME1-vs-$NAME2.nucmer.c500b500l100 -t $THREADS $FASTA2 $FASTA1
delta-filter -m -i 90 -l 1000 $NAME1-vs-$NAME2.nucmer.c500b500l100.delta > $NAME1-vs-$NAME2.nucmer.c500b500l100.mi90l1k.delta
paftools.js delta2paf $NAME1-vs-$NAME2.nucmer.c500b500l100.mi90l1k.delta > $NAME1-vs-$NAME2.nucmer.c500b500l100.mi90l1k.paf

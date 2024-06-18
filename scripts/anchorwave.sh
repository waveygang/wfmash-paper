NAME1=$1
NAME2=$2
THREADS=$3
FASTA1=$4
FASTA2=$5
ANNOTATION=$6

anchorwave gff2seq -i $ANNOTATION -r $FASTA2 -o $NAME2.filter.cds.fa
minimap2 -x splice -t $THREADS -k 12 -a -p 0.4 -N 20 $FASTA2 $NAME2.filter.cds.fa > $NAME2.cds.sam
minimap2 -x splice -t $THREADS -k 12 -a -p 0.4 -N 20 $FASTA1 $NAME2.filter.cds.fa > $NAME1.cds.sam
anchorwave genoAli -IV -t $THREADS -i $ANNOTATION -as $NAME2.filter.cds.fa -r $FASTA2 -a $NAME1.cds.sam -ar $NAME2.cds.sam -s $FASTA1 -n $NAME1.AnchorWave.default.anchors -o $NAME1.AnchorWave.default.maf -f $NAME1.AnchorWave.default.m.maf

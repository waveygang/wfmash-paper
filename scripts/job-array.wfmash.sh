#!/bin/bash

WFMASH=$1
PATH_COMBINATIONS=$2
DIR_ALIGNMENTS=$3
THREADS=$4

combinations=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PATH_COMBINATIONS)

IFS=$'\t' read -r FASTA1 FASTA2 p s l c k hg filter <<< "$combinations"

if [ "$l" == "5s" ]; then
    s_number=${s%k}  # Remove the trailing 'k'
    l_number=$((s_number * 5))
    l="${l_number}k"
else
    $l=0
fi

NAME1=$(basename $FASTA1 .fasta)
NAME2=$(basename $FASTA2 .fasta)
DIR_OUTPUT=$DIR_ALIGNMENTS/wfmash.p${p}.s${s}.l${l}.c${c}.k${k}.hg${hg}.$filter/$NAME1-vs-$NAME2

if [ "$filter" == "yes1to1" ]; then
    fil_ani="--one-to-one"
else
    fil_ani=""
fi


hostname
cd /scratch

\time -v $WFMASH -p $p -s $s -l $l -c $c -k $k --hg-filter-ani-diff $hg $fil_ani -t $THREADS $FASTA2 $FASTA1 > $NAME1-vs-$NAME2.wfmash.p${p}s${s}l${l}c${c}k${k}.hg${hg}${filter}.paf

mkdir -p $DIR_OUTPUT
mv $NAME1-vs-$NAME2.wfmash.p${p}s${s}l${l}c${c}k${k}.hg${hg}${filter}.paf $DIR_OUTPUT

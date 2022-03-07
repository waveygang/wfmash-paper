#!/bin/bash

# The resulting files will be:

# Input
assembly='CHM13_v1.1'
species='Homo sapiens'

#path_input_fasta=/lizardfs/guarracino/pggb_grant/genomes/chm13#chr8.fa
path_input_fasta=/lizardfs/guarracino/forerik/chrIV.fa

divergences=(
  0.001
)

presets=(
  'asm5'
)

lengths=(
  5000000
)

samples=({A..C})

threads=48

# Tools
run_wfmash=/home/guarracino/tools/wfmash/build/bin/wfmash-a7b3215b8d34f23f9865d2ea96dcdb36137cb246
run_fpa=~/tools/fpa/target/release/fpa-a273badb68c429683369051a1d389ce186f0dd6a
run_minimap2=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
run_meryl=~/tools/Winnowmap/bin-b373a7790aed2f5cb68ee63cb5f414ac9d63ec5a/meryl
run_winnomap2=~/tools/Winnowmap/bin-b373a7790aed2f5cb68ee63cb5f414ac9d63ec5a/winnowmap
run_seqwish=~/tools/seqwish/bin/seqwish-88cd0ea5f086cadfaf21c4c363d71536a1a7ea09
run_vg=/home/guarracino/tools/vg
run_rtg=/gnu/store/mriq5x6l7kzz51d1z64cvl5qx3d6ylc9-rtg-tools-3.11/rtg

# Preparation
name_input_fasta=$(basename $path_input_fasta .fa)
num_haplotypes=${#samples[@]}
variant_block=10 # Min. distance between simulated variants
path_input_sdf=${name_input_fasta}.sdf

## Mapping and Alignment

#minimap2
mkdir -p alignments/minimap2

for index in ${!divergences[*]}; do
  divergence=${divergences[$index]}
  preset=${presets[$index]}

  name_input_mt_fasta=${name_input_fasta}_$divergence.fa

  for idx_len in ${!lengths[*]}; do
    len=${lengths[$idx_len]}
    path_input_fa_gz=assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz
    dir_paf=alignments/minimap2
    path_name=${name_input_fasta}+samples_$divergence.l${len}.drop0.paf
    echo $path_input_fa_gz

    # Align and filter short alignments
    sbatch -p workers -c $threads --job-name minimap2 --wrap 'hostname; \time -v '$run_minimap2' '$path_input_fa_gz' '$path_input_fa_gz' -t '$threads' -c -x '$preset' -X 2> '$dir_paf'/'$name_input_fasta'.minimap2.'$divergence'.l'${len}'.log  > '$path_name'; for l in 50000 100000 200000 300000 500000; do cat '$path_name' | '$run_fpa' drop -l $l > '$dir_paf'/'${name_input_fasta}'+samples_'$divergence'.l'${len}'.drop$l.paf; done; mv '$path_name' '$dir_paf
  done
done

#winnomap2
mkdir -p alignments/winnomap2

$run_meryl count k=19 output ${name_input_fasta}.merylDB $path_input_fasta
$run_meryl print greater-than distinct=0.9998 ${name_input_fasta}.merylDB >${name_input_fasta}.repetitive_k19.txt

for index in ${!divergences[*]}; do
  divergence=${divergences[$index]}
  preset=${presets[$index]}

  name_input_mt_fasta=${name_input_fasta}_$divergence.fa

  for idx_len in ${!lengths[*]}; do
    len=${lengths[$idx_len]}
    path_input_fa_gz=assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz
    dir_paf=alignments/winnomap2
    path_name=${name_input_fasta}+samples_$divergence.l${len}.drop0.paf
    echo $path_input_fa_gz

    # Align and filter short alignments
    sbatch -p workers -c $threads --job-name winnomap2 --wrap 'hostname; \time -v '$run_winnomap2' -W '${name_input_fasta}'.repetitive_k19.txt '$path_input_fa_gz' '$path_input_fa_gz' -t '$threads' -c -x '$preset' -X 2> '$dir_paf'/'$name_input_fasta'.winnomap2.'$divergence'.l'${len}'.log  > '$path_name'; for l in 50000 100000 200000 300000 500000; do cat '$path_name' | '$run_fpa' drop -l $l > '$dir_paf'/'${name_input_fasta}'+samples_'$divergence'.l'${len}'.drop$l.paf; done; mv '$path_name' '$dir_paf
  done
done

#wfmash
mkdir -p alignments/wfmash

cwd=$(pwd)

for index in ${!divergences[*]}; do
  divergence=${divergences[$index]}
  #identity=$(echo "(1-$divergence) * 100 - 2" | bc -l)

  name_input_mt_fasta=${name_input_fasta}_$divergence.fa

  for idx_len in ${!lengths[*]}; do
    len=${lengths[$idx_len]}
    path_input_fa_gz=$cwd/assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz
    dir_paf=$cwd/alignments/wfmash

    for s in 5k 10k 20k 50k 100k; do
      for p in 95 98; do
        path_name=${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.paf
        echo $path_name

        sbatch -p workers -c $threads --job-name wfmash --wrap 'hostname; cd /scratch; \time -v '$run_wfmash' '$path_input_fa_gz' '$path_input_fa_gz' -t '$threads' -s '$s' -p '$p' -X -n '$num_haplotypes' 2> '$dir_paf'/'$name_input_fasta'.wfmash.'$divergence'.l'${len}'.log > '$path_name'; mv '$path_name' '$dir_paf
      done
    done
  done
done

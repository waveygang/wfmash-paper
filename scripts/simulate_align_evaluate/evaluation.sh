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
run_seqwish=~/tools/seqwish/bin/seqwish-88cd0ea5f086cadfaf21c4c363d71536a1a7ea09
run_vg=/home/guarracino/tools/vg
run_rtg=/gnu/store/mriq5x6l7kzz51d1z64cvl5qx3d6ylc9-rtg-tools-3.11/rtg

# Preparation
name_input_fasta=$(basename $path_input_fasta .fa)
num_haplotypes=${#samples[@]}
variant_block=10 # Min. distance between simulated variants
path_input_sdf=${name_input_fasta}.sdf

## Graph induction
mkdir -p graphs/wfmash/
mkdir -p graphs/minimap2/
mkdir -p graphs/winnomap2/

cwd=$(pwd)
cd /scratch

for index in ${!divergences[*]}; do
  divergence=${divergences[$index]}
  echo $divergence

  name_output=${name_input_fasta}_$divergence
  path_input_mt_fasta=samples/$name_output.samples.fa

  base_output=${name_input_fasta}+samples_$divergence

  for aligner in wfmash; do
    mkdir -p graphs/$aligner
    for idx_len in ${!lengths[*]}; do
      len=${lengths[$idx_len]}
      path_input_fa_gz="$cwd"/assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz

      for s in 5k 10k 20k 50k 100k; do
        for p in 95 98; do
          path_name=${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.paf

          path_input_paf="$cwd"/alignments/$aligner/${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.paf
          path_gfa="$cwd"/graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.gfa
          echo $path_gfa
          $run_seqwish -t $threads -s $path_input_fa_gz -p $path_input_paf -k 0 -g $path_gfa -B 50M -P
        done
      done
    done
  done
  mv graphs/wfmash/* "$cwd"/graphs/wfmash/

  for aligner in minimap2 winnomap2; do
    mkdir -p graphs/$aligner
    for idx_len in ${!lengths[*]}; do
      len=${lengths[$idx_len]}
      path_input_fa_gz="$cwd"/assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz
      for l in 0 50000 100000 200000 300000 500000; do
        path_input_paf="$cwd"/alignments/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.paf
        path_gfa="$cwd"/graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.gfa
        echo $path_gfa
        $run_seqwish -t $threads -s $path_input_fa_gz -p $path_input_paf -k 0 -g $path_gfa -B 50M -P
      done
    done
  done
  mv graphs/minimap2/* "$cwd"/graphs/minimap2/
  mv graphs/winnomap2/* "$cwd"/graphs/winnomap2/
done

cd $cwd

## Variant calling

for index in ${!divergences[*]}; do
  divergence=${divergences[$index]}
  echo $divergence

  name_output=${name_input_fasta}_$divergence
  path_input_mt_fasta=samples/$name_output.samples.fa

  base_output=${name_input_fasta}+samples_$divergence
  path_mutated_base_fa=base/$name_output.fa

  for aligner in wfmash; do
    mkdir -p variants/$aligner
    for idx_len in ${!lengths[*]}; do
      len=${lengths[$idx_len]}

      for s in 5k 10k 20k 50k 100k; do
        for p in 95 98; do
          path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.gfa
          path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.vcf.gz
          echo $path_vcf_gz
          $run_vg deconstruct -e -a -P 'chr' -H '#' $path_gfa -t $threads | bgzip -c >$path_vcf_gz && tabix $path_vcf_gz
        done
      done
    done
  done

  for aligner in minimap2 winnomap2; do
    mkdir -p variants/$aligner
    for idx_len in ${!lengths[*]}; do
      len=${lengths[$idx_len]}
      for l in 0 50000 100000 200000 300000 500000; do
        path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.gfa
        path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.vcf.gz
        echo $path_vcf_gz
        $run_vg deconstruct -e -a -P 'chr' -H '#' $path_gfa -t $threads | bgzip -c >$path_vcf_gz && tabix $path_vcf_gz
      done
    done
  done
done

#Evaluation:
$run_rtg format -o $path_input_sdf $path_input_fasta

for index in ${!divergences[*]}; do
  divergence=${divergences[$index]}
  echo $divergence

  name_output=${name_input_fasta}_$divergence

  path_mutated_base_fa=base/$name_output.fa

  for aligner in wfmash; do
    mkdir -p variants/$aligner
    for idx_len in ${!lengths[*]}; do
      len=${lengths[$idx_len]}

      for index_s in ${!samples[*]}; do
        sample=${samples[$index_s]}
        samplename=sample$sample
        path_truth_mut_vcf_gz=samples/$name_output.$samplename.vcf.gz

        for s in 5k 10k 20k 50k 100k; do
          for p in 95 98; do
            path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.vcf.gz

            echo $aligner $len $sample
            $run_rtg vcfeval -b $path_truth_mut_vcf_gz -c $path_vcf_gz -t $path_input_sdf --sample "$samplename,$samplename" -o evaluations/$aligner/${name_output}_$samplename.l${len}.s$s.p$p
          done
        done
      done
    done
  done

  for aligner in minimap2 winnomap2; do
    mkdir -p variants/$aligner
    for idx_len in ${!lengths[*]}; do
      for l in 0 50000 100000 200000 300000 500000; do
        len=${lengths[$idx_len]}
        path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.vcf.gz
        for index_s in ${!samples[*]}; do
          sample=${samples[$index_s]}
          samplename=sample$sample
          path_truth_mut_vcf_gz=samples/$name_output.$samplename.vcf.gz
          echo $aligner $len $l $sample
          $run_rtg vcfeval -b $path_truth_mut_vcf_gz -c $path_vcf_gz -t $path_input_sdf --sample "$samplename,$samplename" -o evaluations/$aligner/${name_output}_$samplename.l${len}.drop$l
        done
      done
    done
  done
done

grep None evaluations/*/*/summary.txt | grep sampleA | column -t

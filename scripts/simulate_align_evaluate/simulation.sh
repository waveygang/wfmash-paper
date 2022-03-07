#!/bin/bash

# The resulting files will be:
# - ./assemblies/{INPUT}+samples_{divergence}.l{contig_len}.fa.gz
# - ./samples/{INPUT}_{divergence}.sample{A}.vcf.gz

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
run_splitfa=~/tools/splitfa/target/release/splitfa-b6545722c429888dbe292fa2cdb265c6ea2c7004
run_mutation_simulator_py=/gnu/store/k2cs3b88a1ncvalavn85qqkil5h4d03v-mutation-simulator-2.0.3-1.9cb6bd2/bin/mutation-simulator.py
run_fastix=~/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834
run_bcftools=/gnu/store/ilksi2bdpsqkyf7r72lkcknndayrpjq3-bcftools-1.12/bin/bcftools
run_vcfrandomsample=/gnu/store/nczgf2j577bfw3g84f2b6r774bsmj1a6-vcflib-1.0.2/bin/vcfrandomsample
run_wfmash=/home/guarracino/tools/wfmash/build/bin/wfmash-a7b3215b8d34f23f9865d2ea96dcdb36137cb246
run_fpa=~/tools/fpa/target/release/fpa-a273badb68c429683369051a1d389ce186f0dd6a
run_minimap2=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
run_meryl=~/tools/Winnowmap/bin-b373a7790aed2f5cb68ee63cb5f414ac9d63ec5a/meryl
run_winnomap2=~/tools/Winnowmap/bin-b373a7790aed2f5cb68ee63cb5f414ac9d63ec5a/winnowmap
run_seqwish=~/tools/seqwish/bin/seqwish-88cd0ea5f086cadfaf21c4c363d71536a1a7ea09
run_rtg=/gnu/store/mriq5x6l7kzz51d1z64cvl5qx3d6ylc9-rtg-tools-3.11/rtg


# Preparation
name_input_fasta=$(basename $path_input_fasta .fa)
num_haplotypes=${#samples[@]}
variant_block=10 # Min. distance between simulated variants


# Mutated base generation (SNVs + small INDELs):
mkdir -p base
for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	# SNV/INDEL ratio ~ 9
	snv_rate=$(echo "$divergence * 10 * 0.9" | bc -l)
	ins_rate=$(echo "$divergence * 10 * 0.05" | bc -l)
	del_rate=$(echo "$divergence * 10 * 0.05" | bc -l)

	name_output_pre=${name_input_fasta}_${divergence}_pre
	prefix_mutations_vcf_pre=base/$name_output_pre
	$run_mutation_simulator_py $path_input_fasta --output $prefix_mutations_vcf_pre args --snp $snv_rate -titv 2 \
	  --insert $ins_rate --snpblock $variant_block --insertlength 3 --insertblock $variant_block \
	  --deletion $del_rate --deletionlength 3 --deletionblock $variant_block \
	  --assembly $assembly --species "$species" \
	  --sample ${name_input_fasta}_mt

    # Remove invalid variants (REF == ALT, in case on Ns in the reference)
	name_output=${name_input_fasta}_$divergence
	path_mutations_vcf=base/$name_output.vcf
    cat <(grep '^#' $prefix_mutations_vcf_pre.vcf) <(grep '^#' $prefix_mutations_vcf_pre.vcf -v | awk '{if ($4 != $5){print$0}}') > $path_mutations_vcf

    rm $prefix_mutations_vcf_pre.vcf
    mv base/$name_output_pre.fa base/$name_output.fa
done


# Sample generation:

mkdir -p samples

for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	name_output=${name_input_fasta}_$divergence
	path_mutations_vcf=base/$name_output.vcf

    path_mutated_base_fa=base/$name_output.fa
    samtools faidx $path_mutated_base_fa

    for index_s in ${!samples[*]}; do
      sample=${samples[$index_s]};
      samplename=sample$sample;

      # Clean temporary VCF
      rm samples/$name_output.$samplename.tmp.vcf -f

      cut -f 1 $path_mutated_base_fa.fai | while read seq_name; do
        len_input_fasta=$(grep -P "$seq_name\t" $path_input_fasta.fai | cut -f 2)
      	num_mutations_vcf=$(grep '^#' $path_mutations_vcf -v | grep -P "$seq_name\t" -c);
	    divergence_rate=$(echo "(($len_input_fasta / $variant_block) * $divergence * 10) / $num_mutations_vcf" | bc -l);

	    echo $seq_name $num_mutations_vcf $len_input_fasta $divergence_rate

        $run_vcfrandomsample <(cat <(grep '^#' $path_mutations_vcf) <(grep -P "^$seq_name\t" $path_mutations_vcf)) --rate $divergence_rate | sed "s/FORMAT\t${name_input_fasta}_mt/FORMAT\t${samplename}/g" | bgzip -@ $threads -c > samples/$name_output.$samplename.${seq_name}.vcf.gz
        tabix samples/$name_output.$samplename.${seq_name}.vcf.gz

        zgrep -v '^#' samples/$name_output.$samplename.${seq_name}.vcf >> samples/$name_output.$samplename.tmp.vcf
      done

      cat <(grep '^#' $path_mutations_vcf | sed "s/FORMAT\t${name_input_fasta}_mt/FORMAT\t${samplename}/g") <(grep '^#' -v samples/$name_output.$samplename.tmp.vcf) | bgzip -@ $threads -c > samples/$name_output.$samplename.vcf.gz;
      tabix samples/$name_output.$samplename.vcf.gz;

      rm samples/$name_output.$samplename.tmp.vcf
    done
done

for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}

	name_output=${name_input_fasta}_$divergence

	for index_s in ${!samples[*]}; do
	  sample=${samples[$index_s]};
	  samplename=sample$sample;
	  path_sample_mutations_vcf_gz=samples/$name_output.$samplename.vcf.gz;
	  $run_fastix -p "${samplename}#1#" <(cat $path_input_fasta | $run_bcftools consensus $path_sample_mutations_vcf_gz --sample $samplename) >> samples/$name_output.samples.fa
	done
done


# Assemblies generation:
mkdir -p assemblies

for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	name_output=${name_input_fasta}_$divergence
	path_input_mt_fasta=samples/$name_output.samples.fa

	base_output=${name_input_fasta}+samples_$divergence

	for idx_len in ${!lengths[*]}; do
	  len=${lengths[$idx_len]};
	  cat $path_input_fasta <($run_splitfa $path_input_mt_fasta -l $len-$len -s 1 2> assemblies/samples_$divergence.l${len}.paf) | bgzip -@ $threads -c > assemblies/$base_output.l${len}.fa.gz;
	  samtools faidx assemblies/$base_output.l${len}.fa.gz;
	done;
done

# Check
#zgrep '^>' assemblies/*.fa.gz -c
wc assemblies/*.fa.gz.fai -l

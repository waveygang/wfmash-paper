# Download the T2T assembly
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz


# Extract chromosomes
gunzip chm13.draft_v1.1.fasta.gz
samtools faidx chm13.draft_v1.1.fasta
samtools faidx chm13.draft_v1.1.fasta $(grep chr8 chm13.draft_v1.1.fasta.fai | cut -f 1) > chr8.fa


# Subregions
# Centromere (chr8:42,881,543-47,029,467)
samtools faidx chr8.fa chr8:42300000-47600000 | sed 's/[:-]/_/g' | sed 's/-/_/g' | sed 's/>chr8_42300000_47600000/>chr8_centromere/g' > chr8_centromere.fa

# 7 Mbp beta-defensin locus (chr8: 6,300,000-13,300,000)
samtools faidx chr8.fa chr8:5800000-13800000 | sed 's/[:-]/_/g' | sed 's/-/_/g' | sed 's/>chr8_5800000_13800000/>chr8_b_def/g' > chr8_b_def.fa


# Paths
# Tools
run_splitfa=~/git/fork/splitfa/target/release/splitfa # https://github.com/AndreaGuarracino/splitfa.git, random_segment_length branch
path_mutation_simulator_py=~/git/Mutation-Simulator/mutation-simulator.py
run_rtg=/home/tools/RealTimeGenomics/3.12/rtg
run_winnomap2=~/git/Winnowmap/bin/winnowmap
run_meryl=~/git/Winnowmap/bin/meryl

# Input
assembly='CHM13_v1.1'
species='Homo sapiens'
path_input_fasta=chm13.draft_v1.1.fa
path_input_sdf=chm13.draft_v1.1.sdf


# Variables
header_input_fasta=$(grep '>' $path_input_fasta | sed 's/>//g')
name_input_fasta=$(basename $path_input_fasta .fa)
len_input_fasta=$(cut $path_input_fasta.fai -f2)
threads=14

divergences=(
#  0.001
#  0.01
  0.05
#  0.15
#  0.20
)

presets=(
#  'asm5'
  'asm10'
#  'asm20'
)

lengths=(
	#200000
	#500000
	#1000000
	5000000
	#30000000
	#70000000
	#150000000
)

samples=({A..D})


# Mutated chromosome generation (SNVs + small INDELs)
variant_block=10

mkdir -p chromosomes
for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	# SNV/INDEL ratio = 9
	snv_rate=$(echo "$divergence * 2 * 0.9" | bc -l)
	ins_rate=$(echo "$divergence * 2 * 0.00005" | bc -l)
	del_rate=$(echo "$divergence * 2 * 0.00005" | bc -l)

	name_output_pre=${name_input_fasta}_${divergence}_pre
	prefix_mutations_vcf_pre=chromosomes/$name_output_pre
	python3 $path_mutation_simulator_py $path_input_fasta --output $prefix_mutations_vcf_pre args --snp $snv_rate -titv 2 --insert $ins_rate --snpblock $variant_block --insertlength 5 --insertblock $variant_block --deletion $del_rate --deletionlength 5 --deletionblock $variant_block --assembly $assembly --species "$species" --sample ${name_input_fasta}_mt

  # Remove invalid variants (REF == ALT, in case on Ns in the reference)
	name_output=${name_input_fasta}_$divergence
	path_mutations_vcf=chromosomes/$name_output.vcf

  cat <(grep '^#' $prefix_mutations_vcf_pre.vcf) <(grep '^#' $prefix_mutations_vcf_pre.vcf -v | awk '{if ($4 != $5){print$0}}') > $path_mutations_vcf

  rm $prefix_mutations_vcf_pre.vcf
  mv chromosomes/$name_output_pre.fa chromosomes/$name_output.fa
done


# Samples generation
mkdir -p samples

for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	name_output=${name_input_fasta}_$divergence
	path_mutations_vcf=chromosomes/$name_output.vcf

	num_mutations_vcf=$(grep '^#' $path_mutations_vcf -vc)
	divergence_rate=$(echo "(($len_input_fasta / $variant_block) * $divergence) / $num_mutations_vcf" | bc -l)

	for index_s in ${!samples[*]}; do sample=${samples[$index_s]}; samplename=sample$sample; vcflib vcfrandomsample $path_mutations_vcf --rate $divergence_rate | sed "s/FORMAT\t${name_input_fasta}_mt/FORMAT\t${samplename}/g" | bgzip -c > samples/$name_output.$samplename.vcf.gz; tabix samples/$name_output.$samplename.vcf.gz; done
done

for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}

	name_output=${name_input_fasta}_$divergence

	for index_s in ${!samples[*]}; do sample=${samples[$index_s]}; samplename=sample$sample; path_sample_mutations_vcf_gz=samples/$name_output.$samplename.vcf.gz; cat $path_input_fasta | bcftools consensus $path_sample_mutations_vcf_gz --sample $samplename | sed "s/>$header_input_fasta/>${samplename}#$header_input_fasta/g" >> samples/$name_output.samples.fa
	done
done

# Assemblies generation
mkdir -p seqs

for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	name_output=${name_input_fasta}_$divergence
	path_input_mt_fasta=samples/$name_output.samples.fa

	base_output=${name_input_fasta}+samples_$divergence

	for idx_len in ${!lengths[*]}; do len=${lengths[$idx_len]}; cat $path_input_fasta <($run_splitfa $path_input_mt_fasta -l $len-$len -s 1 2> seqs/samples_$divergence.l${len}.paf | perl /home/guarracino/Desktop/SharedFolder/wfmash-paper/scripts/split_contigs.perl -) | bgzip -@ $threads -c > seqs/$base_output.l${len}.fa.gz; samtools faidx seqs/$base_output.l${len}.fa.gz; sed 's/sample.#chr8_b_def\t/chr8_b_def\t/g' seqs/samples_$divergence.l${len}.paf -i; done;
done

# Check
zgrep '^>' seqs/*.fa.gz -c

# Alignments
num_haplotypes=${#samples[@]}


#minimap2
mkdir -p alignments/minimap2
for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	preset=${presets[$index]}

	name_input_mt_fasta=${name_input_fasta}_$divergence.fa

	for idx_len in ${!lengths[*]}; do len=${lengths[$idx_len]}; path_input_fa_gz=seqs/${name_input_fasta}+samples_$divergence.l${len}.fa.gz; path_paf=alignments/minimap2/${name_input_fasta}+samples_$divergence.l${len}.drop0.paf; echo $path_paf; \time -v minimap2 $path_input_fa_gz $path_input_fa_gz -t $threads -c -x $preset -X 2>> alignments/minimap2/minimap2.log > $path_paf; for l in 100000 400000 800000; do cat $path_paf | fpa drop -l $l > alignments/minimap2/${name_input_fasta}+samples_$divergence.l${len}.drop$l.paf; done; done;
done

#winnomap2
mkdir -p alignments/winnomap2

$run_meryl count k=19 output ${name_input_fasta}.merylDB $path_input_fasta
$run_meryl print greater-than distinct=0.9998 ${name_input_fasta}.merylDB > ${name_input_fasta}.repetitive_k19.txt

for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	preset=${presets[$index]}

	name_input_mt_fasta=${name_input_fasta}_$divergence.fa

	for idx_len in ${!lengths[*]}; do len=${lengths[$idx_len]}; path_input_fa_gz=seqs/${name_input_fasta}+samples_$divergence.l${len}.fa.gz; path_paf=alignments/winnomap2/${name_input_fasta}+samples_$divergence.l${len}.drop0.paf; echo $path_paf; \time -v $run_winnomap2 -W ${name_input_fasta}.repetitive_k19.txt $path_input_fa_gz $path_input_fa_gz -t $threads -c -x $preset -X 2>> alignments/winnomap2/winnomap2.log > $path_paf; for l in 100000 400000 800000; do cat $path_paf | fpa drop -l $l > alignments/winnomap2/${name_input_fasta}+samples_$divergence.l${len}.drop$l.paf; done; done;
done

#wfmash
mkdir -p alignments/wfmash
for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	identity=$(echo "(1-$divergence) * 100 - 1" | bc -l)

	name_input_mt_fasta=${name_input_fasta}_$divergence.fa

	for idx_len in ${!lengths[*]}; do len=${lengths[$idx_len]}; path_input_fa_gz=seqs/${name_input_fasta}+samples_$divergence.l${len}.fa.gz; path_paf=alignments/wfmash/${name_input_fasta}+samples_$divergence.l${len}.paf; echo $path_paf; \time -v wfmash $path_input_fa_gz $path_input_fa_gz -t $threads -s 100k -l 300k -p $identity -X -n $num_haplotypes 2>> alignments/wfmash/wfmash.log > $path_paf; done;
done


# Check
grep 'Elapsed (wall clock)' alignments/*/*log
grep 'Maximum resident set size (kbytes)' alignments/*/*log


# Pangenome graphs induction
for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	name_output=${name_input_fasta}_$divergence
	path_input_mt_fasta=samples/$name_output.samples.fa

	base_output=${name_input_fasta}+samples_$divergence

	for aligner in wfmash; do mkdir -p graphs/$aligner; for idx_len in ${!lengths[*]}; do len=${lengths[$idx_len]}; path_input_fa_gz=seqs/${name_input_fasta}+samples_$divergence.l${len}.fa.gz; path_input_paf=alignments/$aligner/${name_input_fasta}+samples_$divergence.l${len}.paf; path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.gfa; echo $path_gfa; seqwish -t $threads -s $path_input_fa_gz -p $path_input_paf -k 0 -g $path_gfa -B 9000000000; done; done;

	for aligner in minimap2 winnomap2; do mkdir -p graphs/$aligner; for idx_len in ${!lengths[*]}; do len=${lengths[$idx_len]}; path_input_fa_gz=seqs/${name_input_fasta}+samples_$divergence.l${len}.fa.gz; for l in 0 100000 400000 800000; do path_input_paf=alignments/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.paf; path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.gfa; echo $path_gfa; seqwish -t $threads -s $path_input_fa_gz -p $path_input_paf -k 0 -g $path_gfa -B 1000000000; done; done; done;
done


# Calling variants
for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	name_output=${name_input_fasta}_$divergence
	path_input_mt_fasta=samples/$name_output.samples.fa

	base_output=${name_input_fasta}+samples_$divergence

	for aligner in wfmash; do mkdir -p variants/$aligner; for idx_len in ${!lengths[*]}; do len=${lengths[$idx_len]}; path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.gfa; path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.vcf.gz; echo $path_vcf_gz; vg deconstruct -e -a -p $header_input_fasta -A sampleA#$header_input_fasta $path_gfa -t $threads | sed 's/#$header_input_fasta:/#${header_input_fasta}_/g' | bgzip -c > $path_vcf_gz && tabix $path_vcf_gz; done; done;

	for aligner in minimap2 winnomap2; do mkdir -p variants/$aligner; for idx_len in ${!lengths[*]}; do len=${lengths[$idx_len]}; for l in 0 100000 400000 800000; do path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.gfa; path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.vcf.gz; echo $path_vcf_gz; vg deconstruct -e -a -p $header_input_fasta -A sampleA#$header_input_fasta $path_gfa -t $threads | sed 's/#$header_input_fasta:/#${header_input_fasta}_/g' | bgzip -c > $path_vcf_gz && tabix $path_vcf_gz; done; done; done;
done


# Evaluation
$run_rtg format -o $path_input_sdf $path_input_fasta

for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	name_output=${name_input_fasta}_$divergence


	for aligner in wfmash; do mkdir -p variants/$aligner; for idx_len in ${!lengths[*]}; do len=${lengths[$idx_len]}; path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.vcf.gz; for index_s in ${!samples[*]}; do sample=${samples[$index_s]}; samplename=sample$sample; path_truth_mut_vcf_gz=samples/$name_output.$samplename.vcf.gz; $run_rtg vcfeval -b $path_truth_mut_vcf_gz -c $path_vcf_gz -t $path_input_sdf --sample "$samplename,$samplename#$header_input_fasta" -o evaluations/$aligner/${name_output}_$samplename.l${len}; done; done; done;

	for aligner in minimap2 winnomap2; do mkdir -p variants/$aligner; for idx_len in ${!lengths[*]}; do for l in 0 100000 400000 800000; do len=${lengths[$idx_len]}; path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.vcf.gz; for index_s in ${!samples[*]}; do sample=${samples[$index_s]}; samplename=sample$sample; path_truth_mut_vcf_gz=samples/$name_output.$samplename.vcf.gz; $run_rtg vcfeval -b $path_truth_mut_vcf_gz -c $path_vcf_gz -t $path_input_sdf --sample "$samplename,$samplename#$header_input_fasta" -o evaluations/$aligner/${name_output}_$samplename.l${len}.drop$l; done; done; done; done;
done

# Check sampleA
grep None evaluations/*/*/summary.txt | grep sampleA | column -t

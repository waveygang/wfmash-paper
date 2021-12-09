# Simulations

### Tools

Installation:

```shell
mkdir -p ~/tools
cd ~/tools

# splitfa
git clone --recursive https://github.com/AndreaGuarracino/splitfa.git
cd splitfa
git checkout b6545722c429888dbe292fa2cdb265c6ea2c7004
cargo build --release
cd ..

git clone --recursive https://github.com/natir/fpa.git
cd fpa
git checkout a273badb68c429683369051a1d389ce186f0dd6a
cargo build --release
cd ..

guix install bcftools

git clone --recursive https://gitlab.com/genenetwork/guix-bioinformatics.git
cd guix-bioinformatics
# Mutation-Simulator
GUIX_PACKAGE_PATH=. guix install mutation-simulator
# vcflib
GUIX_PACKAGE_PATH=. guix install vcflib
# minimap2
GUIX_PACKAGE_PATH=. guix install minimap2
# rtg vcfeval
GUIX_PACKAGE_PATH=. guix install rtg-tools
```

Versions:

```shell
(echo wfmash bcftools vcfrandomsample minimap2 mutation-simulator.py | tr ' ' '\n') | while read tool; do ls -l $(which $tool); done | cut -f 13 -d ' '

/gnu/store/nkfg1wg76zqaig43qgslkwcag9rb9fzz-wfmash-0.6.0+e9a5b02-17/bin/wfmash
/gnu/store/ilksi2bdpsqkyf7r72lkcknndayrpjq3-bcftools-1.12/bin/bcftools
/gnu/store/nczgf2j577bfw3g84f2b6r774bsmj1a6-vcflib-1.0.2/bin/vcfrandomsample
/gnu/store/fa1ng6dhv2lrqc8pfhg9iprs85ac055x-minimap2-2.23/bin/minimap2
/gnu/store/k2cs3b88a1ncvalavn85qqkil5h4d03v-mutation-simulator-2.0.3-1.9cb6bd2/bin/mutation-simulator.py
```

### Obtain the data

Create the main folder:

```shell
mkdir -p /lizardfs/guarracino/pggb_grant/
cd /lizardfs/guarracino/pggb_grant/
```

Download and prepare the references:

```shell
mkdir -p /lizardfs/guarracino/pggb_grant/genomes
cd /lizardfs/guarracino/pggb_grant/genomes

wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz
gunzip chm13.draft_v1.1.fasta.gz
samtools faidx chm13.draft_v1.1.fasta

cd ..
```

### Paths and variables

```shell
# Tools
run_splitfa=~/tools/splitfa/target/release/splitfa
run_fpa=~/tools/fpa/target/release/splitfa
run_mutation_simulator_py=/gnu/store/k2cs3b88a1ncvalavn85qqkil5h4d03v-mutation-simulator-2.0.3-1.9cb6bd2/bin/mutation-simulator.py
run_fastix=~/tools/fastix/target/release/fastix
run_bcftools=/gnu/store/ilksi2bdpsqkyf7r72lkcknndayrpjq3-bcftools-1.12/bin/bcftools
run_wfmash=/gnu/store/nkfg1wg76zqaig43qgslkwcag9rb9fzz-wfmash-0.6.0+e9a5b02-17/bin/wfmash
run_vcfrandomsample=/gnu/store/nczgf2j577bfw3g84f2b6r774bsmj1a6-vcflib-1.0.2/bin/vcfrandomsample
run_minimap2=/gnu/store/fa1ng6dhv2lrqc8pfhg9iprs85ac055x-minimap2-2.23/bin/minimap2
run_rtg=/gnu/store/mriq5x6l7kzz51d1z64cvl5qx3d6ylc9-rtg-tools-3.11/rtg

# Input
assembly='CHM13_v1.1'
species='Homo sapiens'
path_input_fasta=/lizardfs/guarracino/pggb_grant/genomes/chm13.draft_v1.1.fasta
path_input_sdf=/lizardfs/guarracino/pggb_grant/chm13.draft_v1.1.sdf

# Variables
threads=48
name_input_fasta=$(basename $path_input_fasta .fasta)

divergences=(
#  0.001
  0.01
#  0.05
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
	1000000
	5000000
	10000000
	#20000000
	#50000000
)
  
samples=({A..A})
num_haplotypes=${#samples[@]}

variant_block=10
```

Mutated base generation (SNVs + small INDELs):

```shell
cd /lizardfs/guarracino/pggb_grant/

mkdir -p base
for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	# SNV/INDEL ratio = 9
	snv_rate=$(echo "$divergence * 0.9" | bc -l)
	ins_rate=$(echo "$divergence * 0.05" | bc -l)
	del_rate=$(echo "$divergence * 0.05" | bc -l)

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
```

Sample generation:

```shell
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
      rm samples/$name_output.$samplename.tmp.vcf
      
      cut -f 1 $path_mutated_base_fa.fai | while read seq_name; do
        len_input_fasta=$(grep -P "$seq_name\t" $path_input_fasta.fai | cut -f 2)
      	num_mutations_vcf=$(grep '^#' $path_mutations_vcf -v | grep -P "$seq_name\t" -c);
	    divergence_rate=$(echo "(($len_input_fasta / $variant_block) * $divergence) / $num_mutations_vcf" | bc -l);
	    
	    #echo $seq_name $num_mutations_vcf $len_input_fasta $divergence_rate
	    
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
	  $run_fastix -p "${samplename}#" <(cat $path_input_fasta | $run_bcftools consensus $path_sample_mutations_vcf_gz --sample $samplename) >> samples/$name_output.samples.fa
	done
done
```

Assemblies generation:

```shell
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
	  cat $path_input_fasta <($run_splitfa $path_input_mt_fasta -l $len-$len -s 1 2> assemblies/samples_$divergence.l${len}.paf | perl split_contigs.perl -) | bgzip -@ $threads -c > assemblies/$base_output.l${len}.fa.gz;
	  samtools faidx assemblies/$base_output.l${len}.fa.gz;
	done;
done

# Check
#zgrep '^>' assemblies/*.fa.gz -c
wc assemblies/*.fa.gz.fai -l
```

## Mapping and Alignment

`minimap2`:

```shell
mkdir -p alignments/minimap2

for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	preset=${presets[$index]}

	name_input_mt_fasta=${name_input_fasta}_$divergence.fa

	for idx_len in ${!lengths[*]}; do
	  len=${lengths[$idx_len]};
	  path_input_fa_gz=assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz;
	  path_paf=alignments/minimap2/${name_input_fasta}+samples_$divergence.l${len}.drop0.paf;
	  echo $path_paf;
	  
	  # Align and filter short alignments
	  sbatch -p workers -c $threads --wrap '\time -v '$run_minimap2' '$path_input_fa_gz' '$path_input_fa_gz' -t '$threads' -c -x '$preset' -X 2> alignments/minimap2/minimap2.'$divergence'.l'${len}'.log > '$path_paf'; for l in 50000 100000 200000 300000 500000; do cat '$path_paf' | '$run_fpa' drop -l $l > alignments/minimap2/'${name_input_fasta}'+samples_'$divergence'.l'${len}'.drop'$l'.paf; done;'
	done;
done
```

`winnomap2`:

```shell
# ToDo to sbatcherize
mkdir -p alignments/winnomap2

$run_meryl count k=19 output ${name_input_fasta}.merylDB $path_input_fasta
$run_meryl print greater-than distinct=0.9998 ${name_input_fasta}.merylDB > ${name_input_fasta}.repetitive_k19.txt

for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	preset=${presets[$index]}

	name_input_mt_fasta=${name_input_fasta}_$divergence.fa

	for idx_len in ${!lengths[*]}; do
        len=${lengths[$idx_len]};
        path_input_fa_gz=assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz;
        path_paf=alignments/winnomap2/${name_input_fasta}+samples_$divergence.l${len}.drop0.paf;
        echo $path_paf;
        \time -v $run_winnomap2 -W ${name_input_fasta}.repetitive_k19.txt $path_input_fa_gz $path_input_fa_gz -t $threads -c -x $preset -X 2>> alignments/winnomap2/winnomap2.log > $path_paf;
        
        # Filter short alignments
        for l in 50000 100000 200000 300000 500000; do
          cat $path_paf | fpa drop -l $l > alignments/winnomap2/${name_input_fasta}+samples_$divergence.l${len}.drop$l.paf;
        done;
	done;
done
```

`wfmash`:

```shell
mkdir -p alignments/wfmash
for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	identity=$(echo "(1-$divergence) * 100 - 1" | bc -l)

	name_input_mt_fasta=${name_input_fasta}_$divergence.fa

	for idx_len in ${!lengths[*]}; do len=${lengths[$idx_len]};
        path_input_fa_gz=assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz;
        path_paf=alignments/wfmash/${name_input_fasta}+samples_$divergence.l${len}.paf;
        echo $path_paf;
        sbatch -p highmem -w octopus11 -c $threads --wrap '\time -v '$run_wfmash' '$path_input_fa_gz' '$path_input_fa_gz' -t '$threads' -s 100k -l 300k -p '$identity' -X -n '$num_haplotypes' 2>> alignments/wfmash/wfmash.'$divergence'.l'${len}'.log > '$path_paf;
	done;
done
```

```shell
# Check
grep 'Elapsed (wall clock)' alignments/*/*log
grep 'Maximum resident set size (kbytes)' alignments/*/*log
```

## Graph induction

```shell
for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	name_output=${name_input_fasta}_$divergence
	path_input_mt_fasta=samples/$name_output.samples.fa

	base_output=${name_input_fasta}+samples_$divergence

	for aligner in wfmash; do
	  mkdir -p graphs/$aligner;
	  for idx_len in ${!lengths[*]}; do
	    len=${lengths[$idx_len]};
	    path_input_fa_gz=assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz;
	    path_input_paf=alignments/$aligner/${name_input_fasta}+samples_$divergence.l${len}.paf;
	    path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.gfa;
	    echo $path_gfa;
	    seqwish -t $threads -s $path_input_fa_gz -p $path_input_paf -k 0 -g $path_gfa -B 50M;
	  done;
	done;

	for aligner in minimap2 winnomap2; do
	  mkdir -p graphs/$aligner;
	  for idx_len in ${!lengths[*]}; do
	    len=${lengths[$idx_len]};
	    path_input_fa_gz=assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz;
	    for l in 0 50000 100000 200000 300000 500000; do
	      path_input_paf=alignments/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.paf;
	      path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.gfa;
	      echo $path_gfa;
	      seqwish -t $threads -s $path_input_fa_gz -p $path_input_paf -k 0 -g $path_gfa -B 50M;
	    done;
	  done;
    done;
done
```

## Variant calling

```shell
for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	name_output=${name_input_fasta}_$divergence
	path_input_mt_fasta=samples/$name_output.samples.fa

	base_output=${name_input_fasta}+samples_$divergence
    path_mutated_base_fa=base/$name_output.fa
    cut -f 1 $path_mutated_base_fa.fai | while read seq_name; do
        for aligner in wfmash; do
          mkdir -p variants/$aligner;
          for idx_len in ${!lengths[*]}; do
            len=${lengths[$idx_len]};
            path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.gfa;
            path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.${seq_name}.vcf.gz;
            echo $path_vcf_gz;
            vg deconstruct -e -a -p $seq_name -A sampleA $path_gfa -t $threads | bgzip -c > $path_vcf_gz && tabix $path_vcf_gz;
            done
        done;
        
        for aligner in minimap2 winnomap2; do
          mkdir -p variants/$aligner;
          for idx_len in ${!lengths[*]}; do
            len=${lengths[$idx_len]};
            for l in 0 50000 100000 200000 300000 500000; do
              path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.gfa;
              path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.${seq_name}.drop$l.vcf.gz;
              echo $path_vcf_gz;
              vg deconstruct -e -a -p $seq_name -A sampleA $path_gfa -t $threads | bgzip -c > $path_vcf_gz && tabix $path_vcf_gz;
            done;
          done;
        done;
	done;
done
```

Evaluation:

```shell
$run_rtg format -o $path_input_sdf $path_input_fasta

for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	echo $divergence

	name_output=${name_input_fasta}_$divergence

    path_mutated_base_fa=base/$name_output.fa
    cut -f 1 $path_mutated_base_fa.fai | while read seq_name; do
        for aligner in wfmash; do
          mkdir -p variants/$aligner;
          for idx_len in ${!lengths[*]}; do 
            len=${lengths[$idx_len]}; 
            path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.${seq_name}.vcf.gz; 
            for index_s in ${!samples[*]}; do 
              sample=${samples[$index_s]}; 
              samplename=sample$sample; 
              path_truth_mut_vcf_gz=samples/$name_output.$samplename.${seq_name}.vcf.gz;
              echo $aligner $len $sample
              $run_rtg vcfeval -b $path_truth_mut_vcf_gz -c $path_vcf_gz -t $path_input_sdf --sample "$samplename,$samplename" -o evaluations/$aligner/${name_output}_$samplename.l${len}.${seq_name}; 
            done; 
          done; 
        done;

        for aligner in minimap2 winnomap2; do
          mkdir -p variants/$aligner; 
          for idx_len in ${!lengths[*]}; do
            for l in 0 50000 100000 200000 300000 500000; do
              len=${lengths[$idx_len]}; 
              path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.${seq_name}.drop$l.vcf.gz; 
              for index_s in ${!samples[*]}; do
                sample=${samples[$index_s]}; 
                samplename=sample$sample; 
                path_truth_mut_vcf_gz=samples/$name_output.$samplename.${seq_name}.vcf.gz;
                echo $aligner $len $l $sample
                $run_rtg vcfeval -b $path_truth_mut_vcf_gz -c $path_vcf_gz -t $path_input_sdf --sample "$samplename,$samplename" -o evaluations/$aligner/${name_output}_$samplename.l${len}.${seq_name}.drop$l; 
              done;
            done; 
          done; 
        done;
    done
done


grep None evaluations/*/*/summary.txt | grep sampleA | column -t
```
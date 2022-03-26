# Simulations

### Tools

Installation:

```shell
mkdir -p ~/tools
cd ~/tools

git clone --recursive https://github.com/ekg/wfmash.git
cd wfmash
git checkout 25757d3500a84ec1bea64de4fe95782d0b0eb885
cmake -H. -Bbuild && cmake --build build -- -j 48
mv build/bin/wfmash build/bin/wfmash-25757d3500a84ec1bea64de4fe95782d0b0eb885
cd ..

git clone --recursive https://github.com/ekg/seqwish.git
cd seqwish
git checkout 88cd0ea5f086cadfaf21c4c363d71536a1a7ea09
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/seqwish bin/seqwish-88cd0ea5f086cadfaf21c4c363d71536a1a7ea09
cd ..

git clone --recursive https://github.com/AndreaGuarracino/splitfa.git
cd splitfa
git checkout b6545722c429888dbe292fa2cdb265c6ea2c7004
cargo build --release
mv target/release/splitfa target/release/splitfa-b6545722c429888dbe292fa2cdb265c6ea2c7004
cd ..

git clone --recursive https://github.com/natir/fpa.git
cd fpa
git checkout a273badb68c429683369051a1d389ce186f0dd6a
cargo build --release
mv target/release/fpa target/release/fpa-a273badb68c429683369051a1d389ce186f0dd6a
cd ..

git clone --recursive https://github.com/ekg/fastix.git
cd fastix
git checkout 331c1159ea16625ee79d1a82522e800c99206834
cargo build --release
mv target/release/fastix target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834
cd ..

git clone https://github.com/marbl/Winnowmap.git
cd Winnowmap
git checkout b373a7790aed2f5cb68ee63cb5f414ac9d63ec5a
make -j 48
mv bin/ bin-b373a7790aed2f5cb68ee63cb5f414ac9d63ec5a/
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
(echo bcftools vcfrandomsample minimap2 mutation-simulator.py | tr ' ' '\n') | while read tool; do ls -l $(which $tool); done | cut -f 13 -d ' '
    #/gnu/store/nkfg1wg76zqaig43qgslkwcag9rb9fzz-wfmash-0.6.0+e9a5b02-17/bin/wfmash
    /gnu/store/ilksi2bdpsqkyf7r72lkcknndayrpjq3-bcftools-1.12/bin/bcftools
    /gnu/store/nczgf2j577bfw3g84f2b6r774bsmj1a6-vcflib-1.0.2/bin/vcfrandomsample
    /gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
    /gnu/store/k2cs3b88a1ncvalavn85qqkil5h4d03v-mutation-simulator-2.0.3-1.9cb6bd2/bin/mutation-simulator.py

ls -l $(which vg)
  /home/guarracino/tools/vg

/home/guarracino/tools/vg version | head -n 1
  vg version v1.38.0 "Canossa"
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
mv chm13.draft_v1.1.fasta chm13.draft_v1.1.fa
samtools faidx chm13.draft_v1.1.fa

samtools faidx chm13.draft_v1.1.fa chr8 > chm13#chr8.fa
samtools faidx chm13#chr8.fa

cd ..
```

### Paths and variables

Choose the input fasta:

```shell
assembly='CHM13_v1.1'
species='Homo sapiens'

#path_input_fasta=/lizardfs/guarracino/pggb_grant/genomes/chm13.draft_v1.1.fa

path_input_fasta=/lizardfs/guarracino/pggb_grant/genomes/chm13#chr8.fa
```

Set paths:

```shell
# Tools
run_splitfa=~/tools/splitfa/target/release/splitfa-b6545722c429888dbe292fa2cdb265c6ea2c7004
run_mutation_simulator_py=/gnu/store/k2cs3b88a1ncvalavn85qqkil5h4d03v-mutation-simulator-2.0.3-1.9cb6bd2/bin/mutation-simulator.py
run_fastix=~/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834
run_bcftools=/gnu/store/ilksi2bdpsqkyf7r72lkcknndayrpjq3-bcftools-1.12/bin/bcftools
run_vcfrandomsample=/gnu/store/nczgf2j577bfw3g84f2b6r774bsmj1a6-vcflib-1.0.2/bin/vcfrandomsample
run_wfmash=/home/guarracino/tools/wfmash/build/bin/wfmash-a36ab5fa3d435a3030fd584e653b016ca1e89313
run_fpa=~/tools/fpa/target/release/fpa-a273badb68c429683369051a1d389ce186f0dd6a
run_minimap2=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
run_meryl=~/tools/Winnowmap/bin-b373a7790aed2f5cb68ee63cb5f414ac9d63ec5a/meryl
run_winnomap2=~/tools/Winnowmap/bin-b373a7790aed2f5cb68ee63cb5f414ac9d63ec5a/winnowmap
run_seqwish=~/tools/seqwish/bin/seqwish-88cd0ea5f086cadfaf21c4c363d71536a1a7ea09
run_vg=/home/guarracino/tools/vg
run_rtg=/gnu/store/mriq5x6l7kzz51d1z64cvl5qx3d6ylc9-rtg-tools-3.11/rtg

# Variables
threads=48
name_input_fasta=$(basename $path_input_fasta .fa)

divergences=(
  0.001
#  0.01
#  0.05
#  0.15
#  0.20
)

presets=(
  'asm5'
#  'asm10'
#  'asm20'
)

lengths=(
	#200000
	#500000
	#1000000
	5000000
	#10000000
	#20000000
	#50000000
)
  
samples=({A..A})
num_haplotypes=${#samples[@]}

variant_block=10 # Min. distance between simulated variants
path_input_sdf=/lizardfs/guarracino/pggb_grant/${name_input_fasta}.sdf
```

Mutated base generation (SNVs + small INDELs):

```shell
cd /lizardfs/guarracino/pggb_grant/

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
      path_input_fa_gz=/lizardfs/guarracino/pggb_grant/assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz;
      dir_paf=/lizardfs/guarracino/pggb_grant/alignments/minimap2;
      path_name=${name_input_fasta}+samples_$divergence.l${len}.drop0.paf
	  echo $path_input_fa_gz;
	  
	  # Align and filter short alignments
	  sbatch -p workers -c $threads --job-name minimap2 --wrap 'hostname; \time -v '$run_minimap2' '$path_input_fa_gz' '$path_input_fa_gz' -t '$threads' -c -x '$preset' -X 2> '$dir_paf'/'$name_input_fasta'.minimap2.'$divergence'.l'${len}'.log  > '$path_name'; for l in 50000 100000 200000 300000 500000; do cat '$path_name' | '$run_fpa' drop -l $l > '$dir_paf'/'${name_input_fasta}'+samples_'$divergence'.l'${len}'.drop$l.paf; done; mv '$path_name' '$dir_paf;
	done;
done
```

`winnomap2`:

```shell
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
        path_input_fa_gz=/lizardfs/guarracino/pggb_grant/assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz;
        dir_paf=/lizardfs/guarracino/pggb_grant/alignments/winnomap2;
        path_name=${name_input_fasta}+samples_$divergence.l${len}.drop0.paf
        echo $path_input_fa_gz;
        
        # Align and filter short alignments
        sbatch -p workers -c $threads --job-name winnomap2 --wrap 'hostname; cd /scratch; \time -v '$run_winnomap2' -W '${name_input_fasta}'.repetitive_k19.txt '$path_input_fa_gz' '$path_input_fa_gz' -t '$threads' -c -x '$preset' -X 2> '$dir_paf'/'$name_input_fasta'.winnomap2.'$divergence'.l'${len}'.log  > '$path_name'; for l in 50000 100000 200000 300000 500000; do cat '$path_name' | '$run_fpa' drop -l $l > '$dir_paf'/'${name_input_fasta}'+samples_'$divergence'.l'${len}'.drop$l.paf; done; mv '$path_name' '$dir_paf;
	done;
done
```

`wfmash`:

```shell
mkdir -p alignments/wfmash

for index in ${!divergences[*]};
do
	divergence=${divergences[$index]}
	identity=$(echo "(1-$divergence) * 100 - 2" | bc -l)

	name_input_mt_fasta=${name_input_fasta}_$divergence.fa

	for idx_len in ${!lengths[*]}; do len=${lengths[$idx_len]};
        path_input_fa_gz=/lizardfs/guarracino/pggb_grant/assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz;
        dir_paf=/lizardfs/guarracino/pggb_grant/alignments/wfmash;
        
        for s in 5k 10k 20k 50k 100k; do
          for p in 95 98; do
            path_name=${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.paf;
            echo $path_name;
            
            sbatch -p workers -c $threads --job-name wfmash --wrap 'hostname; cd /scratch; \time -v '$run_wfmash' '$path_input_fa_gz' '$path_input_fa_gz' -t '$threads' -s '$s' -p '$p' -X -n '$num_haplotypes' 2> '$dir_paf'/'$name_input_fasta'.wfmash.'$divergence'.l'${len}'.log > '$path_name'; mv '$path_name' '$dir_paf;
          done
        done
	done;
done
```

```shell
# Check
grep 'Elapsed (wall clock)' alignments/*/*log | column -t
grep 'Maximum resident set size (kbytes)' alignments/*/*log | column -t
```

## Graph induction

```shell
mkdir -p /lizardfs/guarracino/pggb_grant/graphs/wfmash/
mkdir -p /lizardfs/guarracino/pggb_grant/graphs/minimap2/
mkdir -p /lizardfs/guarracino/pggb_grant/graphs/winnomap2/

cwd=$(pwd)
cd /scratch

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
	    
	    for s in 5k 10k 20k 50k 100k; do
          for p in 95 98; do
            path_name=${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.paf;
	    
            path_input_paf=alignments/$aligner/${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.paf;
            path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.gfa;
            echo $path_gfa;
            $run_seqwish -t $threads -s /lizardfs/guarracino/pggb_grant/$path_input_fa_gz -p /lizardfs/guarracino/pggb_grant/$path_input_paf -k 0 -g $path_gfa -B 50M -P;
	      done
	    done
	  done;
	done;
    mv graphs/wfmash/* "$cwd"/graphs/wfmash/

	for aligner in minimap2 winnomap2; do
	  mkdir -p graphs/$aligner;
	  for idx_len in ${!lengths[*]}; do
	    len=${lengths[$idx_len]};
	    path_input_fa_gz=assemblies/${name_input_fasta}+samples_$divergence.l${len}.fa.gz;
	    for l in 0 50000 100000 200000 300000 500000; do
	      path_input_paf=alignments/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.paf;
	      path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.gfa;
	      echo $path_gfa;
	      $run_seqwish -t $threads -s /lizardfs/guarracino/pggb_grant/$path_input_fa_gz -p /lizardfs/guarracino/pggb_grant/$path_input_paf -k 0 -g $path_gfa -B 50M -P;
	    done;
	  done;
    done;
    mv graphs/minimap2/* "$cwd"/graphs/minimap2/
    mv graphs/winnomap2/* "$cwd"/graphs/winnomap2/
done

cd $cwd
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
    
    for aligner in wfmash; do
      mkdir -p variants/$aligner;
      for idx_len in ${!lengths[*]}; do
        len=${lengths[$idx_len]};
        
        for s in 5k 10k 20k 50k 100k; do
          for p in 95 98; do
            path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.gfa;
            path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.vcf.gz;
            echo $path_vcf_gz;
            $run_vg deconstruct -e -a -P 'chr' -H '#' $path_gfa -t $threads | bgzip -c > $path_vcf_gz && tabix $path_vcf_gz;
          done
        done
      done
    done;
    
    for aligner in minimap2 winnomap2; do
      mkdir -p variants/$aligner;
      for idx_len in ${!lengths[*]}; do
        len=${lengths[$idx_len]};
        for l in 0 50000 100000 200000 300000 500000; do
          path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.gfa;
          path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.vcf.gz;
          echo $path_vcf_gz;
          $run_vg deconstruct -e -a -P 'chr' -H '#' $path_gfa -t $threads | bgzip -c > $path_vcf_gz && tabix $path_vcf_gz;
        done;
      done;
    done;
    
#    cut -f 1 $path_mutated_base_fa.fai | while read seq_name; do
#        for aligner in wfmash; do
#          mkdir -p variants/$aligner;
#          for idx_len in ${!lengths[*]}; do
#            len=${lengths[$idx_len]};
#            path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.gfa;
#            path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.${seq_name}.vcf.gz;
#            echo $path_vcf_gz $seq_name;
#            vg deconstruct -e -a -p $seq_name -P -H '#' $path_gfa -t $threads | bgzip -c > $path_vcf_gz && tabix $path_vcf_gz;
#            done
#        done;
#        
#        for aligner in minimap2 winnomap2; do
#          mkdir -p variants/$aligner;
#          for idx_len in ${!lengths[*]}; do
#            len=${lengths[$idx_len]};
#            for l in 50000 100000 200000 300000 500000; do
#              path_gfa=graphs/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.gfa;
#              path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.${seq_name}.drop$l.vcf.gz;
#              echo $path_vcf_gz $seq_name;
#              vg deconstruct -e -a -p $seq_name -H '#' $path_gfa -t $threads | bgzip -c > $path_vcf_gz && tabix $path_vcf_gz;
#            done;
#          done;
#        done;
#	done;
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
    
    for aligner in wfmash; do
      mkdir -p variants/$aligner;
      for idx_len in ${!lengths[*]}; do 
        len=${lengths[$idx_len]};
        
        for index_s in ${!samples[*]}; do 
          sample=${samples[$index_s]}; 
          samplename=sample$sample;
          path_truth_mut_vcf_gz=samples/$name_output.$samplename.vcf.gz;
          
          for s in 5k 10k 20k 50k 100k; do
            for p in 95 98; do
              path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.s$s.p$p.vcf.gz;
          
              echo $aligner $len $sample
              $run_rtg vcfeval -b $path_truth_mut_vcf_gz -c $path_vcf_gz -t $path_input_sdf --sample "$samplename,$samplename" -o evaluations/$aligner/${name_output}_$samplename.l${len}.s$s.p$p; 
            done; 
          done;
        done;
      done; 
    done;

    for aligner in minimap2 winnomap2; do
      mkdir -p variants/$aligner; 
      for idx_len in ${!lengths[*]}; do
        for l in 0 50000 100000 200000 300000 500000; do
          len=${lengths[$idx_len]}; 
          path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.drop$l.vcf.gz; 
          for index_s in ${!samples[*]}; do
            sample=${samples[$index_s]}; 
            samplename=sample$sample; 
            path_truth_mut_vcf_gz=samples/$name_output.$samplename.vcf.gz;
            echo $aligner $len $l $sample
            $run_rtg vcfeval -b $path_truth_mut_vcf_gz -c $path_vcf_gz -t $path_input_sdf --sample "$samplename,$samplename" -o evaluations/$aligner/${name_output}_$samplename.l${len}.drop$l; 
          done;
        done; 
      done; 
    done;
    
#    cut -f 1 $path_mutated_base_fa.fai | while read seq_name; do
#        for aligner in wfmash; do
#          mkdir -p variants/$aligner;
#          for idx_len in ${!lengths[*]}; do 
#            len=${lengths[$idx_len]}; 
#            path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.${seq_name}.vcf.gz; 
#            for index_s in ${!samples[*]}; do 
#              sample=${samples[$index_s]}; 
#              samplename=sample$sample;
#              path_truth_mut_vcf_gz=samples/$name_output.$samplename.${seq_name}.vcf.gz;
#              echo $aligner $len $sample
#              $run_rtg vcfeval -b $path_truth_mut_vcf_gz -c $path_vcf_gz -t $path_input_sdf --sample "$samplename,$samplename" -o evaluations/$aligner/${name_output}_$samplename.l${len}.${seq_name}; 
#            done; 
#          done; 
#        done;
#
#        for aligner in minimap2 winnomap2; do
#          mkdir -p variants/$aligner; 
#          for idx_len in ${!lengths[*]}; do
#            for l in 0 50000 100000 200000 300000 500000; do
#              len=${lengths[$idx_len]}; 
#              path_vcf_gz=variants/$aligner/${name_input_fasta}+samples_$divergence.l${len}.${seq_name}.drop$l.vcf.gz; 
#              for index_s in ${!samples[*]}; do
#                sample=${samples[$index_s]}; 
#                samplename=sample$sample; 
#                path_truth_mut_vcf_gz=samples/$name_output.$samplename.${seq_name}.vcf.gz;
#                echo $aligner $len $l $sample
#                $run_rtg vcfeval -b $path_truth_mut_vcf_gz -c $path_vcf_gz -t $path_input_sdf --sample "$samplename,$samplename" -o evaluations/$aligner/${name_output}_$samplename.l${len}.${seq_name}.drop$l; 
#              done;
#            done; 
#          done; 
#        done;
#    done
done


grep None evaluations/*/*/summary.txt | grep sampleA | column -t
```
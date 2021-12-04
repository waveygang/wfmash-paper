# Simulations

### Installation

```shell
mkdir -p ~/tools
cd ~/tools

# Mutation-Simulator
git clone --recursive https://github.com/mkpython3/Mutation-Simulator.git
cd Mutation-Simulator
git checkout 9cb6bd2acf8201151bc610be14963e65b41d8899
cd ..

# splitfa
git clone --recursive https://github.com/AndreaGuarracino/splitfa.git
cd splitfa
git checkout b6545722c429888dbe292fa2cdb265c6ea2c7004
cargo build --release
cd ..

# On octopus01
git clone --recursive https://gitlab.com/genenetwork/guix-bioinformatics.git

# minimap2
cd guix-bioinformatics
GUIX_PACKAGE_PATH=. guix install minimap2

# rtg vcfeval
GUIX_PACKAGE_PATH=. guix install rtg-tools

# Mutation-Simulator dependencies)
GUIX_PACKAGE_PATH=. guix install python-tqdm python-pandas python-numpy python-pyfaidx
```

### Versions

```shell
(echo wfmash minimap2 | tr ' ' '\n') | while read tool; do ls -l $(which $tool); done | cut -f 13 -d ' '

/gnu/store/nkfg1wg76zqaig43qgslkwcag9rb9fzz-wfmash-0.6.0+e9a5b02-17/bin/wfmash
/gnu/store/1nhgd86l7r2j4r5f420jhyd65kvc52a6-minimap2-2.18/bin/minimap2
```

Create the main folder:

```shell
mkdir -p /lizardfs/guarracino/pggb_grant/
cd /lizardfs/guarracino/pggb_grant/
```

### Obtain the data

Download and prepare the references:

```shell
mkdir -p /lizardfs/guarracino/pggb_grant/genomes
cd /lizardfs/guarracino/pggb_grant/genomes

wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz
gunzip chm13.draft_v1.1.fasta.gz
samtools faidx chm13.draft_v1.1.fasta

cd ..
```

### Paths

```shell
# Tools
run_splitfa=~/tools/splitfa/target/release/splitfa
path_mutation_simulator_py=~/tools/Mutation-Simulator/mutation-simulator.py
run_wfmash=/gnu/store/nkfg1wg76zqaig43qgslkwcag9rb9fzz-wfmash-0.6.0+e9a5b02-17/bin/wfmash
run_minimap2=/gnu/store/1nhgd86l7r2j4r5f420jhyd65kvc52a6-minimap2-2.18/bin/minimap2
run_rtg=/gnu/store/mriq5x6l7kzz51d1z64cvl5qx3d6ylc9-rtg-tools-3.11/rtg

# Input
assembly='CHM13_v1.1'
species='Homo sapiens'
path_input_fasta=/lizardfs/guarracino/pggb_grant/genomes/chm13.draft_v1.1.fasta
path_input_sdf=/lizardfs/guarracino/pggb_grant/chm13.draft_v1.1.sdf
```

```shell
# Variables

name_input_fasta=$(basename $path_input_fasta .fasta)
len_input_fasta=$(cut $path_input_fasta.fai -f2)

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
	#10000000
	#20000000
	#50000000
)
```

```shell
# Mutated chromosome generation (SNVs + small INDELs)
cd /lizardfs/guarracino/pggb_grant/

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
	python3 $path_mutation_simulator_py $path_input_fasta --output $prefix_mutations_vcf_pre args --snp $snv_rate -titv 2 \
	  --insert $ins_rate --snpblock $variant_block --insertlength 5 --insertblock $variant_block \
	  --deletion $del_rate --deletionlength 5 --deletionblock $variant_block \
	  --assembly $assembly --species "$species" \
	  --sample ${name_input_fasta}_mt

    # Remove invalid variants (REF == ALT, in case on Ns in the reference)
	name_output=${name_input_fasta}_$divergence
	path_mutations_vcf=chromosomes/$name_output.vcf

    cat <(grep '^#' $prefix_mutations_vcf_pre.vcf) <(grep '^#' $prefix_mutations_vcf_pre.vcf -v | awk '{if ($4 != $5){print$0}}') > $path_mutations_vcf

    rm $prefix_mutations_vcf_pre.vcf
    mv chromosomes/$name_output_pre.fa chromosomes/$name_output.fa
done
```
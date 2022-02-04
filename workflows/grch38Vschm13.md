# GRCh38 vs CHM13

## Tools

### Octopus

```shell
#(echo wfmash | tr ' ' '\n') | while read tool; do ls -l $(which $tool); done | cut -f 13 -d ' '

mkdir -p ~/tools $$ cd ~/tools

git clone --recursive https://github.com/ekg/wfmash.git
cd wfmash

git checkout 09e73eb3fcf24b8b7312b8890dd0741933f0d1cd
cmake -H. -Bbuild && cmake --build build -- -j 48
mv build/bin/wfmash build/bin/wfmash-09e73eb3fcf24b8b7312b8890dd0741933f0d1cd

git checkout 948f1683d14927745aef781cdabeb66ac6c7880b
cmake -H. -Bbuild && cmake --build build -- -j 48
mv build/bin/wfmash build/bin/wfmash-948f1683d14927745aef781cdabeb66ac6c7880b

cd ..

git clone --recursive https://github.com/mrvollger/rustybam.git
cd rustybam/
git checkout 63a8ab437f4f04c95aa5c0f6683f3171d47aefd6
cargo build --release
```

### BSC

```shell
module load intel mkl gsl jemalloc htslib cmake gcc/10.2.0
LIBRARY_PATH=$LIBRARY_PATH:/apps/JEMALLOC/5.2.1/INTEL/lib

# Copy wfmash repository
scp -r wfmash bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/
cd wfmash
git checkout 6f4a9248f34470dfe36196e96e018fe1f526a8c8
cmake -H. -Bbuild && cmake --build build -- -j 128
```


### Obtain the data

Create the main folder:

```shell
mkdir -p /lizardfs/guarracino/vgp/grch38_vs_chm13
cd /lizardfs/guarracino/vgp/grch38_vs_chm13
```

Download and prepare the references:

```shell
mkdir -p /lizardfs/guarracino/vgp/grch38_vs_chm13/genomes
cd /lizardfs/guarracino/vgp/grch38_vs_chm13/genomes
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz
TO DO: rename sequences to obtain: chm13.fa and grch38.fa

samtools faidx chm13.fa
samtools faidx grch38.fa
```

Download and prepare the annotation:

```shell
mkdir -p /lizardfs/guarracino/vgp/grch38_vs_chm13/annotation
cd /lizardfs/guarracino/vgp/grch38_vs_chm13/annotation

wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
zgrep '^#' -v  gencode.v38.annotation.gtf.gz | awk -v OFS='\t' '$3 == "gene" {print $1,$4,$5,$10,".",$7}' | sed 's/"//g' | sed 's/;//g' > gencode.v38.annotation.genes.bed
```


# TODO: masking simulated/modeled regions
 - Modeled centromeres and heterochromatin regions
Modeled centromeres and heterochromatin regions (from https://www.ncbi.nlm.nih.gov/grc/human): https://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/38/Modeled_regions_for_GRCh38.tsv

 - Regions (with false duplications) that are getting masked to improve the alignment on chr21: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_GRC_exclusions.bed

 - Info (but broken link to a BED): http://genomeref.blogspot.com/2021/07/one-of-these-things-doest-belong.html also contains a BED file that contains the
extra falsely duplicated regions on 21p, however the link to the bed file seems broken


### Mapping evaluation

Map. We use `-c < num_of_threads` to allow SLURM to run multiple processes on the same node.
Indeed, the mapping is limited by the number of chromosomes.

```shell
cd /lizardfs/guarracino/vgp/grch38_vs_chm13/
mkdir /lizardfs/guarracino/vgp/grch38_vs_chm13/mappings

query=chm13.fa
target=grch38.fa

for p in 99 98 95 90 85 80 75; do
    for s in 2M 1M 900k 800k 700k 600k 500k 450k 400k 350k 300k 250k 200k 150k 100k 50k 20k 10k; do
        l=0
        for n in 1 5 10 20 50 100 200 500 1000; do
            for w in 0; do
                prefix=$query-$target.s$s.l$l.p$p.n$n.w$w
        
                if [ ! -f mappings/$prefix.approx.paf ]; then
                    sbatch -p 386mem -c 12 --wrap '\time -v ~/tools/wfmash/build/bin/wfmash-09e73eb3fcf24b8b7312b8890dd0741933f0d1cd genomes/'${target}' genomes/'${query}' -t 12 -s '$s' -l '$l' -p '$p' -n '$n' -w '$w' -m >mappings/'$prefix'.approx.paf'
                fi
            done
        done
    done
done
```

Evaluate the mappings:

```shell
path_gencode_genes_target=annotation/gencode.v38.annotation.genes.bed

# Number of GENCODE genes in the target
total_genes_in_target=$(cat $path_gencode_genes_target | wc -l)

echo query target total_genes_in_target s l p n w present_gene_ratio_target missing_genes_in_target | tr ' ' '\t' >gencode_evaluation.mapping.tsv
for p in 99 98 95 90 85 80 75; do
    for s in 2M 1M 900k 800k 700k 600k 500k 450k 400k 350k 300k 250k 200k 150k 100k 50k 20k 10k; do
        l=0
        for n in 1 5 10 20 50 100 200 500 1000; do
            for w in 0; do
                prefix=$query-$target.s$s.l$l.p$p.n$n.w$w
                echo $prefix
                
                if [ -f mappings/$prefix.approx.paf ]; then
                    #cat mappings/$prefix.approx.paf | awk -v OFS='\t' '{print $1, $3, $4, "", "", $5}' >mappings/$prefix.approx.$query.bed
                    cat mappings/$prefix.approx.paf | awk -v OFS='\t' '{print $6, $8, $9, "", "", "+"}' | sed 's/grch38#//g' >mappings/$prefix.approx.$target.bed
            
                    # Number of GENCODE genes not entirely covered
                    missing_genes_in_target=$(bedtools subtract -a $path_gencode_genes_target -b mappings/$prefix.approx.$target.bed | cut -f 4 | sort | uniq | wc -l)
            
                    missing_gene_ratio_target=$(echo "scale=4; 1 - $missing_genes_in_target / $total_genes_in_target" | bc)
            
                    echo $query $target $total_genes_in_target $s $l $p $n $w $missing_gene_ratio_target $missing_genes_in_target | tr ' ' '\t' >>gencode_evaluation.mapping.tsv
                fi
            done
      done
    done
done
```

## Alignment evaluation

```shell
cd /lizardfs/guarracino/vgp/grch38_vs_chm13/
mkdir /lizardfs/guarracino/vgp/grch38_vs_chm13/alignment

query=chm13.fa
target=grch38.fa

s=500k
l=0
p=75
n=200
w=0

prefix=$query-$target.s$s.l$l.p$p.n$n.w$w
echo $prefix

\time -v ~/tools/wfmash/build/bin/wfmash-948f1683d14927745aef781cdabeb66ac6c7880b genomes/${target} genomes/${query} -t 48 -s $s -l $l -p $p -n $n -w $w -i mappings/$prefix.approx.paf > alignment/$prefix.paf
#\time -v ~/tools/wfmash/build/bin/wfmash-948f1683d14927745aef781cdabeb66ac6c7880b genomes/'${target}' genomes/'${query}' -t 48 -s '$s' -l '$l' -p '$p' -n '$n' -w '$w' -i mappings/'$prefix'.approx.paf > alignment/'$prefix'.paf


path_gencode_genes_target=annotation/gencode.v38.annotation.genes.bed

# Number of GENCODE genes in the target
total_genes_in_target=$(cat $path_gencode_genes_target | wc -l)

#cat alignment/$prefix.paf | awk -v OFS='\t' '{print $1, $3, $4, "", "", $5}' >alignment/$prefix.$query.bed
cat alignment/$prefix.paf | awk -v OFS='\t' '{print $6, $8, $9, "", "", "+"}' | sed 's/grch38#//g' >alignment/$prefix.$target.bed

# Number of GENCODE genes not entirely covered
missing_genes_in_target=$(bedtools subtract -a $path_gencode_genes_target -b alignment/$prefix.$target.bed | cut -f 4 | sort | uniq | wc -l)

missing_gene_ratio_target=$(echo "scale=4; 1 - $missing_genes_in_target / $total_genes_in_target" | bc)

echo $query $target $total_genes_in_target $s $l $p $n $w $missing_gene_ratio_target $missing_genes_in_target | tr ' ' '\t' >>gencode_evaluation.mapping.tsv
```

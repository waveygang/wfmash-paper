# GRCh38 vs CHM13

### Octopus

```
(echo wfmash | tr ' ' '\n') | while read tool; do ls -l $(which $tool); done | cut -f 13 -d ' '

/gnu/store/nkfg1wg76zqaig43qgslkwcag9rb9fzz-wfmash-0.6.0+e9a5b02-17/bin/wfmash
```

### BSC

```
module load intel mkl gsl jemalloc htslib cmake gcc/10.2.0
LIBRARY_PATH=$LIBRARY_PATH:/apps/JEMALLOC/5.2.1/INTEL/lib

# Copy wfmash repository
scp -r wfmash bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/
cd wfmash
git checkout 6f4a9248f34470dfe36196e96e018fe1f526a8c8
cmake -H. -Bbuild && cmake --build build -- -j 128
```

Create the main folder:

```
mkdir -p /lizardfs/guarracino/vgp/grch38_vs_chm13
cd /lizardfs/guarracino/vgp/grch38_vs_chm13
```

### Obtain the data

Download and prepare the references:

```
mkdir -p /lizardfs/guarracino/vgp/grch38_vs_chm13/genomes
cd /lizardfs/guarracino/vgp/grch38_vs_chm13/genomes
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz
TO DO: rename sequences to obtain: chm13.fa and grch38.fa

samtools faidx chm13.fa
samtools faidx grch38.fa
```

Download and prepare the annotation:

```
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
#TODO: clean the gene name in the 4-th column of the BED ($10)
zgrep '^#' -v  gencode.v38.annotation.gtf.gz | awk -v OFS='\t' '$3 == "gene" {print $1,$4,$5,$10,".",$7}' > gencode.v38.annotation.genes.bed

```

### Mapping evaluation

TODO: USE sbatch -c 24 and wfmash -t 24, to allow SLURM to run multiple processes on the same node. The mapping is
limited by the number of chromosomes TODO: -p 75 and 70

Map:

```
mkdir /lizardfs/guarracino/vgp/grch38_vs_chm13/mappings

query=chm13.fa
target=grch38.fa

for p in 99 98 95 90 85 80; do
    for s in 2M 1M 900k 800k 700k 600k 500k 450k 400k 350k 300k 250k 200k 150k 100k 50k 20k 10k; do
        l=0
        for n in 1 5 10 20 50 100 200 500 1000; do
            for w in 0; do
                prefix=$query-$target.s$s.l$l.p$p.n$n.w$w
        
                if [ ! -f mappings/$prefix.approx.paf ]; then
                    sbatch -p lowmem -c 24 --wrap '\time -v /gnu/store/nkfg1wg76zqaig43qgslkwcag9rb9fzz-wfmash-0.6.0+e9a5b02-17/bin/wfmash genomes/'${target}' genomes/'${query}' -t 24 -s '$s' -l '$l' -p '$p' -n '$n' -w '$w' -m >mappings/'$prefix'.approx.paf'
                fi
            done
        done
    done
done
```

Evaluate the mappings:

```
path_gencode_genes_target=gencode.v38.annotation.genes.bed

# Number of GENCODE genes in the target
total_genes_in_target=$(cat $path_gencode_genes_target | wc -l)

echo query target total_genes_in_target s l p n w present_gene_ratio_target missing_genes_in_target | tr ' ' '\t' >gencode_evaluation.mapping.tsv
for p in 99 98 95 90 85 80; do
    for s in 2M 1M 900k 800k 700k 600k 500k 450k 400k 350k 300k 250k 200k 150k 100k 50k 20k 10k; do
        l=0
        for n in 1 5 10 20 50 100 200 500 1000; do
            for w in 0; do
                prefix=$query-$target.s$s.l$l.p$p.n$n.w$w
        
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


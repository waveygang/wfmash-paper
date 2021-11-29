# 45 fish

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

#### Octopus

```
mkdir -p /lizardfs/guarracino/vgp/45_fish
cd /lizardfs/guarracino/vgp/45_fish
```

#### BSC

```
mkdir -p /gpfs/projects/bsc18/bsc18995/vgp/45_fish
cd /gpfs/projects/bsc18/bsc18995/vgp/45_fish
```

### Obtain and explore the genomes

#### Octopus

Put in `/lizardfs/guarracino/vgp/45_fish` the `45_fish_alignment.fixed.xlsx` file. In the original
file (`45_fish_alignment.xlsx`) there were swapped IDs.

Put in `/lizardfs/guarracino/vgp/45_fish` the `45_fish_alignment.download.py` file.

```
python3 45_fish_alignment.download.py 45_fish_alignment.fixed.xlsx genomes

# Decompress genomes
cd /lizardfs/guarracino/vgp/45_fish/genomes
ls /lizardfs/guarracino/vgp/45_fish/genomes/*.fna.gz | while read f; do gunzip $f; done
```

#### BSC (where there is no internet connection)

```
mkdir -p /gpfs/projects/bsc18/bsc18995/vgp/45_fish/genomes

#Log in on Octopus
scp -r /lizardfs/guarracino/vgp/45_fish/genomes bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/vgp/45_fish/genomes
```

Compute the mash distances:

#### Octopus

```
# guix install mash

cd /lizardfs/guarracino/vgp/45_fish/genomes
ls *.fna | while read f; do mash sketch $f; done
mash triangle *.fna >45_fish_alignment.mash_triangle.txt
```

Check the nucleotide composition:

#### Octopus

```
# https://github.com/lh3/seqtk/issues/47
# guix install seqtk

ls *.fna | while read f; do seqtk comp $f; done
```

Use the `mash_triangle_heatmap.R` to produce a heatmap for the `45_fish_alignment.mash_triangle.txt` file

### Obtain BUSCO genes

```
mkdir /lizardfs/guarracino/vgp/45_fish/busco_genes/
cd /lizardfs/guarracino/vgp/45_fish/busco_genes/
```

Prepare the `vertebrata_odb10` database:

#### Octopus

```
wget -c https://busco-data.ezlab.org/v5/data/lineages/vertebrata_odb10.2021-02-19.tar.gz
tar -xvzf vertebrata_odb10.2021-02-19.tar.gz && rm vertebrata_odb10.2021-02-19.tar.gz
```

#### BSC (where there is no internet connection)

```
mkdir -p /gpfs/projects/bsc18/bsc18995/45_fish/busco_genes

#Log in on Octopus
wget -c https://busco-data.ezlab.org/v5/data/lineages/vertebrata_odb10.2021-02-19.tar.gz
scp vertebrata_odb10.2021-02-19.tar.gz bsc18995@amdlogin.bsc.es:/gpfs/projects/bsc18/bsc18995/45_fish/busco_genes
rm vertebrata_odb10.2021-02-19.tar.gz

#Log in on BSC
cd /gpfs/projects/bsc18/bsc18995/45_fish/busco_genes
tar -xvzf vertebrata_odb10.2021-02-19.tar.gz && rm vertebrata_odb10.2021-02-19.tar.gz
```

Identify BUSCO genes in each genome:

#### Octopus

```
#git clone --recursive https://gitlab.com/genenetwork/guix-bioinformatics.git
#cd guix-bioinformatics
#GUIX_PACKAGE_PATH=. guix install busco

sbatch -p lowmem -c 48 --wrap 'cd /lizardfs/guarracino/vgp/45_fish/busco_genes/ && ls /lizardfs/guarracino/vgp/45_fish/genomes/*.fna | while read path_genome; do name_genome="$(basename "$path_genome")"; if [[ ! -s /lizardfs/guarracino/vgp/45_fish/busco_genes/output_$name_genome/run_vertebrata_odb10/short_summary.txt ]]; then busco -i $path_genome --lineage_dataset /lizardfs/guarracino/vgp/45_fish/busco_genes/vertebrata_odb10 -o output_$name_genome --mode genome --cpu 48 --offline --metaeuk_parameters="--remove-tmp-files=1" --metaeuk_rerun_parameters="--remove-tmp-files=1"; fi; done'

# `busco` doesn't work on /scratch
#ls /lizardfs/guarracino/vgp/45_fish/genomes/*.fna | while read path_genome; do name_genome="$(basename "$path_genome")"; sbatch -p lowmem -c 48 --wrap 'cd /scratch && busco -i '"$path_genome"' --lineage_dataset /lizardfs/guarracino/vgp/45_fish/busco_genes/vertebrata_odb10 -o '"$name_genome"' --mode genome --cpu 48 --offline --augustus && mv '"$name_genome"' /lizardfs/guarracino/vgp/45_fish/busco_genes/ '; done

# `augustus_specie` doesn't work
busco -i /lizardfs/guarracino/vgp/45_fish/genomes/GCA_007364275.2_fArcCen1_genomic.fna --lineage_dataset /lizardfs/guarracino/vgp/45_fish/busco_genes/vertebrata_odb10 -o xxx --mode genome --cpu 48 --offline --augustus_specie fish
```

#### BSC

```
module load anaconda
source activate busco # version 5.2.2

dir_genomes=/gpfs/projects/bsc18/bsc18995/vgp/45_fish/genomes
dir_busco_genes=/gpfs/projects/bsc18/bsc18995/vgp/45_fish/busco_genes
path_vertebrata_odb10_db=/gpfs/projects/bsc18/bsc18995/45_fish/busco_genes/vertebrata_odb10/

ls $dir_genomes/*.fna | while read path_genome; do
    name_genome="$(basename "$path_genome")";
    dir_busco_genes_current=$dir_busco_genes/output_$name_genome;
    if [[ ! -s $dir_busco_genes_current/run_vertebrata_odb10/short_summary.txt ]]; then
        sbatch -c 128 --wrap 'TMPFOLDER=/scratch/tmp/$SLURM_JOBID; cd $TMPFOLDER; busco -i '$path_genome' --lineage_dataset '$path_vertebrata_odb10_db' -o output_'$name_genome' --mode genome --cpu 128 --offline --metaeuk_parameters="--remove-tmp-files=1" --metaeuk_rerun_parameters="--remove-tmp-files=1"; mv output_$name_genome '$dir_busco_genes;
    fi;
done
```

# ToDo Cleaning and directory renaming

Prepare BED files for the completely identified BUSCO genes:

# ToDo on Octopus

```
cd /lizardfs/guarracino/vgp/45_fish/busco_genes/
ls */full_table.tsv | while read f; do
dir_parent=$(dirname "$f")
grep Complete "$f" | awk -v OFS='\t' '{print $3, $4, $5, $1, $7, $6}' >"$dir_parent"/$dir_parent.busco_genes.complete.bed
done
cd /lizardfs/guarracino/vgp/45_fish
```

### Mapping evaluation


```
mkdir /lizardfs/guarracino/vgp/45_fish/mappings

query=GCA_904848185.1_fAcaLat1.1
target=GCA_900880675.2_fSpaAur1.2

path_complete_busco_genes_query=busco_genes/"$query"/$query.busco_genes.complete.bed
path_complete_busco_genes_target=busco_genes/"$target"/$target.busco_genes.complete.bed

# Number of complete BUSCO genes in the genomes
total_genes_in_query=$(cat $path_complete_busco_genes_query | wc -l)
total_genes_in_target=$(cat $path_complete_busco_genes_target | wc -l)

echo query target total_genes_in_query total_genes_in_target s l p w present_gene_ratio_query present_gene_ratio_target missing_genes_in_query missing_genes_in_target | tr ' ' '\t' >busco_evaluation.mapping.tsv
for s in 20k 50k 100k 150k 200k 250k 300k 350k 400k 450k 500k 600k 700k 800k 900k 1M; do
  l=$s
  for p in 85 80 75 70; do
    for w in 256 0; do
      prefix=$query-$target.s$s.l$l.p$p.w$w

      if [ ! -f $prefix.approx.paf ]; then
        wfmash genomes/${target}_genomic.fna genomes/${query}_genomic.fna -t 48 -s $s -l $l -p $p -w $w -m >mappings/$prefix.approx.paf
      fi

      cat mappings/$prefix.approx.paf | awk -v OFS='\t' '{print $1, $3, $4, "", "", $5}' >mappings/$prefix.approx.$query.bed
      cat mappings/$prefix.approx.paf | awk -v OFS='\t' '{print $6, $8, $9, "", "", "+"}' >mappings/$prefix.approx.$target.bed

      # Number of complete BUSCO genes not entirely covered
      missing_genes_in_query=$(bedtools subtract -a $path_complete_busco_genes_query -b mappings/$prefix.approx.$query.bed | cut -f 4 | sort | uniq | wc -l)
      missing_genes_in_target=$(bedtools subtract -a $path_complete_busco_genes_target -b mappings/$prefix.approx.$target.bed | cut -f 4 | sort | uniq | wc -l)

      missing_gene_ratio_query=$(echo "scale=4; 1 - $missing_genes_in_query / $total_genes_in_query" | bc)
      missing_gene_ratio_target=$(echo "scale=4; 1 - $missing_genes_in_target / $total_genes_in_target" | bc)
      #ToDo compute the average?

      echo $query $target $total_genes_in_query $total_genes_in_target $s $l $p $w $missing_gene_ratio_query $missing_gene_ratio_target $missing_genes_in_query $missing_genes_in_target | tr ' ' '\t' >>busco_evaluation.mapping.tsv
    done
  done
done
```

# ???

```

sbatch -p 1tbmem -c 48 --wrap 'cd /scratch && \time -v /home/guarracino/wfmash/build/bin/wfmash /lizardfs/guarracino/vgp/45_fish/genomes/GCA_900880675.2_fSpaAur1.2_genomic.fna.gz /lizardfs/guarracino/vgp/45_fish/genomes/GCA_904848185.1_fAcaLat1.1_genomic.fna.gz -t 48 -s 50k -l 50k -p 70 -i /lizardfs/guarracino/vgp/45_fish/GCA_904848185.1_fAcaLat1.1-GCA_900880675.2_fSpaAur1.2.s50k.l50k.p70.approx.paf > GCA_904848185.1_fAcaLat1.1-GCA_900880675.2_fSpaAur1.2.s50k.l50k.p70.paf && mv GCA_904848185.1_fAcaLat1.1-GCA_900880675.2_fSpaAur1.2.s50k.l50k.p70.paf /lizardfs/guarracino/vgp/45_fish/'

seq 0 30 | while read i; do sbatch -p lowmem -c 30 --wrap 'cd /scratch && \time -v /home/guarracino/wfmash/build/bin/wfmash /lizardfs/guarracino/vgp/45_fish/genomes/GCA_900880675.2_fSpaAur1.2_genomic.fna.gz /lizardfs/guarracino/vgp/45_fish/genomes/GCA_904848185.1_fAcaLat1.1_genomic.fna.gz -t 30 -s 50k -l 50k -p 70 -i /lizardfs/guarracino/vgp/45_fish/GCA_904848185.1_fAcaLat1.1-GCA_900880675.2_fSpaAur1.2.s50k.l50k.p70.approx.paf.chunk_'$i' > GCA_904848185.1_fAcaLat1.1-GCA_900880675.2_fSpaAur1.2.s50k.l50k.p70.chunk_'$i'.paf'; done
```
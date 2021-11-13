num_threads=48

mkdir -p /lizardfs/guarracino/vgp/45_fish
cd /lizardfs/guarracino/vgp/45_fish

# Put in /lizardfs/guarracino/vgp/45_fish the 45_fish_alignment.fixed.xlsx file
# In the original file (45_fish_alignment.xlsx) there were swapped IDs

# Put in /lizardfs/guarracino/vgp/45_fish the 45_fish_alignment.download.py file
python3 45_fish_alignment.download.py 45_fish_alignment.fixed.xlsx genomes

# Compute mash distances
# guix install mash
cd /lizardfs/guarracino/vgp/45_fish/genomes
ls *.fna.gz | while read f; do mash sketch $f; done
mash triangle *.fna.gz >45_fish_alignment.mash_triangle.txt

# Use the mash_triangle_heatmap.R to produce a heatmap for the 45_fish_alignment.mash_triangle.txt file

# Check nucleotide composition
# https://github.com/lh3/seqtk/issues/47
# guix install seqtk
ls *.fna.gz | while read f; do seqtk comp $f; done


# Get BUSCO genes
mkdir /lizardfs/guarracino/vgp/45_fish/busco_genes/
# ToDo on Octopus
# busco -i GCA_900880675.2_fSpaAur1.2_genomic.fna --lineage_dataset vertebrata_odb10 -o output_GCA_900880675.2_fSpaAur1.2 --mode genome --cpu 16
# ToDo Cleaning and directory renaming


# Prepare BED files for the completely identified BUSCO genes
# ToDo on Octopus
cd /lizardfs/guarracino/vgp/45_fish/busco_genes/
ls */full_table.tsv | while read f; do
  dir_parent=$(dirname "$f")
  grep Complete "$f" | awk -v OFS='\t' '{print $3, $4, $5, $1, $7, $6}' >"$dir_parent"/$dir_parent.busco_genes.complete.bed
done
cd /lizardfs/guarracino/vgp/45_fish


# Mapping and evaluate
mkdir /lizardfs/guarracino/vgp/45_fish/mappings

query=GCA_904848185.1_fAcaLat1.1
target=GCA_900880675.2_fSpaAur1.2

path_complete_busco_genes_query=busco_genes/"$query"/$query.busco_genes.complete.bed
path_complete_busco_genes_target=busco_genes/"$target"/$target.busco_genes.complete.bed

# Number of complete BUSCO genes in the genomes
total_genes_in_query=$(cat $path_complete_busco_genes_query | wc -l)
total_genes_in_target=$(cat $path_complete_busco_genes_target | wc -l)

echo query target total_genes_in_query total_genes_in_target s l p w present_gene_ratio_query present_gene_ratio_target missing_genes_in_query missing_genes_in_target | tr ' ' '\t' > busco_evaluation.mapping.tsv
for s in 20k 50k 100k 150k 200k 250k 300k 350k 400k 450k 500k 600k 700k 800k 900k 1M
do
  l=$s
  for p in 85 80 75 70
  do
    for w in 256 0
    do
      prefix=$query-$target.s$s.l$l.p$p.w$w

      if [ ! -f $prefix.approx.paf ]; then
        wfmash genomes/${target}_genomic.fna.gz genomes/${query}_genomic.fna.gz -t $num_threads -s $s -l $l -p $p -w $w -m >mappings/$prefix.approx.paf
      fi

      cat mappings/$prefix.approx.paf | awk -v OFS='\t' '{print $1, $3, $4, "", "", $5}' >mappings/$prefix.approx.$query.bed
      cat mappings/$prefix.approx.paf | awk -v OFS='\t' '{print $6, $8, $9, "", "", "+"}' >mappings/$prefix.approx.$target.bed

      # Number of complete BUSCO genes not entirely covered
      missing_genes_in_query=$(bedtools subtract -a $path_complete_busco_genes_query -b mappings/$prefix.approx.$query.bed | cut -f 4 | sort | uniq | wc -l)
      missing_genes_in_target=$(bedtools subtract -a $path_complete_busco_genes_target -b mappings/$prefix.approx.$target.bed | cut -f 4 | sort | uniq | wc -l)

      missing_gene_ratio_query=$(echo "scale=4; 1 - $missing_genes_in_query / $total_genes_in_query" | bc)
      missing_gene_ratio_target=$(echo "scale=4; 1 - $missing_genes_in_target / $total_genes_in_target" | bc)

      echo $query $target $total_genes_in_query $total_genes_in_target $s $l $p $w $missing_gene_ratio_query $missing_gene_ratio_target $missing_genes_in_query $missing_genes_in_target | tr ' ' '\t' >> busco_evaluation.mapping.tsv
    done
  done
done

num_threads=16

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


# Prepare BED files for the BUSCO genes identified completely in the genomes
# ToDo on Octopus
cd /lizardfs/guarracino/vgp/45_fish/busco_genes/
ls */full_table.tsv | while read f; do
  dir_parent=$(dirname "$f")
  grep Complete "$f" | awk -v OFS='\t' '{print $3, $4, $5, $1, $7, $6}' >"$dir_parent"/$dir_parent.busco_genes.complete.bed
done
cd /lizardfs/guarracino/vgp/45_fish


# Align and evaluate
query=GCF_904848185.1_fAcaLat1.1
target=GCA_900880675.2_fSpaAur1.2
s=500k
l=500k
p=80

prefix=$query-$target.s300k.l300k.p80
path_complete_busco_genes_query=busco_genes/"$query"/$query.busco_genes.complete.bed
path_complete_busco_genes_target=busco_genes/"$target"/$target.busco_genes.complete.bed

wfmash ${target}_genomic.fna ${query}_genomic.fna -t $num_threads -s $s -l $l -p $p -m >$prefix.approx.paf
cat $prefix.approx.paf | awk -v OFS='\t' '{print $1, $3, $4, "", "", "+"}' >$prefix.approx.$query.bed
cat $prefix.approx.paf | awk -v OFS='\t' '{print $6, $8, $9, "", "", $5}' >$prefix.approx.$target.bed

# Number of complete BUSCO genes in the genomes
total_genes_in_query=$(cat $path_complete_busco_genes_query | wc -l)
total_genes_in_target=$(cat $path_complete_busco_genes_target | wc -l)

# Number of complete BUSCO genes not entirely covered
missing_genes_in_query=$(bedtools subtract -a $path_complete_busco_genes_query -b $prefix.approx.$query.bed | cut -f 4 | sort | uniq | wc -l)
missing_genes_in_target=$(bedtools subtract -a $path_complete_busco_genes_target -b $prefix.approx.$target.bed | cut -f 4 | sort | uniq | wc -l)

echo "scale=4; 1 - $missing_genes_in_query / $total_genes_in_query" | bc
echo "scale=4; 1 - $missing_genes_in_target / $total_genes_in_target" | bc

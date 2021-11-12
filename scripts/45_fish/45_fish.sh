
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
mash triangle *.fna.gz > 45_fish_alignment.mash_triangle.txt

# Use the mash_triangle_heatmap.R to produce an heatmap for the 45_fish_alignment.mash_triangle.txt



# https://github.com/lh3/seqtk/issues/47
# guix install seqtk
ls *.fna.gz | while read f; do seqtk comp $f; done

# Read in a mash output lower triangle file as a matrix object
read_mash_tri <- function(file) {
  x <- scan(file, what = 'character', skip = 1)   # read in as vector
  dims <- floor(sqrt(length(x) * 2))              # get dimensions
  m <- matrix(NA, dims, dims)                     # construct matrix
  m[upper.tri(m, diag = TRUE)] <- x               # fill in values from vector to upper tri
  m <- t(m)                                       # transpose to get lower tri (better way?)
  rownames(m) <- m[, 1]                            # add rownames from first col
  m <- m[, -1]                                     # remove column containing rownames
  extracol <- c(rep(NA, nrow(m) - 1), '0')          # construct vector for missing column
  m <- cbind(m, extracol)                         # bind missing column
  diag(m) <- 0                                    # set diagonal to zero
  colnames(m) <- rownames(m)                      # add column names
  m <- as.dist(m)                                 # coerce to dist type object
  return(m)                                       # return dist object
}

D <- read_mash_tri(file = '45_fish_alignment.mash_triangle.txt')


# Rename rows and columns
D <- as.matrix(D)
rownames(D) <- gsub("_genomic.*", "", rownames(D))
colnames(D) <- gsub("_genomic.*", "", colnames(D))
D <- as.dist(D)


library(ade4) #install.packages("ape")
temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg = 0, clabel.row = .9, clabel.col = .9) #darker shades of gray mean a larger distance # you can also make cool color plots but they're much more complicated because they use the image() function


pheatmap::pheatmap(
  #cor(as.matrix(D), method = "spearman"),
  as.matrix(D),
  fontsize = 12, fontsize_row = 12, height = 20,
  #annotation = neoant_meta_and_data_subset %>%
  #  dplyr::select(SampleName, `Cell line`, Ploidy) %>%
  #  column_to_rownames(var = "SampleName"),
  #annotation_colors=list(
  #  Ploidy = c(Hyperdiploid = '#F8766D', `Near-to-diploid`='#00BFC4'),
  #  `Cell line` = c(`17` = '#D89000', `30` = '#DA62DB', `32` = '#55CC7D')
)
# 45_fish_alignment.mash_triangle.heatmap.11x11.5.pdf
# 45_fish_alignment.mash_triangle.heatmap.1550x1500.png


# Other trials
library(ape) #install.packages("ape") #https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html
tre <- nj(D)
class(tre) #all trees created using {ape} package will be of class phylo

tre <- ladderize(trw)
tre # tells us what the tree will look like but doesn't show the actual construction
plot(tre, cex = 0.6)
title("45_fish_alignment.mash_triangle")

# or

h_cluster <- hclust(D, method = "average", members = NULL) # method = average is used for UPGMA, members can be equal to NULL or a vector with a length of size D
plot(h_cluster, cex = 0.6)

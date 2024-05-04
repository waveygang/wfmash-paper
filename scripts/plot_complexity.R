library(ggplot2)
library(dplyr)

options(scipen = 9)

# Read the BED files
path_sattools_bed_chm13v2 <- '/home/guarracino/Desktop/Garrison/wfmash-paper/chm13v2.sattols.bed'
path_sattools_bed_hg002v101 <- '/home/guarracino/Desktop/Garrison/wfmash-paper/hg002v101.sattols.bed'

# Read and process chm13v2 data
data_chm13v2 <- read.table(path_sattools_bed_chm13v2, header = FALSE, sep = "\t",
                           col.names = c("chrom", "start", "end", "seq.complexity"))
data_chm13v2$position <- (data_chm13v2$start + data_chm13v2$end) / 2
data_chm13v2$source <- "chm13v2"
data_chm13v2$chrom_number <- data_chm13v2$chrom
data_chm13v2$parent <- NA

# Read and process hg002v101 data
data_hg002v101 <- read.table(path_sattools_bed_hg002v101, header = FALSE, sep = "\t", comment.char = '?',
                             col.names = c("chrom", "start", "end", "seq.complexity"))
data_hg002v101$position <- (data_hg002v101$start + data_hg002v101$end) / 2
data_hg002v101 <- data_hg002v101 %>%
  mutate(chrom_number = ifelse(grepl("chrEBV|chrM$", chrom), 
                               sub(".*#(chrEBV|chrM)$", "\\1", chrom), 
                               sub(".*#.*#(.*?)_.*", "\\1", chrom)),
         parent = ifelse(grepl("#P#", chrom), "PATERNAL", ifelse(grepl("#M#", chrom), "MATERNAL", NA)),
         source = paste0("hg002v101.", tolower(parent))
         )
data_hg002v101 <- data_hg002v101 %>%
  mutate(
    parent = ifelse(grepl("#P#", chrom), "PATERNAL", ifelse(grepl("#M#", chrom), "MATERNAL", NA)),
    source = paste0("hg002v101.", tolower(parent))
  )

# Merge the two datasets
data <- rbind(data_chm13v2, data_hg002v101)

# Define the desired order of chromosomes for faceting
chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                 "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
                 "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM", "chrEBV")

# Convert 'chrom_number' to a factor with the desired order
data$chrom_number <- factor(data$chrom_number, levels = chrom_order)

# Filter out 'chrM'
data <- data %>%
  filter(!chrom_number %in% c('chrM', 'chrEBV'))

scale <- 1000000
x_axis_label <- 'Position (Mbp)'

# Create the ggplot
p <- ggplot(data, aes(x = position / scale, y = seq.complexity, color = source)) +
  geom_line(alpha = 0.7) +  # Set alpha for transparency
  facet_wrap(~ chrom_number, ncol = 1, scales = "fixed") +
  labs(
    x = paste(x_axis_label),
    y = "Sequence Complexity",
    title = "Sequence complexity across chromosomes",
    color = "Source"
  ) +
  ylim(0, 1.0) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

height <- max(8, length(unique(data$chrom_number)) * 5)
path_output <- 'merged_plot.png'
width <- 70

ggsave(plot = p, path_output, width = width, height = height, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)
































library(ggplot2)
library(dplyr)

options(scipen = 9)

# Read the BED file
path_sattools_bed <- '/home/guarracino/Desktop/Garrison/wfmash-paper/chm13v2.sattols.bed'
path_sattools_bed <- '/home/guarracino/Desktop/Garrison/wfmash-paper/hg002v101.sattols.bed'
data <- read.table(path_sattools_bed, header = FALSE, sep = "\t", comment.char = '?',
                   col.names = c("chrom", "start", "end", "seq.complexity"))

# Create a new column 'position' as the midpoint of start and end
data$position <- (data$start + data$end) / 2

# Extract the chromosome number and parent from the 'chrom' column
data <- data %>%
  mutate(chrom_number = sub(".*#(chr\\d+)_.*", "\\1", chrom),
         parent = sub(".*_(MATERNAL|PATERNAL)", "\\1", chrom))
# Filter the data based on the specified parent
parent_to_plot <- "PATERNAL"
data <- data %>%
  filter(parent == parent_to_plot)
data$chrom <- data$chrom_numer


# Define the desired order of chromosomes for faceting
chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
                 "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
                 "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")

# Convert 'chrom' to a factor with the desired order
data$chrom <- factor(data$chrom, levels = chrom_order)

data <- data %>%
  filter(! chrom %in% c('chrM'))

scale <- 1000000
x_axis_label <- 'Position (Mbp)'

# Create the ggplot
p <- ggplot(data, aes(x = position / scale, y = seq.complexity)) +
  geom_line() +
  facet_wrap(~ chrom, ncol = 1, scales = "fixed") +
  labs(
    x = paste(x_axis_label),
    y = "Sequence Complexity",
    title = "Sequence complexity across chromosomes") +
  ylim(0, 1.0) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))
height <- max(8, length(unique(data$chrom))*5)
path_output <- 'hg002v101.paternal.png'
width <- 60
ggsave(plot = p, path_output, width = width, height = height, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)


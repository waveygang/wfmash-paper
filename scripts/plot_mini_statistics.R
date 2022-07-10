# Input from:
#  grep None *-penalties/chr*/*/*/summary.txt  | cut -f 1,2,3,4,5 -d '/' | tr -s " " | sed 's/summary.txt: //g' | tr '/' ' ' | cut -f 1,2,3,4,8,9,10,11,12 -d ' ' | tr ' ' '\t' > mini.tsv 

library(ggplot2)
library(ggrepel)

x <- read.table('/home/guarracino/mini.tsv', sep = '\t', header = F)
names(x) <- c('branch', 'chromosome', 'sample', 'region', 'FP', 'FN', 'precision', 'sensitivity', 'F1.score')
x$chromosome <- factor(x$chromosome, levels = unique(x[order(as.integer(gsub("[^0-9]", "", x$chromosome))), ]$chromosome))

a <- x[(x$chromosome != 'chrX' & x$chromosome != 'chrY'),]
ggplot(a, aes(x = chromosome, y = F1.score, fill=branch, label = sample)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(size=0.4, alpha=0.9, aes(color=branch)) +
  #ylim(0.98, 0.996) +
  #geom_text_repel(
  #  size=2,
  #  max.iter=100,
  #  max.time=2,
  #  show.legend  = FALSE, # to hide the `a` from the legend
  #  max.overlaps=Inf
  #) +
  facet_wrap (~region, ncol = 1, scales = "free") + 
  theme(legend.position = "top")


xy <- x[(x$chromosome == 'chrX' | x$chromosome == 'chrY'),]
ggplot(xy, aes(x = chromosome, y = F1.score, fill=branch, label = sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  geom_text_repel(
    size=3,
    max.iter=100,
    max.time=2,
    show.legend  = FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  ) +
  facet_wrap (~region, ncol = 2, scales = "free") + 
  theme(legend.position = "top")

library(tidyr)
library(dplyr)
library(ggplot2)

x <- read.table('/home/guarracino/mini_tomato.stats.tsv', header = T, sep = '\t', comment.char = '?')
x$chromosome <- factor(x$chromosome, levels = unique(x[order(as.integer(gsub("[^0-9]", "", x$chromosome))), ]$chromosome))

################################################################################
# Sort branches as needed
x$branch <- factor(x$branch, levels = c('fixed-0-4-6-1', 'fixed-0-7-11-1', 'fixed-0-11-17-1'))
################################################################################

xx <- x %>%
  dplyr::select(-tp.baseline,-tp.call,-fp,-fn) %>%
  pivot_longer(., precision:f1.score,"metric")
xx$metric <- factor(xx$metric, levels = c('precision', 'recall', 'f1.score'))

# If there are just 2 branches
if (FALSE) {
  xxx <- xx %>% 
    group_by(chromosome,sample,vcf,metric) %>% 
    mutate(delta=diff(value))
  
  ggplot(xxx %>% filter(vcf == 'vg'), aes(x = chromosome, y = delta, fill=metric, label = sample)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size=1, alpha=0.9, color='black') +
    ylim(-max(abs(min(xxx$delta)), abs(max(xxx$delta))), +max(abs(min(xxx$delta)), abs(max(xxx$delta)))) +
    #ggrepel::geom_text_repel(
    #  size=2,
    #  max.iter=10,
    #  max.time=2,
    #  show.legend=FALSE, # to hide the `a` from the legend
    #  max.overlaps=Inf
    #) +
    facet_grid(~metric, scales = "free") + 
    theme(legend.position = "top") +
    theme_bw()
  
}


ggplot(xx %>% filter(vcf == 'vg'), aes(x = branch, y = value, fill=branch, label = sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=1, alpha=0.9, aes(color=branch)) +
  ylim(min(xx$value), 1) +
  #ggrepel::geom_text_repel(
  #  size=2,
  #  max.iter=10,
  #  max.time=2,
  #  show.legend=FALSE, # to hide the `a` from the legend
  #  max.overlaps=Inf
  #) +
  facet_grid(metric~chromosome, scales = "free") + 
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )

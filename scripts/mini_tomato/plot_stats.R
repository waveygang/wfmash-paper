x <- read.table('mini_tomato.stats.tsv', header = T, sep = '\t', comment.char = '?')
x$chromosome <- factor(x$chromosome, levels = unique(x[order(as.integer(gsub("[^0-9]", "", x$chromosome))), ]$chromosome))

xx <- x %>%
  select(-tp.baseline,-tp.call,-fp,-fn) %>%
  tidyr::pivot_longer(., precision:f1.score,"metric")
xx$metric <- factor(xx$metric, levels = c('precision', 'recall', 'f1.score'))

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

ggplot(xx %>% filter(vcf == 'vg'), aes(x = branch, y = value, fill=branch, label = sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=1, alpha=0.9, aes(color=branch)) +
  ylim(min(y$value), 1) +
  ggrepel::geom_text_repel(
    size=2,
    max.iter=10,
    max.time=2,
    show.legend=FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  ) +
  facet_grid(metric~chromosome, scales = "free") + 
  theme(legend.position = "top") +
  theme_bw()

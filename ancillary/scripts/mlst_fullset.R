library(tidyverse)
library(RColorBrewer)


setwd('~/Projects/bactocap/datasets/mlst/')
mapstats <- read_delim('idxstats.txt', delim = '\t')
targstats <- read_delim('targstats.txt', delim = '\t')
metadata <- read_delim('mlst_metadata.tsv', delim = '\t')

#remove garbage from end of sample name
mapstats$sample_id <- gsub("\\_.*","",mapstats$sample_id)

#join mapping stats to metadata
#inflates the metadata to contain one entry w/metadata per-target per-sample
dfa <- metadata %>% left_join(mapstats, by = (c('IDInBactocap' = 'sample_id')))
dfb <- dfa %>% left_join(targstats, by = c('locus' = 'locus_name'))
#drop useless columns
dfb <- dfb %>% select(-c(X1,`base_pl byte_pl`,byte_index, PolyomicsList, BATCH))

#add per-sample frac mapped and frac unmapped
dfc <- dfb %>% 
  group_by(IDInBactocap) %>% 
  summarise(total_mapped = sum(mapped_reads), total_reads = sum(mapped_reads, unmapped_reads), unmapped_total = sum(unmapped_reads)) %>% 
  mutate(frac_unmapped = unmapped_total/total_reads) %>% 
  mutate(frac_mapped = total_mapped/total_reads)

#split locus name to create a column containing the target species
dfb$targ_species <- gsub("^.*_", "", dfb$locus)
#calculate per-target doc (no.mapped reads * 75 (rd length) / target length)
dfb$targ_mean_doc <- (dfb$mapped_reads*75)/dfb$sequence_length


palette = brewer.pal(n = 5, name = "YlGnBu")
dfb %>% select(IDInBactocap, mapped_reads, targ_species, sample_composition) %>% 
  filter(!(targ_species == '*')) %>% 
  ggplot(., aes(fill=targ_species, y=mapped_reads, x=IDInBactocap)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid( ~ sample_composition, scales = "free", space='free') +
  scale_fill_manual(values=palette)

library(tidyverse)

cov <- read.csv('~/Projects/bactocap/datasets/mlst/analysis/mlst-mlstonly-cov.csv', col.names = c("sample_name", "target_name", "target_len", "reads_mapped", "unmapped_frags"))
cov <- cov %>% extract(target_name, into = c("sample_pref", "target_org"), "(.*)_([^_]+)$")

#get mean d.o.c over lepto targets and summarise as table
cov %>% 
  filter(target_org == 'Lepto') %>% 
  mutate(meandoc = (reads_mapped * 75)/target_len) %>% select(sample_name, meandoc) %>% 
  ggplot(aes(x=sample_name, y = meandoc)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal()

#return per-sample mapped vs unmapped reads
cov %>% 
  select(sample_name, reads_mapped, unmapped_frags) %>%
  group_by(sample_name) %>% 
  summarise(sum_map = sum(reads_mapped), sum_unmap = sum(unmapped_frags)) %>% 
  mutate(Fraction_Of_Reads_Mapped = sum_map / (sum_map + sum_unmap))
  
  
#create lepto/notlepto target observation, then aggregate by whether reads mapped to a given target are mapped to a lepto target or not
#"the number of on vs off-target reads

tib <- cov %>% mutate(.,
               target = case_when(
                 target_org == 'Lepto' ~ 'On-target',
                 target_org != 'Lepto' ~ 'Off-target',
                 TRUE ~ NA_character_)) %>% 
  select(sample_name, target, reads_mapped) %>% 
  group_by(sample_name, target) %>% 
  summarise(reads = sum(reads_mapped)) 

tab <- tib %>% pivot_wider(names_from = target, values_from = reads) %>% 
  mutate(Fraction_On_Target = `On-target`/(`On-target`+`Off-target`)) %>% 
  select(sample_name, `Off-target`, `On-target`, Fraction_On_Target)


tob <- cov %>% filter(target_org == 'Lepto') %>% 
  mutate(meandoc = (reads_mapped * 75)/target_len) %>% select(sample_name, sample_pref, meandoc) %>% 
  pivot_wider(names_from = sample_name, values_from = meandoc)



#now, get d.o.c over each target per-sample and filter based on the lepto only ones
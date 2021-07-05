#ls pkg
pkg = c("tidyverse", "cowplot", "RColorBrewer", "cowplot", "pheatmap")
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#



setwd('~/Projects/bactocap/datasets/mlst/')
mapstats <- read_delim('/Users/tristanpwdennis/Projects/bactocap/ancillary/metadata/idxstats_udref.txt', delim = '\t', col_names = FALSE)
colnames(mapstats) = c('sample_id', 'locus', 'length', 'mapped', 'unmapped')
mapstats$sample_id = sub("/[^/]+$", "", mapstats$sample_id) 
mapstats$sample_id = sub("_[^_]+$", "", mapstats$sample_id) 
#fix mapstats sample id
targstats <- read_delim('~/Projects/bactocap/datasets/mlst/metadata/target_stats.txt', delim = '\t')
metadata <- read_delim('~/Projects/bactocap/datasets/mlst/metadata/mlst_metadata.tsv', delim = '\t')
mlstmapping <- read_csv('~/Projects/bactocap/datasets/mlst/metadata/mapstatstable_mlst.csv')



dir <- c("~/Projects/bactocap/datasets/mlst/results/")
new_tbl <- NULL 
files <- list.files(path = dir,pattern = ".sample_interval_summary", full.names = T)
for (file in files) {
  the_data <- read_csv(file)
  the_data<-the_data %>% mutate(filename = file)
  the_data <- the_data %>% select(1,2,3,9,10)
  colnames(the_data) <- c("target", "total", "average", "percent_above_15","filename")
  new_tbl<-rbind(new_tbl,the_data)
}

#clean out the filename
new_tbl$filename <- str_replace(new_tbl$filename, "/Users/tristanpwdennis/Projects/bactocap/datasets/mlst/results//", "")
new_tbl$filename <- str_replace(new_tbl$filename, ".sample_interval_summary", "")


s= mapstats %>% select(sample_id) %>% group_by(sample_id) %>% count()
s= mapstats %>% filter(sample_id == 'Bacto-16')



mapstatstable = mapstats %>% dplyr::select(sample_id, locus, mapped) %>% pivot_wider(names_from = 'sample_id', values_from = 'mapped')
write_csv(mapstatstable, "/Users/tristanpwdennis/Projects/bactocap/datasets/mlst/results/mapstatstable_mlst.csv")


mapstats %>% filter(X1 > 60)

?pivot_wider
#separate the target column into organism, locus, type etc
 
new_tbl$target<-sub(":", "_", new_tbl$target)
new_tbl$target<-sub("_[^_]+$", "", new_tbl$target)
new_tbl$targetname <- new_tbl$target
new_tbl <- new_tbl %>% extract(target, into = c("target", "target_organism"), "(.*)_([^_]+)$") 
new_tbl <- new_tbl %>% extract(target, into = c("locus", "locus_type"), "(.*)_([^_]+)$")
#remove polyomics suffix
new_tbl$filename <- sub("_[^_]+$", "", new_tbl$filename)


#select sample composition from metadata and join to new_tbl
ud_tbl <- metadata %>% select(IDInBactocap, sample_composition) %>% left_join(., new_tbl, by = c("IDInBactocap" = "filename"))


target_annotations <- ud_tbl %>% select(targetname, target_organism) %>% distinct() %>% remove_rownames %>% column_to_rownames(var="targetname")
sample_annotations <- ud_tbl %>% select(IDInBactocap, sample_composition) %>% distinct() %>% remove_rownames %>% column_to_rownames(var="IDInBactocap")
heatmapdata <- ud_tbl %>% select(targetname, IDInBactocap, percent_above_15, target_organism)
heatmapdata <- heatmapdata %>% pivot_wider(names_from=IDInBactocap, values_from = percent_above_15)
heatmapdata <- heatmapdata %>% arrange(., target_organism)
heatmapdata <- heatmapdata %>% remove_rownames %>% column_to_rownames(var="targetname")
heatmapdata <- heatmapdata %>% select(-target_organism)

#make nice heatmap
pheatmap(heatmapdata,cluster_rows=FALSE, cluster_cols=FALSE, annotation_col = sample_annotations, annotation_row = target_annotations)



#join mapping stats to metadata
#inflates the metadata to contain one entry w/metadata per-target per-sample
dfb <- metadata %>% 
  left_join(mapstats, by = (c('IDInBactocap' = 'sample_id'))) %>% 
  left_join(targstats, by = c('locus' = 'locus_name')) %>% 
  select(-c(X1,PolyomicsList, BATCH))

#add per-sample frac mapped and frac unmapped
dfc <- dfb %>% 
  group_by(IDInBactocap) %>% 
  summarise(total_mapped = sum(mapped_reads), total_reads = sum(mapped_reads, unmapped_reads), unmapped_total = sum(unmapped_reads)) %>% 
  mutate(frac_unmapped = unmapped_total/total_reads) %>% 
  mutate(frac_mapped = total_mapped/total_reads) %>% 
  left_join(., dfb)
dfc$frac_mapped_to_targ = dfc$mapped_reads / dfc$total_mapped
#split locus name to create a column containing the target species
dfc$targ_species <- gsub("^.*_", "", dfb$locus)
#calculate per-target doc (no.mapped reads * 75 (rd length) / target length)
dfc$targ_mean_doc <- (dfb$mapped_reads*75)/dfb$sequence_length
dfd <- dfc %>% mutate(on_or_off_targ = case_when(
  targ_species == 'Lepto' & sample_composition == 'lepto' ~ 'on_target',
  targ_species == 'Coxiella' & sample_composition == 'cox' ~ 'on_target',
  targ_species == 'Bartonella' & sample_composition == 'bart' ~ 'on_target',
  targ_species == 'Brucella' & sample_composition == 'bru' ~ 'on_target',
  sample_composition == 'bru_cox' & (targ_species == 'Brucella' | targ_species == 'Coxiella') ~ 'on_target',
  sample_composition == 'bru_cox' & (targ_species == 'Brucella' | targ_species == 'Coxiella') ~ 'on_target',
  sample_composition == 'cox_bart' & (targ_species == 'Bartonella' | targ_species == 'Coxiella') ~ 'on_target',
  TRUE ~ as.character('off_target')
)) 


#calculate number of 'on target' targets for each sample
dfe <- dfd %>% filter(on_or_off_targ == 'on_target') %>% 
  select(locus, IDInBactocap) %>% group_by(IDInBactocap) %>% tally() %>% rename(num_targs = n) %>% 
  left_join(., dfd)

#define function for calculating the number of targets with mean doc > 1 the threshold defined
gettargscov <- function(df, covthreshold) {
  df %>% filter(on_or_off_targ == 'on_target') %>% 
    filter(field_control == 'field' )%>% 
    mutate(mean_doc = (mapped_reads*75)/sequence_length) %>% 
    filter(mean_doc > covthreshold) %>% 
    group_by(IDInBactocap, num_targs, sample_composition) %>% tally() %>% mutate(doc = covthreshold)
}


#create stepwise d.o.c sequence
x <- seq(1,10000,10)

#calculate decreasing d.o.c over all targs for each sample
outputdf <- NULL
for (num in x) {
  s <- gettargscov(dfe, num)
  outputdf <- rbind(outputdf, s)
}


#plot of single vs double cap - don't really know how to interpret this
single_double_line <- dfd %>% filter(field_control == 'control') %>%
  group_by(`Sample type`, Capture) %>% 
  summarise(frac_mapped = sum(mapped_reads)/(sum(mapped_reads, unmapped_reads))) %>% 
  select(`Sample type`, Capture, frac_mapped) %>% 
  pivot_wider(names_from = 'Capture', values_from = 'frac_mapped') %>% 
  ggplot(aes(x=Double, y=Single)) +
  geom_point() +
  geom_abline(slope =1, color = '#e31a1c') +
  theme_minimal()
#dotted lined boxplot  of single vs double cap - don't really know how to interpret this
single_double_boxplot <- dfd %>% filter(field_control == 'control') %>%
  group_by(`Sample type`, Capture) %>% 
  summarise(frac_mapped = sum(mapped_reads)/(sum(mapped_reads, unmapped_reads))) %>% 
  ggplot(aes(x=Capture, y = frac_mapped, fill =Capture)) +
  geom_boxplot(width = 0.3) +
  ylim(0,1) +
  geom_point() +
  geom_line(aes(group = `Sample type`), alpha = 0.5) +
  scale_fill_manual(values = c("#577590", "#F94144")) +
  theme_minimal() +
  ylab("Proportion of Reads Mapped")
single_double_boxplot
#remove garbage
t<-dfd %>% select(IDInBactocap, frac_mapped, sample_composition, Capture, `Host species`) %>% distinct()
#models
mx <- m0 <- MASS::glm.nb(data=t, frac_mapped)
m0 <- MASS::glm.nb(data=t, frac_mapped ~Capture)
m1 <- MASS::glm.nb(data=t, frac_mapped ~Capture*sample_composition)
m2 <- MASS::glm.nb(data=t, frac_mapped ~Capture*`Host species`)
summary(m0)
summary(m1)
summary(m2)
#plot models
interactions::cat_plot(m0, pred = Capture, interval = TRUE, plot.points = TRUE, line.thickness = 0.5)


dilution <- dfd %>% filter(!(is.na(Dilution_factor))) %>% 
  group_by(`Sample type`, Dilution_factor, Diluant, Capture)%>% 
  summarise(frac_mapped = sum(mapped_reads)/(sum(mapped_reads, unmapped_reads))) %>% 
  ggplot(aes(x=Dilution_factor, y = frac_mapped, color = Diluant)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Capure) +
  theme_minimal()+
  scale_fill_manual("legend", values = c("#F3722C", "#43AA8B"))+
  ylab("Proportion of Reads Mapped")





#get proportion mapped to unmapped for the control samples
proportion_mapped_control <- dfd %>% group_by(IDInBactocap) %>% 
  filter(field_control == 'control') %>% 
  filter(sample_composition != 'bart') %>% 
  filter(!(grepl("Pos Control Bru Serial Dilution", SelectionRationale))) %>% 
  group_by(IDInBactocap, sample_composition) %>% 
  summarise(unmapped_sum = sum(unmapped_reads), mapped_sum = sum(mapped_reads)) %>% 
  pivot_longer(cols = c('unmapped_sum', 'mapped_sum'), names_to = 'mapped_or_unmapped', values_to = 'reads') %>% 
  ggplot(aes(x=IDInBactocap, y=reads, fill = mapped_or_unmapped)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+  
  facet_grid( ~ sample_composition, scales = "free", space='free') +
  scale_fill_manual("legend", values = c("#F3722C", "#577590"))
proportion_mapped_control

#get proportion mapped to unmapped for the control samples
proportion_mapped_field <- dfd %>% group_by(IDInBactocap) %>% 
  filter(field_control == 'field') %>% 
  group_by(IDInBactocap, sample_composition) %>% 
  summarise(unmapped_sum = sum(unmapped_reads), mapped_sum = sum(mapped_reads)) %>% 
  pivot_longer(cols = c('unmapped_sum', 'mapped_sum'), names_to = 'mapped_or_unmapped', values_to = 'reads') %>% 
  ggplot(aes(x=IDInBactocap, y=reads, fill = mapped_or_unmapped)) +
  geom_bar(position="fill", stat="identity") +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
  strip.background = element_blank(),
strip.text.x = element_blank())+
  facet_grid( ~ sample_composition, scales = "free", space='free')+
  ylab("Proportion_Of_Total_Reads_Mapped")

t<-dfd %>% group_by(IDInBactocap) %>% 
  filter(field_control == 'control') %>% 
  group_by(IDInBactocap, sample_composition) %>% 
  summarise(unmapped_sum = sum(unmapped_reads), mapped_sum = sum(mapped_reads)) %>% 
  mutate(proportion_mapped = mapped_sum/(mapped_sum+unmapped_sum))

#plot stacked bar plot of composition of all the samples, for control samples
u <- dfd %>% select(IDInBactocap, mapped_reads, targ_species, sample_composition, field_control) %>% 
  filter(!(targ_species == '*')) %>% 
  filter(field_control == 'field') %>% 
  ggplot(., aes(fill=targ_species, y=mapped_reads, x=IDInBactocap)) + 
  geom_bar(position="fill", stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+  
  facet_grid( ~ sample_composition, scales = "free", space='free') 

#plot above for field samples
sample_composition_stacked <- dfc %>% select(IDInBactocap, mapped_reads, targ_species, sample_composition, field_control) %>% 
  filter(!(targ_species == '*')) %>% 
  filter(field_control == 'field') %>% 
  ggplot(., aes(fill=targ_species, y=mapped_reads, x=IDInBactocap)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+  
  facet_grid( ~ sample_composition, scales = "free", space='free') +
  ylab("Proportion_Reads_Mapped_To_Organism_Targets") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#plot breakdown of hte above, faceted by organism
sampel_composition_windowed <- dfc %>% select(IDInBactocap, frac_mapped_to_targ, targ_species, sample_composition) %>% 
  filter(!(targ_species == '*')) %>% 
  ggplot(., aes(fill=targ_species, y=frac_mapped_to_targ, x=IDInBactocap)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap( ~ targ_species, ncol = 1, nrow = 4) 

#samples by % on target mapping
on_target_mapping_boxplot <- dfe %>%
  filter(locus != '*') %>%  
  group_by(IDInBactocap, on_or_off_targ, sample_composition, field_control) %>% 
  summarise(readonorofftarg = sum(mapped_reads)/total_mapped) %>% 
  filter(on_or_off_targ == 'on_target') %>% distinct() %>% 
  ggplot(aes(x=sample_composition, y = readonorofftarg, colour= sample_composition))+
  geom_boxplot(outlier.colour = NA, alpha = 0.5, width = 0.5) +
  geom_jitter() +
  theme_minimal()+
  facet_wrap(~field_control)+
  theme(legend.position = "none") +
  ylab('% Reads On-Target') +
  xlab('PCR +ve organism in sample')

on_target_mapping_boxplot
#fraction of reads mapped per-sample
frac_per_sample_map_boxplot <- dfe %>% 
  distinct(frac_mapped, sample_composition, field_control) %>% 
  ggplot(aes(x=sample_composition, y = frac_mapped, colour = sample_composition)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5, width = 0.5) +
  geom_jitter()+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  facet_wrap(~field_control)+
  ylim(0.0,1)+
  ylab('% Of Total Reads Mapped') +
  xlab('PCR +ve organism in sample')

legend <- get_legend(
  # create some space to the left of the legend
  frac_per_sample_map_boxplot + theme(legend.box.margin = margin(0, 0, 0, 12))
)

#fraction of total readset mapped
dfc %>% select(IDInBactocap, frac_mapped, sample_composition) %>% distinct() %>%  
  ggplot(aes(x = IDInBactocap, y = frac_mapped)) +
  geom_jitter() +
  theme_minimal() +
  facet_grid( ~ sample_composition, scales = "free", space='free') +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('% Of Total Reads Mapped') +
  xlab('ID')


###########################################################################
#Plotting coverage
dfc %>% 
  filter(targ_mean_doc > 100) %>% 
  ggplot(aes(x = IDInBactocap, y = targ_mean_doc)) +
  geom_jitter() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid( ~ sample_composition, scales = "free", space='free') 

#create cowplots

#for proportion mapping based on single and double cap
toprowcaptures <- plot_grid(single_double_boxplot, labels = 'C')
lab_plots <- plot_grid(proportion_mapped_control, dilution, toprowcaptures, labels = c('A', 'B', ''), label_size = 12, ncol = 1)
lab_plots

#plot composition plots
plotslist <- list(proportion_mapped_field, sample_composition_stacked)
n_x = sapply(plotslist, function(p) {
  max(layer_data(p, 1)$x)
})
plot_grid(proportion_mapped_field, sample_composition_stacked, align = 'v', ncol = 1, rel_widths = n_x)
proportion_mapped_field
#plot dec cov plots
bottomrow <- plot_grid(frac_mapped, frac_mapped_lim, align = 'h',labels = c('C', 'D'))
toprow <- plot_grid(frac_per_sample_map_boxplot,on_target_mapping_boxplot, align = 'h', labels = c('A', 'B'))

plot_grid(frac_per_sample_map_boxplot, on_target_mapping_boxplot, align = 'v', ncol = 1, axis = 'b', labels = c('A', 'B'))


proportion_mapped_control
single_double_boxplot
frac_mapped
bottomrow <- plot_grid(single_double_boxplot, frac_mapped, align = 'h')
plot_grid(proportion_mapped_control, bottomrow, align = 'v', ncol = 1)
single_double_boxplot


legend

plot_grid(single_double_boxplot, frac_per_sample_map_boxplot, dilution, NULL, ncol = 2, align = 'h', axis = 'b', labels = c("A", "B", "C"))

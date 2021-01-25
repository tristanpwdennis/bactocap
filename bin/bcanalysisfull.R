library(tidyverse)
library(data.table)
library(sjPlot)
library(cowplot)

#func to read csvs and take filename as a column
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = basename(flnm))
}


#################
#Per-sample fraction of baits covered to a given d.o.c

anthrax_dir <- c("~/Projects/bactocap/datasets/anthrax/results/")
myco_dir <- c("~/Projects/bactocap/datasets/mycoplasma/results/")
mlst_dir <- c("~/Projects/bactocap/datasets/mlst/results/")

dirs <- c(anthrax_dir, myco_dir)

#read csvs from list of files

new_tbl <- NULL 
for (dir in dirs){
tbl <- list.files(path = dir,pattern = ".sample_interval_statistics", full.names = T) %>% 
  map_df(~read_plus(.)) 
tbl$organism = paste0(dir)
new_tbl <- rbind(new_tbl, tbl)
}

#replace gatk's gubbins from the header
colnames(new_tbl) <- str_replace(colnames(new_tbl), "depth>=", "")
#remove useless column (Number of Sources)
new_tbl <- new_tbl %>% select(-Number_of_sources)
new_tbl$organism

new_tbl$filename <- str_replace(new_tbl$filename, ".sample_interval_statistics", "")
new_tbl$organism <- str_replace(new_tbl$organism, "~/Projects/bactocap/datasets/", "")
new_tbl$organism <- str_replace(new_tbl$organism, "/results/", "")


#pivot into longer
new_tbl <- new_tbl %>% pivot_longer(cols = c(0:501))
new_tbl$frac <- fifelse(new_tbl$organism == 'mycoplasma', new_tbl$value/24444, new_tbl$value/148729)
new_tbl$name <- as.numeric(new_tbl$name)
new_tbl <- new_tbl %>% select(-filename, -organism)
new_tbl$name <- as.numeric(new_tbl$name)
new_tbl %>% ggplot(aes(x=name, y=value)) +
  geom_line()+
  theme_minimal() +
  facet_wrap(~sample_composition) +
  theme(legend.position = "none") 


##############################
#persample summary data

sum_tbl <- NULL 
for (dir in dirs){
  tbl <- list.files(path = dir,pattern = ".sample_summary", full.names = T) %>% 
    map_df(~read_plus(.)) %>% filter(sample_id != "Total")
  tbl$organism = paste0(dir)
  sum_tbl <- rbind(sum_tbl, tbl)
}


sum_tbl$filename <- str_replace(sum_tbl$filename, ".sample_summary", "")
sum_tbl$organism <- str_replace(sum_tbl$organism, "~/Projects/bactocap/datasets/", "")
sum_tbl$organism <- str_replace(sum_tbl$organism, "/results/", "")

#plot mean_d_o_c
sum_tbl %>% 
  ggplot(aes(x = organism, y=mean, fill=organism)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("#577590", "#F94144")) +
  theme_minimal() +
  ylab("Mean Depth-of-Coverage") +
  ylab("Organism")

metadata <- rbind(read.csv("~/Projects/bactocap/metadata/anthrax-metadata.csv") %>% mutate(organism = 'anthrax'),
read.csv("~/Projects/bactocap/metadata/mycoplasma-metadata.csv") %>% mutate(organism = 'mycoplasma'))

mappingdata  <- rbind(read.csv("~/Projects/bactocap/datasets/mycoplasma/results/myco_mapping.csv") %>% mutate(organism = 'mycoplasma'), 
read.csv("~/Projects/bactocap/datasets/anthrax/results/anth_mapping.csv") %>% mutate(organism = 'anthrax'))



total_tbl

#################
#let's take a look at our metadata
tb1 <- left_join(sum_tbl, metadata, by = c('sample_id'='polyom_id'))
total_tbl <- left_join(tb1, mappingdata, by = c('sample_id'='sample_id'))
total_tbl <- total_tbl %>% mutate(frac_mapped = mapped/total.y)
total_tbl$sample_conc <- as.numeric(total_tbl$sample_conc)
#################
#frac  mapped vs ct value
total_tbl %>% 
  ggplot(aes(x=max_ct, y=frac_mapped, colour=organism)) +
  geom_point()+
  scale_colour_manual(values = c( "#F94144", "#577590")) +
  #geom_smooth(method=lm) +
  theme_minimal() +
  labs(y='Fraction of Reads Mapped', x='Max. Ct Value')
#################
#total  mapped vs ct value
total_tbl %>% 
  ggplot(aes(x=max_ct, y=mapped, colour=organism)) +
  geom_point()+
  scale_colour_manual(values = c( "#F94144", "#577590")) +
  #geom_smooth(method=lm) +
  theme_minimal() +
  labs(y='Reads Mapped', x='Max. Ct Value')

#what kind of distribution are the data coming from...?
total_tbl %>% ggplot(aes(x=frac_mapped, fill=organism)) +
  geom_histogram(position = "dodge")+
  theme_minimal()
total_tbl %>% ggplot(aes(x=max_ct, fill=organism)) +
  geom_histogram(binwidth = 1, position = "dodge")+
  theme_minimal()
total_tbl %>% ggplot(aes(x=cap_lib_conc, fill=organism)) +
  geom_histogram(position = "dodge")+
  theme_minimal()

#fit a model to frac_mapped  (binomial glm - if we consider frac mapped as ratio of successes to failures)
m0 <- glm(data=total_tbl, frac_mapped ~ max_ct, family=binomial(link="logit"))
m1 <- glm(data=total_tbl, frac_mapped ~ organism*max_ct, family=binomial(link="logit"))
m2 <- glm(data=total_tbl, frac_mapped ~ organism+max_ct, family=binomial(link="logit"))
m3 <- glm(data=total_tbl, frac_mapped ~ organism*I(max_ct^2), family=binomial(link="logit"))
m4 <- glm(data=total_tbl, frac_mapped ~ organism*max_ct*cap_lib_conc, family=binomial(link="logit"))
m5 <- glm(data=total_tbl, frac_mapped ~ organism*max_ct+cap_lib_conc, family=binomial(link="logit"))
summary(m0)
summary(m1)
summary(m2)
summary(m3)
summary(m4)
#which is lowest?
min(m1$aic, m0$aic, m2$aic, m3$aic)
#m1 fits best (organism*max_ct)

#plot the data and model predictions
interactions::interact_plot(m1, pred = max_ct, modx = 'organism', interval = TRUE, plot.points = TRUE, line.thickness = 0.5)
interactions::interact_plot(m1, pred = cap_lib_conc, modx = 'organism', interval = TRUE, plot.points = TRUE, line.thickness = 0.5)
#looks like the model is not fitting very well to the anthrax data (underestimating the lower ct and overestimating at the higher ct)

#let's try fitting to the raw mapped data
#the frac mapped shouldn't matter as we are modelling organism as an interaction again
#so it will consider mapped for each organism and whether the organisms are significantly different to one another

m0 <- MASS::glm.nb(data=total_tbl, mapped ~ max_ct)
m1 <- MASS::glm.nb(data=total_tbl, mapped ~ organism*max_ct)
m2 <- MASS::glm.nb(data=total_tbl, mapped ~ organism+max_ct)
min(m1$aic, m0$aic, m2$aic)
summary(m0)
summary(m1)
summary(m2)

predict(m1)
#m1 fites best, let's plot
modelplot <- interactions::interact_plot(m1, pred = max_ct, modx = 'organism', interval = TRUE, plot.points = TRUE, line.thickness = 0.5)
modelplot
#intuitively this looks like the model fits the data better than the frac_mapped 
#let's stick with this guy
#let's plot with the coverage plot

#now let's get the summary data
total_tbl <- total_tbl %>% mutate(frac_duplicates = duplicates / mapped)
total_tbl %>% if_else(organism == 'mycoplasma', genome_size = 1, genome_size = 2)

total_tbl <-  total_tbl %>% mutate(genome_size = case_when(
  organism == 'mycoplasma' ~ 1029022,
  organism == 'anthrax'~ 5227419
))
total_tbl$mean_doc = (total_tbl$mapped*75)/total_tbl$genome_size

tbldata <- total_tbl %>% group_by(organism) %>% summarise(median_reads = median(total.x), min_reads = min(total.x), max_reads = max(total.x), min_max_ct = min(max_ct), max_max_ct = max(max_ct), median_max_ct = median(max_ct), min_conc = min(cap_lib_conc), max_conc = max(cap_lib_conc),median_conc = median(cap_lib_conc), min_dups = min(frac_duplicates), max_dups = max(frac_duplicates) ,median_duplicate_proportion = median(frac_duplicates), min_mapped = min(frac_mapped), max_mapped = max(frac_mapped), median_mapped_proportion = median(frac_mapped), median_mean_doc = median(mean_doc), min_doc = min(mean_doc), max_doc = max(mean_doc))


tbl <- flextable::flextable(tbldata)
tbl <- flextable::set_header_labels(tbl, median_reads = "Median Reads", median_max_ct = "Median Ct", median_conc = "Median Conc (ng/ul)", median_duplicate_proportion = 'Median Proportion Duplicates', median_mapped_proportion = "Median Proportion Mapped", median_mean_doc = "Median Avg. Depth of Coverage", organism = "Organism")

toprow <- plot_grid(covplot, modelplot, base_asp = 1.618, labels = c('A', 'B'))

total_reads <- total_tbl %>% ggplot(aes(x=organism, y=total.x)) +
  geom_boxplot(width = 0.5)+
  geom_jitter() +
  theme_minimal() +
  labs(x="Organism", y="Total Reads")
frac_mapped <- total_tbl %>% ggplot(aes(x=organism, y=frac_mapped)) +
  geom_boxplot(width = 0.5)+
  geom_jitter() +
  theme_minimal() +
  labs(x="Organism", y="Proportion Mapped Reads")
ct <- total_tbl %>% ggplot(aes(x=organism, y=max_ct)) +
  geom_boxplot(width = 0.5)+
  geom_jitter() +
  theme_minimal() +
  labs(x="Organism", y="Ct")
duplicate_reads <- total_tbl %>% ggplot(aes(x=organism, y=frac_duplicates)) +
  geom_boxplot(width = 0.5)+
  geom_jitter() +
  theme_minimal() +
  labs(x="Organism", y="Proportion of Reads Duplicate")
total_tbl$`%_bases_above_15` <- as.numeric(total_tbl$`%_bases_above_15`)
abovefifteen <- total_tbl %>% ggplot(aes(x=organism, y=`%_bases_above_15`)) +
  geom_boxplot(width = 0.5)+
  geom_jitter() +
  theme_minimal() +
  labs(x="Organism", y="Percent Bases Above 15X")
conc <- total_tbl %>% ggplot(aes(x=organism, y=cap_lib_conc)) +
  geom_boxplot(width = 0.5)+
  geom_jitter() +
  theme_minimal() +
  labs(x="Organism", y="Concentration (ng/ul)")
doc <- total_tbl %>% ggplot(aes(x=organism, y=mean_doc)) +
  geom_boxplot(width = 0.5)+
  geom_jitter() +
  theme_minimal() +
  labs(x="Organism", y="Mean DOC")

doc <- total_tbl %>% ggplot(aes(x=`%_bases_above_15`, y=mean_doc, color=organism)) +
  geom_point()
doc
sequencingstats <-   plot_grid(total_reads, frac_mapped, abovefifteen, duplicate_reads, labels = c('A', 'B', 'C', 'D'))
labstats <- plot_grid(ct,conc, labels = c("A", "B"))
labstats

total_tbl %>% filter(organism == 'mycoplasma') %>% count()

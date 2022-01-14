###################
#Bactocap Mycoplasma Anthrax
#In support of CITATION, DATE
#Tristan Dennis, August 2021

#load and.or install packages we need
pkg = c("tidyverse", "data.table", "sjPlot", "cowplot", "RColorBrewer", "cowplot", "DHARMa", "lme4", "jtools", "huxtable")
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#


#######
#define some functions

#func to read GATK DepthOfCoverage Summary csvs and take filename as a column
read_summarycov <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = basename(flnm))
}

##func to read GATK DepthOfCoverage Interval (quite big, hence fread) csvs and take filename as a column
read_intervalcov <- function(flnm) {
  data.table::fread(flnm) %>% 
    dplyr::select(9) %>% 
    filter(. > 80) %>% #remove
    count() %>% 
    mutate(filename = basename(flnm))
}


#function for getting rid of rubbish from fields (path, suffix, etc)
remove_rubbish <- function(table, suffix) {
  table$filename <- str_replace(table$filename, suffix, " ") #remove suffix from filename
  table$organism =  str_replace(table$organism, "~/Projects/bactocap/datasets/", "") #remove path from organism
  table$organism =  str_replace(table$organism, "/results/", "") #remove rest of path
 # table$sample_id = gsub('(.*)_\\w+', '\\1', table$sample_id)
  return(table)
}


###############
#Read sample data and generate sample metadata dataframes 
#this uses ST1 and ST2

#get our dirs
anthrax_dir <- c("~/Projects/bactocap/datasets/anthrax/results/")
myco_dir <- c("~/Projects/bactocap/datasets/mycoplasma/results/")
dirs <- c(anthrax_dir, myco_dir)

#read in anthrax and mycoplasma metadata (sample data, etc) - ST1
anth_metadata = read.csv('~/Projects/bactocap/metadata/anthrax_sample_data.csv') %>% mutate(organism = 'anthrax')
#formyco we don't have amp cycles yet so add an na column for that field - ST2
myco_metadata = read.csv('~/Projects/bactocap/metadata/myco_sample_data.csv') %>% dplyr::select(organism, sample_id, cap_lib_conc, init_lib_conc, max_ct) %>% add_column(bc_amp_cycles = NA)
#read myco mapping data
myco_mapping = read.csv("~/Projects/bactocap/metadata/myco_mapping.csv")
#readd anth_mapping data
anth_mapping = read.csv("~/Projects/bactocap/metadata/anth_mapping.csv")

#collect mapping data
anth_metadata = anth_metadata %>% left_join(anth_mapping)
myco_metadata = myco_metadata %>% left_join(myco_mapping)
#myco_mapping$


###############
#Collect coverage information from GATK DoC output files, anc join to mapping/flagstat information
sum_tbl <- NULL 
for (dir in dirs){
  tbl <- list.files(path = dir,pattern = ".sample_summary", full.names = T) %>% 
    map_df(~read_summarycov(.)) %>% filter(sample_id != "Total")
  tbl$organism = paste0(dir)
  sum_tbl <- rbind(sum_tbl, tbl)
}

interval_tbl <- NULL 
for (dir in dirs){
  tbl <- list.files(path = dir,pattern = ".sample_interval_statistics", full.names = T) %>% 
    map_df(~read_summarycov(.)) 
  tbl$organism = paste0(dir)
  interval_tbl <- rbind(interval_tbl, tbl)
}

##remove rubbish
sum_tbl = remove_rubbish(sum_tbl, ".sample_summary")
interval_tbl = remove_rubbish(sum_tbl, ".sample_interval_statistics")
#sum_tbl$sample_id <- sub("_[^_]+$", "", sum_tbl$sample_id)

anth_metadata = left_join(anth_metadata, sum_tbl %>% filter(organism == 'anthrax'), by = c('sample_id' = 'sample_id'))
myco_metadata = left_join(myco_metadata, sum_tbl %>% filter(organism == 'mycoplasma'), by = c('sample_id' = 'sample_id'))



anth_metadata = anth_metadata %>% mutate(amisuccessful= case_when(as.numeric(`%_bases_above_15`) > 80 ~ 'yes',
                                                                  as.numeric(`%_bases_above_15`) < 80 ~ 'no'))

myco_metadata = myco_metadata %>% mutate(amisuccessful= case_when(as.numeric(`%_bases_above_15`) > 80 ~ 'yes',
                                                                  as.numeric(`%_bases_above_15`) < 80 ~ 'no'))

anth_metadata = anth_metadata %>% mutate(frac_mapped = mapped/total.x)
anth_metadata = anth_metadata %>% mutate(frac_duplicates = duplicates/total.x)

myco_metadata = myco_metadata %>% mutate(frac_mapped = mapped/total.x)
myco_metadata = myco_metadata %>% mutate(frac_duplicates = duplicates/total.x)


write_csv(anth_metadata, "~/Projects/bactocap/metadata/anthrax_sample_sequencing_coverage_data_st1.csv")
write_csv(myco_metadata, "~/Projects/bactocap/metadata/mycoplasma_sample_sequencing_coverage_data_st2.csv")

#x = anth_metadata %>% select(-total.y, -secondary, -supplementary, -paired, -read1, -read2, -properly_paired, -with_itself_and_mate_mapped, -mate_on_diff_chr, -mate_on_diff_chrover5, -singletons, total.y, organism.y)
#y = myco_metadata %>% select(-total.y, -secondary, -supplementary, -paired, -read1, -read2, -properly_paired, -with_itself_and_mate_mapped, -mate_on_diff_chr, -mate_on_diff_chrover5, -singletons,total.y)



x = anth_metadata %>% select(sample_id, organism.x, total.x, duplicates,mapped, mean, max_ct, bc_amp_cycles, cap_lib_conc, amisuccessful, `%_bases_above_15`)
y = myco_metadata %>% select(sample_id, organism.x, total.x, duplicates,mapped, mean, max_ct, bc_amp_cycles, cap_lib_conc, amisuccessful, `%_bases_above_15`)
total_tbl = rbind(x,y)
total_tbl$cap_lib_conc = as.numeric(total_tbl$cap_lib_conc)
total_tbl$`%_bases_above_15` = as.numeric(total_tbl$`%_bases_above_15`)

total_tbl %>% filter(organism.x == 'anthrax' & sample_id != 'not-sequenced') %>% group_by(bc_amp_cycles, amisuccessful) %>% summarise(countf = n())
total_tbl %>% filter(max_ct > 30 & sample_id != 'not-sequenced') %>% group_by(organism.x, amisuccessful) %>% summarise(n = n())


#join mapping data to coverage data to make ST3 final
#total_tbl = left_join(total_tbl, covtable, by = c('sample_id' = 'sample_id'))

#create fracmapped column
total_tbl$frac_mapped = total_tbl$mapped/total_tbl$total.x

#define successful and unsuccessful samples as frac genome bases over 15X > 80%
total_tbl = total_tbl %>% mutate(amisuccessful= case_when(`%_bases_above_15` > 0.8 ~ 'yes',
                                          `%_bases_above_15` < 0.8 ~ 'no'))
#wo
write_csv(total_tbl, '~/Projects/bactocap/metadata/all_metadata.csv')

#############
#PLOTS

#ct / cap eff
points = total_tbl %>% dplyr::select(max_ct, organism.x, frac_mapped) %>% 
  ggplot(aes(x=max_ct, y=frac_mapped, colour=organism.x))+
  geom_point()+
  scale_color_manual(values=c("deepskyblue1", "darkorange"), labels = c('B. anthracis','M. amphoriforme')) +
  theme_minimal()+
  theme(legend.text = element_text(face = "italic"))+
  labs(x='Ct', y='Capture Efficiency (Proportion of Mapped Reads)', color = 'Organism')

#plot mean doc by organism
meandoc <- total_tbl %>% dplyr::select(organism.x, mean) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y=mean, fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c('B. anthracis','M. amphoriforme'))+
  theme(axis.text.x = element_text(face = "italic"))+
  ylab("Mean Depth-of-Coverage") +
  xlab("Organism")
meandoc

#baited bases above 15
genomeabovefifteen<- total_tbl %>% dplyr::select(organism.x, `%_bases_above_15`) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y= as.numeric(`%_bases_above_15`), fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha=0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c('B. anthracis','M. amphoriforme'))+
  theme(axis.text.x = element_text(face = "italic"))+
  ylab("Proportion Of Genome Above 15X") +
  xlab("Organism")
genomeabovefifteen

#duplicates
duplicates <- total_tbl %>% dplyr::select(organism.x, duplicates, total.x) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y=duplicates/total.x, fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha=0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  ylim(0,0.3) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c('B. anthracis','M. amphoriforme'))+
  theme(axis.text.x = element_text(face = "italic"))+
  ylab("Proportion Duplicates") +
  xlab("Organism")
duplicates

#cap eff
fracmapped = total_tbl %>% dplyr::select(organism.x, frac_mapped, total.x) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y=frac_mapped, fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha=0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c('B. anthracis','M. amphoriforme'))+
  theme(axis.text.x = element_text(face = "italic"))+
  ylab("Proportion Of Total Reads Mapped") +
  xlab("Organism")
fracmapped


#############
#model inference

#fit a model to frac_mapped  (Jess says this is basically logistic regression)
m0 <- glm(data=total_tbl, frac_mapped ~ max_ct, family=binomial(link="logit"))
m1 <- glm(data=total_tbl, frac_mapped ~ organism.x*max_ct, family=binomial(link="logit"))
m2 <- glm(data=total_tbl, frac_mapped ~ organism.x+max_ct, family=binomial(link="logit"))
m3 <- glm(data=total_tbl, frac_mapped ~ organism.x*I(max_ct^2), family=binomial(link="logit"))
m4 <- glm(data=total_tbl, frac_mapped ~ organism.x*max_ct+organism.x*cap_lib_conc, family=binomial(link="logit"))

#create 'cap eff' from a vector of proportion column (successes v failures)
total_tbl$cap_eff = cbind(total_tbl$mapped, (total_tbl$total.x - total_tbl$mapped))

m5 <- glm(data=total_tbl, cap_eff ~ organism.x*max_ct+organism.x*cap_lib_conc, family=binomial(link="logit"))
m6 <- glm(data=total_tbl, cap_eff ~ organism.x*max_ct+organism.x*cap_lib_conc, family=binomial)
m7 <- glm(data=total_tbl, cap_eff ~ organism.x*max_ct+organism.x*cap_lib_conc, family=quasibinomial)
m7a <- glm(data=total_tbl, cap_eff ~ organism.x*max_ct+cap_lib_conc, family=quasibinomial)
m7b <- glm(data=total_tbl, cap_eff ~ organism.x+max_ct+cap_lib_conc, family=quasibinomial)
m7c <- glm(data=total_tbl, cap_eff ~ max_ct+cap_lib_conc, family=quasibinomial)
m7d <- lmer(data=total_tbl, cap_eff ~ max_ct+cap_lib_conc^2, family=quasibinomial)

summary(m7c)
drop1(m7d, test="Chisq")
summary(total_tbl$cap_lib_conc)

effect_plot(m7c, pred = max_ct, interval = TRUE, plot.points = TRUE)
jtools::summ( m7d)

plot_model(m7c, type='pred', terms=c('max_ct'), show.values =T)+theme_minimal()+labs(y='Capture Efficiency', x='Ct', title = 'Model Predictions for Ct ~ Capture Efficiency' )


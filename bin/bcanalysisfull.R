###################
#Bactocap Mycoplasma Anthrax
#In support of CITATION, DATE
#Tristan Dennis, August 2021
#r version R version 4.0.4 
#load and.or install packages we need
pkg = c("tidyverse", "sjPlot", "cowplot", "cowplot", "DHARMa", "lme4", "MuMIn")
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

#setwd
#setwd(getSrcDirectory()[1])
#if running interactively
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


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
anthrax_dir <- c("../datasets/anthrax/results/")
myco_dir <- c("..//datasets/mycoplasma/results/")
dirs <- c(anthrax_dir, myco_dir)

#read in anthrax and mycoplasma metadata (sample data, etc) - ST1
anth_metadata = read.csv('../ancillary/metadata/anthrax_sample_data.csv') %>% mutate(organism = 'anthrax')
#formyco we don't have amp cycles yet so add an na column for that field - ST2
myco_metadata = read.csv('../ancillary/metadata/myco_sample_data.csv') %>% dplyr::select(organism, sample_id, pooled, cap_lib_conc, init_lib_conc, max_ct) %>% add_column(bc_amp_cycles = NA)
#read myco mapping data
myco_mapping = read.csv("../ancillary/metadata/myco_mapping.csv")
#readd anth_mapping data
anth_mapping = read.csv("../ancillary/metadata/anth_mapping.csv")
#read coverage data/mapping from baited regions
sum_tbl = read.csv('../ancillary/metadata/anth_myco_baited_mapping.csv')

#collect mapping data
anth_metadata = anth_metadata %>% left_join(anth_mapping)
myco_metadata = myco_metadata %>% left_join(myco_mapping)

#join species specific metadata to species specific coverage info
anth_metadata = left_join(anth_metadata, sum_tbl %>% filter(organism == 'anthrax'), by = c('sample_id' = 'sample_id'))
myco_metadata = left_join(myco_metadata, sum_tbl %>% filter(organism == 'mycoplasma'), by = c('sample_id' = 'sample_id'))

#get rid of crap cols/select useful ones
x = anth_metadata %>% select(sample_id, organism.x, total.x, duplicates,mapped, mean, max_ct, bc_amp_cycles, pooled, cap_lib_conc, amisuccessful, percent_bases_above_15, mapped_in_baited_region)
y = myco_metadata %>% select(sample_id, organism.x, total.x, duplicates,mapped, mean, max_ct, bc_amp_cycles, pooled, cap_lib_conc, amisuccessful, percent_bases_above_15, mapped_in_baited_region)
#bind
total_tbl = rbind(x,y)

#create frac dup cols and coerce some vals to numeric
total_tbl$frac_mapped = total_tbl$mapped/total_tbl$total.x
total_tbl$frac_duplicates = total_tbl$duplicates/total_tbl$total.x
total_tbl$cap_lib_conc = as.numeric(total_tbl$cap_lib_conc)
total_tbl$`%_bases_above_15` = as.numeric(total_tbl$percent_bases_above_15)

#join mapping data to coverage data to make ST3 final
#total_tbl = left_join(total_tbl, covtable, by = c('sample_id' = 'sample_id'))

#create fracmapped column
total_tbl$frac_mapped = total_tbl$mapped_in_baited_region/total_tbl$total.x

#define successful and unsuccessful samples as frac genome bases over 15X > 80%
total_tbl = total_tbl %>% mutate(amisuccessful= case_when(percent_bases_above_15 > 0.8 ~ 'yes',
                                                          percent_bases_above_15 < 0.8 ~ 'no'))

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

#baited bases above 15
genomeabovefifteen<- total_tbl %>% dplyr::select(organism.x, percent_bases_above_15) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y= as.numeric(percent_bases_above_15), fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha=0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c('B. anthracis','M. amphoriforme'))+
  theme(axis.text.x = element_text(face = "italic"))+
  ylab("Proportion Of Genome Above 15X") +
  xlab("Organism")

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


#plot figure 1
plots = cowplot::plot_grid(meandoc, genomeabovefifteen, fracmapped, nrow=1, labels = c('A', 'B', 'C'))
ggsave(filename = 'fig1plots.tiff', plot = plots, device = grDevices::tiff, path = '../figures_and_tables/', width = 7, height = 2.8)

#############
#model inference

#update dataframe for modelling: add column with 'yes or no'  factor for pooled:unpooled
total_tbl = total_tbl %>% mutate(pooledyesno = (case_when(pooled==1 ~ 'yes', TRUE ~ 'no')))
#change factor labelling for our species to look nicer in plots
total_tbl$organism.x = factor(total_tbl$organism.x, levels = c("anthrax", "mycoplasma"),
       labels = c("B. anthracis", "M. amphoriforme")
)

#create 'cap eff' from a vector of proportion column (successes v failures)
total_tbl$cap_eff = cbind(total_tbl$mapped, (total_tbl$total.x - total_tbl$mapped_in_baited_region))

#add observation level variable to account for high individual variation w/big absolute values
total_tbl$rowid = seq(1, nrow(total_tbl))
total_tbl$rowid = as.factor(total_tbl$rowid)

#final model
m9 = glmer(data=total_tbl, cap_eff ~ max_ct + cap_lib_conc + as.factor(pooledyesno) + (1|rowid), family = binomial)
#drop1 - which covariates are needed? 
drop1(m9, test='Chisq')
#take a look at the model
summary(m9)
#r2 of our main model
MuMIn::r.squaredGLMM(m9)

#run sub-models with dropped covariates to see how informative each covariate is
m9a = glmer(data=total_tbl, cap_eff ~ as.factor(pooledyesno) + (1|rowid), family = binomial)
m9b = glmer(data=total_tbl, cap_eff ~  cap_lib_conc  + (1|rowid), family = binomial)
m9c = glmer(data=total_tbl, cap_eff ~ max_ct + (1|rowid), family = binomial)
#r2 of our manual drop1
MuMIn::r.squaredGLMM(m9a)
MuMIn::r.squaredGLMM(m9b)
MuMIn::r.squaredGLMM(m9c)

modelplot = plot_model(m9, type='pred', terms=c('max_ct'), show.values =T)+theme_minimal()+labs(y='Capture Efficiency', x='Ct', title = 'Model Predictions for Ct ~ Capture Efficiency' )
fig2plot = cowplot::plot_grid(points, modelplot, labels=c('A', 'B'), rel_widths = c(1.3, 1))
ggsave(filename = 'fig2plots.tiff', plot = fig2plot, device = grDevices::tiff, path = '../figures_and_tables/', width = 7, height = 2.8)


####summary stats table
total.sum <- total_tbl %>% filter(sample_id != 'not-sequenced' & sample_id != 'AN16-149-1-T_S26') %>% 
  rename(Total.Reads = total.x,Mapped.reads.in.baited.region = mapped_in_baited_region,Cap.Eff = frac_mapped,Mean.DOC = mean) %>% 
  group_by(organism.x) %>% 
  select(Total.Reads, Mapped.reads.in.baited.region, Cap.Eff, Mean.DOC, `%_bases_above_15`) %>% # select variables to summarise
  summarise_each(funs(min = min, 
                      q25 = quantile(., 0.25), 
                      median = median, 
                      q75 = quantile(., 0.75), 
                      max = max,
                      mean = mean, 
                      sd = sd))

write.csv(total.sum, file = '../figures_and_tables/sumstats.csv')

     
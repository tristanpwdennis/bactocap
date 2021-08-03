###################
#Bactocap Mycoplasma Anthrax
#In support of CITATION, DATE
#Tristan Dennis, August 2021

#load and.or install packages we need
pkg = c("tidyverse", "data.table", "sjPlot", "cowplot", "RColorBrewer", "cowplot", "DHARMa", "interactions", "jtools")
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
myco_metadata = read.csv('~/Projects/bactocap/metadata/myco-full-metadata.csv') %>% dplyr::select(organism, sample_id, cap_lib_conc, init_lib_conc, max_ct) %>% add_column(bc_amp_cycles = NA)
#read myco mapping data
myco_mapping = read.csv("~/Projects/bactocap/datasets/mycoplasma/results/myco_mapping.csv")
#readd anth_mapping data
anth_mapping = read.csv("~/Projects/bactocap/datasets/anthrax/results/anth_mapping.csv")

#collect mapping data
anth_metadata = anth_metadata %>% left_join(anth_mapping)
myco_metadata = myco_mapping %>% left_join(myco_mapping)



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

anth_metadata %>% select(-total.y, -secondary, -supplementary, -paired, -read1, -read2, -properly_paired, -with_itself_and_mate_mapped, -mate_on_diff_chr, -mate_on_diff_chrover5, -singletons, `X`, total.y, organism.y)
myco_metadata %>% select(-total.y, -secondary, -supplementary, -paired, -read1, -read2, -properly_paired, -with_itself_and_mate_mapped, -mate_on_diff_chr, -mate_on_diff_chrover5, -singletons, X, total.y)


#join mapping data to coverage data to make ST3 final
total_tbl = left_join(mappingdata, covtable, by = c('sample_id' = 'sample_id'))

#create fracmapped column
total_tbl$frac_mapped = total_tbl$mapped/total_tbl$total.x

#check everything is ok (e.g. we have the required number of samples for anthrac (93) and myco (56))
total_tbl %>%  group_by(organism) %>% count()

#join to sample metadata
s = left_join(metadata, total_tbl, by =c('sample_id' = 'sample_id'))

#define successful and unsuccessful samples as frac genome bases over 15X > 80%
s = s %>% mutate(amisuccessful= case_when(frac_genome_bases_over_15 > 0.8 ~ 'yes',
                                          frac_genome_bases_over_15 < 0.8 ~ 'no'))

s %>% group_by(organism.x, amisuccessful) %>% count()


write_csv(s, '~/Projects/bactocap/metadata/all_metadata.csv')

#############
#PLOTS


#plot mean doc by organism
meandoc <- s %>% dplyr::select(organism.x, mean) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y=mean, fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  theme(legend.position = "none") +
  ylab("Mean Depth-of-Coverage") +
  xlab("Organism")
meandoc

##plot genome bases above 15X
#baitsabovefifteen<- s %>% dplyr::select(organism.x, frac_genome_bases_over_15) %>% drop_na() %>% 
#  ggplot(aes(x = organism.x, y=as.numeric(frac_genome_bases_over_15), fill=organism.x)) +
#  geom_boxplot(width = 0.5, alpha = 0.7) +
#  geom_jitter(width = 0.1, alpha=0.5) +
#  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
#  theme_minimal() +
#  theme(legend.position = "none") +
#  ylab("Proportion Baits Covered Above 15X") +
#  xlab("Organism")
#baitsabovefifteen

#genome bases above 15
genomeabovefifteen<- s %>% dplyr::select(organism.x, frac_genome_bases_over_15) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y= as.numeric(frac_genome_bases_over_15), fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha=0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  theme(legend.position = "none") +
  ylab("Proportion Baited Genome Above 15X") +
  xlab("Organism")
genomeabovefifteen


#duplicates
duplicates <- s %>% dplyr::select(organism.x, duplicates, total.x) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y=duplicates/total.x, fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha=0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  ylim(0,0.3) +
  theme(legend.position = "none") +
  ylab("Proportion Duplicates") +
  xlab("Organism")
duplicates


cowplot::plot_grid(meandoc, genomeabovefifteen, duplicates,  NULL, labels = c("A", "B", "C"), ncol = 2)


##collect summary stats - FIX
#statstab = s %>%  
#  mutate(frac_duplicates = duplicates/total.x) %>% 
#  select(organism.x,frac_duplicates, frac_genome_bases_over_15, mean) %>% drop_na()
#
#statstab %>% group_by(organism.x) %>% summarise(
#  median_dups = median(as.numeric(frac_genome_bases_over_15)),
#  median
#)  
  
  
#duplicates
d0 <- MASS::glm.nb(data=s, duplicates ~ max_ct)
d1 <- MASS::glm.nb(data=s, duplicates ~ organism.x*max_ct)
summary(d1)
summary(d0)

#fit a model to frac_mapped  (Jess says this is basically logistic regression)
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
summary(m5)

#m1 fits best (organism*max_ct)

#plot the data and model predictions
fracmapped = interactions::interact_plot(m1, pred = max_ct, modx = 'organism', interval = TRUE, plot.points = TRUE, line.thickness = 0.5)
#looks like the model is not fitting very well to the anthrax data (underestimating the lower ct and overestimating at the higher ct)
fracmapped
#let's try fitting to the raw mapped data
#the frac mapped shouldn't matter as we are modelling organism as an interaction again
#so it will consider mapped for each organism and whether the organisms are significantly different to one another

m0 <- MASS::glm.nb(data=s, mapped ~ max_ct)
m1 <- MASS::glm.nb(data=s, mapped ~ organism.x*max_ct)
m2 <- MASS::glm.nb(data=s, mapped ~ organism.x+max_ct)
m3 <- MASS::glm.nb(data=s, mapped ~ organism.x*cap_lib_conc)
m4 <- MASS::glm.nb(data=s, mapped ~ organism.x*max_ct+cap_lib_conc)
min(m1$aic, m0$aic, m2$aic)
summary(m0)
summary(m1)
summary(m2)
summary(m3)
summary(m4)


ggplot(s, aes(x=as.factor(bc_amp_cycles), y=duplicates))+
  geom_jitter()
ggplot(s, aes(x=as.factor(bc_amp_cycles), y=total.x))+
  geom_jitter()

ggplot(atable, aes(x=duplicates, y=mapped))+
  geom_point()

atable %>% ggplot(x=bc_amp_cycles, y=mapped)+
  geom_point()


atable %>% mutate(frac_duplicates = duplicates/total.x) %>% ggplot(aes(x=as.factor(bc_amp_cycles), y=frac_duplicates)) + geom_jitter()
atab

interactions::interact_plot(y1, pred = max_ct, modx = 'lp_amp_cycles', interval = TRUE, plot.points = TRUE, line.thickness = 0.5)

#m1 fites best, let's plot
modelplot <- interactions::interact_plot(m1, pred = max_ct, modx = 'organism', interval = TRUE, plot.points = TRUE, line.thickness = 0.5)
modelplot + labs(x="Max Ct", y ="Mapped Reads", color="Organism")
totalmapped = modelplot
modelplot$data %>% 
  ggplot(aes(x=max_ct, y=mapped))+
  geom_point()
modelplot$layers
totalmapped

antbl = total_tbl %>% filter(organism == 'anthrax')
anm0 <- MASS::glm.nb(data=total_tbl, mapped ~ max_ct)

jtools::effect_plot(anm0, pred = max_ct, interval = TRUE, plot.points = TRUE, colors = "#f58d42", point.color = '#f58d42')
jtools::effect_plot(y1)

plot1
plot1 + theme_sjplot()

#for assessing fit of nb
#install.packages("DHARMa")
simulationOutput <- DHARMa::simulateResiduals(fittedModel = m1, plot = T)
hist(residuals(simulationOutput))

s_tbl<-total_tbl
s_tbl<-rename(s_tbl, fifteen = `%_bases_above_15`)
s_tbl$fifteen <- s_tbl$fifteen/100
d0  <- glm(data=s_tbl,fifteen ~ max_ct, family=binomial(link="logit"))
d1  <- glm(data=s_tbl,fifteen ~ organism*max_ct, family=binomial(link="logit"))
d2  <- glm(data=s_tbl,fifteen ~ organism*I(max_ct^2), family=binomial(link="logit"))

summary(d2)
fifteen<-interactions::interact_plot(d2, pred = max_ct, modx = 'organism',interval = TRUE, plot.points = TRUE, line.thickness = 0.5)

#how many reads to cover genome to 15X
newdf <- data.frame(
  organism = c("anthrax", "mycoplasma"),
  max_ct = c(39, 34)
)

predict(m1, newdf, type = "response", interval = 'prediction')


legend <- get_legend(
  # create some space to the left of the legend
  totalmapped + theme(legend.box.margin = margin(0, 0, 0, 12))
)
legend

totalmapped = totalmapped + xlab("Max Ct") + ylab("Total Mapped Reads") 
fracmapped =  fracmapped  + xlab("Max Ct")  + ylab("Proportion Mapped Reads") +theme_minimal()+ theme(legend.position = "none")
fifteen = fifteen + xlab("Max Ct") + ylab("Proportion Baited Region > 15X") +theme_minimal()+ theme(legend.position = "none")
totalmapped


total_anth = filter(total_tbl, organism == 'anthrax')
total_anth %>% ggplot(aes(x=bc_amp_cycles, y=duplicates))+geom_boxplot()




#plot pool/frac mapped reads
tlfmdata = read.csv('~/Projects/bactocap/metadata/tlf_mdata.csv')
tlfmdata$sample_id = paste0(tlfmdata$Carcass.ID, "-", toupper(substring(tlfmdata$sample.type, 1, 1)))
anthmdata = read.csv("~/Projects/bactocap/metadata/anthrax-metadata.csv")
full_anth_mdata = left_join(tlfmdata, anthmdata)

anmapping = read.csv("~/Projects/bactocap/datasets/anthrax/results/anth_mapping.csv")
final_full_mdata = left_join(full_anth_mdata, anmapping, by = c('polyom_id' = 'sample_id'))
write.csv(final_full_mdata, '~/Projects/bactocap/metadata/full_anth_metadata_ud.csv')
mdata = read.csv('~/Projects/bactocap/metadata/full_anth_metadata.csv')
mdata %>% mutate(frac_mapped = mapped/total) %>% ggplot(aes(x=Pooled, y=frac_mapped, fill=Pooled))+
  geom_jitter()+
  geom_boxplot(width = 0.3) +
  theme_minimal()+
  scale_fill_manual(values = c("#577590", "#F94144")) +
  labs(x='Pooled Yes/No', y='Proportion Mapped Reads')

ggplot(mdata, aes(x=mdata$Pooled, y=frac_mapped))


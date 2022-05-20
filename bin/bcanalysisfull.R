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

#read metadata
total_tbl = read.csv('~/Projects/bactocap/ancillary/metadata/all_metadata.csv')


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
ggsave(filename = 'fig1plots.tiff', plot = plots, device = grDevices::tiff, path = '../figures_and_tables/', width = 10, height = 3.4)

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

m9a = glmer(data=total_tbl, cap_eff ~ as.factor(pooledyesno) + (1|rowid), family = binomial)
m9b = glmer(data=total_tbl, cap_eff ~  cap_lib_conc  + (1|rowid), family = binomial)
m9c = glmer(data=total_tbl, cap_eff ~ max_ct + (1|rowid), family = binomial)
#r2 of our manual drop1
MuMIn::r.squaredGLMM(m9a)
MuMIn::r.squaredGLMM(m9b)
MuMIn::r.squaredGLMM(m9c)

modelplot = plot_model(m9, type='pred', terms=c('max_ct'), show.values =T)+theme_minimal()+labs(y='Capture Efficiency', x='Ct', title = 'Model Predictions for Ct ~ Capture Efficiency' )
fig2plot = cowplot::plot_grid(points, modelplot, labels=c('A', 'B'), rel_widths = c(1.3, 1))
ggsave(filename = 'fig2plots.tiff', plot = fig2plot, device = grDevices::tiff, path = '../figures_and_tables/', width = 10, height = 4.2)



####summary stats table
total.sum <- total_tbl %>% filter(sample_id != 'not-sequenced' & sample_id != 'AN16-149-1-T_S26') %>% 
  rename(Total.Reads = total.x,Mapped.reads.in.baited.region = mapped_in_baited_region,Cap.Eff = frac_mapped,Mean.DOC = mean) %>% 
  group_by(organism.x) %>% 
  select(Total.Reads, Mapped.reads.in.baited.region, Cap.Eff, Mean.DOC, `percent_bases_above_15`) %>% # select variables to summarise
  summarise_each(funs(min = min, 
                      q25 = quantile(., 0.25), 
                      median = median, 
                      q75 = quantile(., 0.75), 
                      max = max,
                      mean = mean, 
                      sd = sd))

write.csv(total.sum, file = '../figures_and_tables/sumstats.csv')

####plot ct vs genome coverage
sf_ct_menadoc = total_tbl %>% dplyr::select(organism.x, mean,max_ct ) %>% drop_na() %>% 
  ggplot(aes(x=max_ct, y=mean, colour=organism.x))+
  geom_point()+
  scale_colour_manual(values=c("deepskyblue1", "darkorange"), name = "Organism")+
  theme_minimal() +
  ylab("Mean Depth-of-Coverage") +
  xlab("Ct")

ggsave(filename = 'fig_s3.tiff', plot = sf_ct_menadoc, device = grDevices::tiff, path = '../figures_and_tables/', width = 5, height = 3)

#ploit effect of pooling by organism on cap eff
pooledplot = ggplot(total_tbl, aes(x=pooledyesno, y=frac_mapped))+
  facet_wrap(~organism.x)+
  geom_boxplot()+
  geom_jitter()+
  theme_minimal()+
  labs(x='Was the sample pooled?', y='Capture efficiency')
ggsave(filename = 'fig_s2.tiff', plot = pooledplot, device = grDevices::tiff, path = '../figures_and_tables/')

#plot model predictions for cap eff vs cap lib conc
capeff_clc_plot = plot_model(m9, type='pred', terms=c('cap_lib_conc'), show.values =T)+theme_minimal()+labs(y='Capture Efficiency', x='Captured Library Concentration', title = 'Captured Library Concentration ~ Capture Efficiency' )

ggsave(filename = 'fig_s1.tiff', plot = capeff_clc_plot, device = grDevices::tiff, path = '../figures_and_tables/')


#parse kraken reports and turn into supp tabs
mycokrakenreps = list.files('../datasets/mycoplasma/results/', pattern = 'flatkrakrep')
ankrakenreps = list.files('../datasets/anthrax/results/', pattern = 'flatkrakrep')
mycokrakdata = lapply(seq_along(mycokrakenreps), function(i){fread(paste0('../datasets/mycoplasma/results/',mycokrakenreps[[i]])) %>% mutate(fn=mycokrakenreps[[i]])}) #fread all, adding filename as col
anthkrakdata = lapply(seq_along(ankrakenreps), function(i){fread(paste0('../datasets/anthrax/results/',ankrakenreps[[i]])) %>% mutate(fn=ankrakenreps[[i]])}) #fread all, adding filename as col
mycokrakdata = do.call(rbind, mycokrakdata) 
anthkrakdata = do.call(rbind, anthkrakdata) 
anthkrakdata$organism = 'B. anthracis'
mycokrakdata$organism = 'M. amphoriforme'
allkrakdata = rbind(mycokrakdata, anthkrakdata)
allkrakdata = allkrakdata[V1 > 1]
allkrakdata = allkrakdata[V4 == 'G' | V4 == 'U']
allkrakdata$fn = gsub('.flatkrakrep', '', allkrakdata$fn)
antab = allkrakdata %>% dplyr::select(organism, fn, V1, V6) %>% filter(organism == 'B. anthracis')
mycotab = allkrakdata %>% dplyr::select(organism, fn, V1, V6) %>% filter(organism == 'M. amphoriforme')
write.csv(mycotab, file = '../figures_and_tables/st6.csv')


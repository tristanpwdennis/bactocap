###################
#Bactocap Mycoplasma Anthrax
#In support of CITATION, DATE
#Tristan Dennis, August 2021

#load and.or install packages we need
pkg = c("tidyverse", "data.table", "cowplot", "RColorBrewer", "cowplot", "DHARMa", "interactions", "jtools")
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

anth = read.csv('~/Projects/bactocap/metadata/anthrax_sample_sequencing_coverage_data_st1.csv')
myco = read.csv('~/Projects/bactocap/metadata/mycoplasma_sample_sequencing_coverage_data_st2.csv')

anth = anth %>% select(sample_id, mean, max_ct, total.x, init_lib_conc, cap_lib_conc, X._bases_above_15, frac_duplicates, frac_mapped, duplicates, mapped, amplification_cycles_bait_capture, amisuccessful) %>% mutate(organism.x = 'anthrax')
myco = myco %>% select(sample_id, mean, max_ct, total.x, init_lib_conc, cap_lib_conc, X._bases_above_15, frac_duplicates, frac_mapped, duplicates, mapped, amisuccessful) %>% mutate(amplification_cycles_bait_capture = NA) %>% mutate(organism.x = 'mycoplasma')

s = rbind(anth, myco)

xtab = s %>% filter(sample_id != 'AN16-149-1-T_S26') %>% select(frac_duplicates, X._bases_above_15, mean, frac_mapped, mapped, 'total.x', organism.x) %>% group_by(organism.x) %>% drop_na() %>% 
summarise(edian_total = median('total.x'),
          ax_total = max('total.x'),
          in_total = min('total.x'),
          ean_total = mean('total.x'),
          edian_map = median(mapped),
          ax_map = max(mapped),
          in_map = min(mapped),
          ean_map = mean(mapped),
          edian_fmap = median(frac_mapped),
          ax_fmap = max(frac_mapped),
          in_fmap = min(frac_mapped),
          ean_fmap = mean(frac_mapped),
          edian_meandoc = median(mean),
          ax_meandoc = max(mean),
          in_meandoc = min(mean),
          ean_meandoc = mean(mean),
          qr_meandoc = IQR(mean),
          pper_quartile_meandoc = quantile(mean, 0.75),
          ower_quartile_meandoc = quantile(mean, 0.25),
          median_fabove15 = median(X._bases_above_15),
          max_fabove15 = max(X._bases_above_15),
          min_fabove15 = min(X._bases_above_15),
          mean_fabove15 = mean(X._bases_above_15),
          iqr_fabove15 = IQR(X._bases_above_15),
          upper_quartile_fabove15 = quantile(X._bases_above_15, 0.75),
          lower_quartile_fabove15 = quantile(X._bases_above_15, 0.25),
          median_dup = median(frac_duplicates),
          max_dup = max(frac_duplicates),
          min_dup = min(frac_duplicates),
          mean_dup = mean(frac_duplicates),
          iqr_dup = IQR(frac_duplicates),
          upper_quartile_dup = quantile(frac_duplicates, 0.75),
          lower_quartile_dup = quantile(frac_duplicates, 0.25),
          median_fracmap = median(frac_mapped),
          max_fracmap = max(frac_mapped),
          min_fracmap = min(frac_mapped),
          mean_fracmap = mean(frac_mapped),
          iqr_fracmap = IQR(frac_mapped),
          upper_quartile_fracmap = quantile(frac_mapped, 0.75),
          lower_quartile_fracmap = quantile(frac_mapped, 0.25)
          
          
          )

write_csv(xtab, "~/Projects/bactocap/metadata/supp_summ_stats.csv")
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

#genome bases above 15
genomeabovefifteen<- s %>% dplyr::select(organism.x, percent_bases_above_15) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y= as.numeric(percent_bases_above_15)/100, fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha=0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  theme(legend.position = "none") +
  ylab("Proportion") +
  xlab("Organism")
genomeabovefifteen

#duplicates
duplicates <- s %>% dplyr::select(organism.x, frac_duplicates) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y=frac_duplicates, fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha=0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  ylim(0,1) +
  theme(legend.position = "none") +
  ylab("Proportion") +
  xlab("Organism")
duplicates

fracmap <- s %>% dplyr::select(organism.x, frac_mapped) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y=frac_mapped, fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha=0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  ylim(0,1) +
  theme(legend.position = "none") +
  ylab("Proportion") +
  xlab("Organism")
fracmap

#plot final
cowplot::plot_grid(meandoc, genomeabovefifteen, duplicates, fracmap, NULL, labels = c("A", "B", "C", "D"), ncol = 2)



s %>% ggplot(aes(x=amisuccessful,y=as.numeric(cap_lib_conc),color=organism.x))+
  geom_point()+
  facet_wrap(~organism.x)+
  theme_minimal()


s %>% select(organism.x, max_ct) %>% group_by(organism.x) %>% summarise(min(max_ct), max(max_ct))

###############
#MODELS

#duplicates
d0 <- MASS::glm.nb(data=s, duplicates ~ max_ct)
d1 <- MASS::glm.nb(data=s, duplicates ~ organism.x*max_ct)
summary(d1)
summary(d0)

#fit a model to frac_mapped  (Jess says this is basically logistic regression)
m0 <- glm(data=s, frac_mapped ~ max_ct, family=binomial(link="logit"))
m1 <- glm(data=s, frac_mapped ~ organism.x*max_ct, family=binomial(link="logit"))
m2 <- glm(data=s, frac_mapped ~ organism.x+max_ct, family=binomial(link="logit"))
m3 <- glm(data=s, frac_mapped ~ organism.x*I(max_ct^2), family=binomial(link="logit"))
m4 <- glm(data=s, frac_mapped ~ organism.x*max_ct*cap_lib_conc, family=binomial(link="logit"))
m5 <- glm(data=s, frac_mapped ~ organism.x*max_ct+cap_lib_conc, family=binomial(link="logit"))

summary(m0)
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)

#m1 fits best (organism*max_ct)

#plot the data and model predictions
fracmapped = interactions::interact_plot(m1, pred = max_ct, modx = 'organism.x', interval = TRUE, plot.points = TRUE, line.thickness = 0.5)
#looks like the model is not fitting very well to the anthrax data (underestimating the lower ct and overestimating at the higher ct)
fracmapped + labs(x="Max Ct", y ="Proportion Mapped Reads", color="Organism") + theme(legend.position = "none")#add labels
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

#m1 fites best, let's plot interxns
modelplot <- interactions::interact_plot(m1, pred = max_ct, modx = 'organism.x', interval = TRUE, plot.points = TRUE, line.thickness = 0.5)
modelplot + labs(x="Max Ct", y ="Mapped Reads", color="Organism") + theme(legend.position = "none")#add labels


#for assessing fit of nb
#install.packages("DHARMa")
simulationOutput <- DHARMa::simulateResiduals(fittedModel = m1, plot = T)
hist(residuals(simulationOutput))

#####bc cycles and mapped reads
s %>% ggplot(aes(x=as.factor(amplification_cycles_bait_capture), y=mapped)) + geom_jitter()
y0 = MASS::glm.nb(data=s, mapped ~ as.factor(amplification_cycles_bait_capture))
summary(y0)

s %>% filter(sample_id != 'not_sequenced') %>% 
ggplot(., aes(x=as.numeric(cap_lib_conc), y=as.numeric(init_lib_conc), color=organism.x))+
  geom_point()+
  theme_minimal()
  




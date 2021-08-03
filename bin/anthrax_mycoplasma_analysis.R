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

anth = anth %>% select(sample_id, mean, max_ct, cap_lib_conc, organism.x, X._bases_above_15, frac_duplicates, frac_mapped, duplicates, mapped, amplification_cycles_bait_capture)
myco = myco %>% select(sample_id, mean, max_ct, cap_lib_conc, organism.x, X._bases_above_15, frac_duplicates, frac_mapped, duplicates, mapped) %>% mutate(amplification_cycles_bait_capture = NA)

s = rbind(anth, myco)

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
genomeabovefifteen<- s %>% dplyr::select(organism.x, X._bases_above_15) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y= as.numeric(X._bases_above_15), fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha=0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  theme(legend.position = "none") +
  ylab("Proportion Baited Genome Above 15X") +
  xlab("Organism")
genomeabovefifteen

#duplicates
duplicates <- s %>% dplyr::select(organism.x, frac_duplicates) %>% drop_na() %>% 
  ggplot(aes(x = organism.x, y=frac_duplicates, fill=organism.x)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha=0.5) +
  scale_fill_manual(values=c("deepskyblue1", "darkorange")) +
  theme_minimal() +
  ylim(0,0.3) +
  theme(legend.position = "none") +
  ylab("Proportion Duplicates") +
  xlab("Organism")
duplicates

#plot final
cowplot::plot_grid(meandoc, genomeabovefifteen, duplicates,  NULL, labels = c("A", "B", "C"), ncol = 2)


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

#m1 fites best, let's plot interxns
modelplot <- interactions::interact_plot(m1, pred = max_ct, modx = 'organism.x', interval = TRUE, plot.points = TRUE, line.thickness = 0.5)
modelplot + labs(x="Max Ct", y ="Mapped Reads", color="Organism") #add labels
modelplot

#for assessing fit of nb
#install.packages("DHARMa")
simulationOutput <- DHARMa::simulateResiduals(fittedModel = m1, plot = T)
hist(residuals(simulationOutput))

#####bc cycles and mapped reads
s %>% ggplot(aes(x=as.factor(amplification_cycles_bait_capture), y=mapped)) + geom_jitter()
y0 = MASS::glm.nb(data=s, mapped ~ as.factor(amplification_cycles_bait_capture))
summary(y0)
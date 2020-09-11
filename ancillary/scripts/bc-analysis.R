#####examining coverage over bait sequences for mycoplasma and anthrax
library(viridis)
library(tidyverse)
library(GGally)
library(reticulate)
#get metadata for mycoplasma and anthrax
mycodata <- read.csv("/Users/tristanpwdennis/Projects/bactocap-data/mycoplasma/metadata/myco-seq-sample-data.csv", na.strings = c("", "NA"))
anthdata <- read.csv("/Users/tristanpwdennis/Projects/bactocap-data/anthrax/metadata/anthrax-metadata-1.csv",na.strings = c("", "NA"))




#plot covariance
bcdata %>% 
  select(total_reads, highest_ct, captured_lib_conc, rate_mapped_reads, pcnt_30X_captures,bc_camp_cycles, lp_amp_cycles ) %>% 
  sample_frac(0.5, replace = TRUE) %>% 
  ggpairs(axislabels = "none") +
  theme_bw()

mycodata$frac = (24444 - mycodata$baits_covered_.30)/24444
my$sampleid <-sub("_[^_]+$", "", my$sampleid)
mycodata$frac_mapped_reads = mycodata$mapped/mycodata$total
mycodat <- mycodata %>% select(total, frac_mapped_reads, captured_amamplified_library_conc_ngul, udg_PCR_run1, frac)
colnames(mycodat) <- c('total_reads', 'rate_mapped_reads', 'captured_lib_conc', 'highest_ct', 'pcnt_30X_captures')


mycodat %>% select(total_reads, highest_ct, captured_lib_conc, rate_mapped_reads, pcnt_30X_captures) %>% 
  sample_frac(0.5, replace = TRUE) %>% 
  ggpairs(axislabels = "none") +
  theme_bw()


  
summary(mycodata$frac_mapped_reads)
summary(mycodata$dupfrac)




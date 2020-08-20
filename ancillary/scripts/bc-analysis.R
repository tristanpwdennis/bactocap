#####examining coverage over bait sequences for mycoplasma and anthrax

library(viridis)
library(data.table)
library(tidyverse)
library(GGally)

#get metadata for mycoplasma and anthrax
mycodata <- read.csv("/Users/tristanpwdennis/Projects/bactocap-data/mycoplasma/metadata/myco-seq-sample-data.csv", na.strings = c("", "NA"))
anthdata <- read.csv("/Users/tristanpwdennis/Projects/bactocap-data/anthrax/metadata/anthrax-metadata-1.csv",na.strings = c("", "NA"))

#remove trailing underscore
#mycoseqdata$sample_id <- sub("_[^_]+$", "", mycoseqdata$sample_id)

#load paths for coverage info
anthraxpath <- c("/Users/tristanpwdennis/Projects/bactocap/datasets/anthrax/results/0to400covstats.txt")
mycopath <- c("/Users/tristanpwdennis/Projects/bactocap/datasets/mycoplasma/covstatstotal.txt")

#function to get our files and derive the fraction of bait captures for each coverage level
get_bait_cap_frac <- function(fp, numbaits) {
  cov <- read.csv(fp, sep  = ' ', header = F)
  colnames(cov) <- c('no_baits', 'sampleid', 'coverage')
  cov2 <- mutate(cov, frac = (numbaits - no_baits)/numbaits*100)
  return(cov2)
}

an <- get_bait_cap_frac(anthraxpath, 148729)
my <- get_bait_cap_frac(mycopath, 24444)

#plot fraction
ggplot(my, aes(x = coverage, y = frac, color = sampleid)) +
  geom_line(show.legend = FALSE) +
  theme(legend.position = "none") +
  theme_minimal() 

#ok now we combine into one big dataframe
my$organism <- paste0('mycoplasma')
my$ninds <- nrow(mycodata)
an$organism <- paste0('anthrax')
an$ninds <- nrow(anthdata)

#join anthrax bait frac data with metadata, clean some dud values, change percentages to decimals
bcdata <- left_join(an, anthdata, by = c('sampleid' = 'sample_id')) %>% filter(coverage == 30) %>% select(sampleid, highest_ct, captured_lib_conc, reads, frac_mapped_reads, frac)
bcdata$frac_mapped_reads <- as.numeric(sub("%", "",bcdata$frac_mapped_reads,fixed=TRUE))/100
bcdata <- bcdata %>% na_if('x') %>% na_if('ND')

#plot covariance
bcdata %>% 
  select(frac, highest_ct, captured_lib_conc, reads, frac_mapped_reads) %>% 
  sample_frac(0.5, replace = TRUE) %>% 
  ggpairs(axislabels = "none") +
  theme_bw()

mcdata <- left_join(my, mycodata, by = c('sampleid' = 'ID_on_remote')) %>% filter(coverage == 30) %>% select(sampleid, highest.Ct, captured.lib.conc, reads, Percentage.mapped.reads, frac)
bcdata$Percentage.mapped.reads <- as.numeric(sub("%", "",bcdata$Percentage.mapped.reads,fixed=TRUE))/100
bcdata <- bcdata %>% na_if('x') %>% na_if('ND')

mycodata$frac = mycodata$baits_covered_.30/24444
mycodata$sample_id
my$sampleid <-sub("_[^_]+$", "", my$sampleid)
mycodata$frac_mapped_reads = mycodata$mapped/mycodata$total

mycodata %>% select(total, mapped, captured_amamplified_library_conc_ngul, udg_PCR_run1, frac) %>% 
  sample_frac(0.5, replace = TRUE) %>% 
  ggpairs(axislabels = "none") +
  theme_bw()

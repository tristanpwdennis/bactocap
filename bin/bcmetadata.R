library(tidyverse)

######################################################
#metadata builder for bactocap
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

anthflag <- read.csv('../datasets/anthrax/results/anth_mapping.csv')
mycoflag <- read.csv('../datasets/mycoplasma/results/myco_mapping.csv')
#mlstflag <- read.csv('../datasets/mlst/results/myco_mapping.csv')

anthflag <- anthflag %>% mutate(dataset = 'anthrax')
mycoflag <- mycoflag %>% mutate(dataset = 'mycoplasma')

anthmycoflag <- rbind(anthflag, mycoflag)
anthmycoflag$sample_id<-sub("_[^_]+$", "", anthmycoflag$sample_id)

anmeta <- read.csv('../metadata/anthrax-metadata.csv')
mymeta <- read.csv('../metadata/mycoplasma-metadata.csv')

anmymeta <- rbind(anmeta, mymeta)

t<- anmymeta %>% left_join(anthmycoflag)  %>% filter(!is.na(dataset))
t %>% filter(dataset == 'anthrax') %>% count()



metadatabuild <- function() {
}
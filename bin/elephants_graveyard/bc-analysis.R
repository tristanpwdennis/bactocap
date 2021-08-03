#####examining coverage over bait sequences for mycoplasma and anthrax

library(dplyr)
library(ggplot2)
library(viridis)
library(data.table)

#get metadata for mycoplasma and anthrax
mycodata <- read.csv("/Users/tristanpwdennis/Projects/bactocap-data/mycoplasma/metadata/myco-seq-sample-data.csv")
anthdata <- read.csv("/Users/tristanpwdennis/Projects/bactocap-data/anthrax/metadata/anthrax-metadata-1.csv")

#remove trailing underscore
#mycoseqdata$sample_id <- sub("_[^_]+$", "", mycoseqdata$sample_id)

#load paths for coverage info
anthraxpath <- c("/Users/tristanpwdennis/Projects/bactocap/datasets/anthrax/results/covstats.txt")
mycopath <- c("/Users/tristanpwdennis/Projects/bactocap/datasets/mycoplasma/covstats.txt")
anthraxhi <- c("/Users/tristanpwdennis/Projects/bactocap/datasets/anthrax/results/200to400covstats.txt")
mycohi <- c("/Users/tristanpwdennis/Projects/bactocap/datasets/mycoplasma/filtered_200to400covstats.txt")

#function to get our files and derive the fraction of bait captures for each coverage level
get_bait_cap_frac <- function(fp, numbaits) {
  cov <- read.csv(fp, sep  = ' ', header = F)
  colnames(cov) <- c('no_baits', 'sampleid', 'coverage')
  cov2 <- mutate(cov, frac = (numbaits - no_baits)/numbaits*100)
  return(cov2)
}

an <- bind_rows(get_bait_cap_frac(anthraxpath, 148729), get_bait_cap_frac(anthraxhi, 148729))
my <- bind_rows(get_bait_cap_frac(mycopath, 24444), get_bait_cap_frac(mycohi, 24444))

ggplot(my, aes(x = coverage, y = frac, color = sampleid)) +
  geom_line(show.legend = FALSE) +
  theme(legend.position = "none") +
  theme_minimal() 




#call function and concatenate results for both anthrax and mycoplasma into single df
concatenated_frac_cov_ind <- rbind(get_noinds_cov(an, anthdata, 'anthrax'), get_noinds_cov(my, mycodata, 'mycoplasma'))

#weird result for 200 as I split the calc over two windows - will rerun and plot again
t<-filter(concatenated_frac_cov_ind, capture_rate != 200)

#plot
ggplot(t, aes(x = capture_rate, y = prop_ind, color = organism)) +
  geom_line()+
  xlab("Depth of coverage over 90% of the captured baits") + ylab("Fraction of Samples") +
  theme_minimal() 



my$organism <- paste0('mycoplasma')
my$ninds <- nrow(mycodata)
an$organism <- paste0('anthrax')
an$ninds <- nrow(anthdata)

comb <- rbind(an, my)

#find the cartesian product of the coverage and genome fraction ranges (every % genome covered at every d.o.c), and organism
covfrac <- expand.grid(seq(10,400, 10), seq (0.01, 1, 0.01), as.numeric(c('56', '117')))

testfunc <- function(se, fracse, numin) {
  x <- nrow(filter(comb, coverage == se & frac > fracse & ninds == numin))/numin
}

results <- as.data.frame(mapply(testfunc, covfrac[,1], covfrac[2], covfrac[,3]))
t <- cbind(results, covfrac)























#percent of samples, per anthrax and mycoplasma, covered to each X

get_noinds_cov <- function(df, metadata, orgname) {
  prop <- data.frame()
  for (i in seq(10,400, 10)) {
    x <- nrow(filter(df, coverage == i & frac > 0.9))/nrow(metadata)
    prop = rbind(prop, x)
  }
  prop <- cbind(prop, as.data.frame(seq(10,400, 10)))
  colnames(prop)<-c('prop_ind', 'capture_rate')
  prop$organism <- paste0(orgname)
  return(prop)
}


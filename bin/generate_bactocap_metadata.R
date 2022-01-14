###################
#Bactocap Mycoplasma Anthrax
#In support of CITATION, DATE
#Tristan Dennis, August 2021

#load and.or install packages we need
pkg = c("tidyverse", "data.table")
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
  #table$sample_id = gsub('(.*)_\\w+', '\\1', table$sample_id)
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
myco_metadata = read.csv('~/Projects/bactocap/metadata/myco_sample_data.csv') 
#read myco mapping data
myco_mapping = read.csv("~/Projects/bactocap/metadata/myco_mapping.csv") 
#myco_mapping$sample_id = gsub('(.*)_\\w+', '\\1', myco_mapping$sample_id)
#readd anth_mapping data
anth_mapping = read.csv("~/Projects/bactocap/metadata/anth_mapping.csv")

t = as.data.frame(cbind(myco_mapping$sample_id, gsub('(.*)_\\w+', '\\1', myco_mapping$sample_id)))
colnames(t) = c('sample_id', 'short_id')
myco_metadata = left_join(t, myco_metadata, by = c('short_id' = 'sample_id'))
myco_metadata = myco_metadata %>% select(-short_id)

#collect mapping data
anth_metadata = anth_metadata %>% left_join(anth_mapping)
#anth_metadata$sample_id = gsub('(.*)_\\w+', '\\1', anth_metadata$sample_id)
myco_metadata = myco_metadata %>% left_join(myco_mapping, by = c('sample_id' = 'sample_id'))



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
interval_tbl = remove_rubbish(interval_tbl, ".sample_interval_statistics")




#join to cov tbl
anth_metadata = left_join(anth_metadata, sum_tbl %>% filter(organism == 'anthrax'), by = c('sample_id' = 'sample_id'))
myco_metadata = left_join(myco_metadata, sum_tbl %>% filter(organism == 'mycoplasma'), by = c('sample_id' = 'sample_id'))


###define successful samples (80% genome over 15X)
anth_metadata = anth_metadata %>% mutate(amisuccessful= case_when(as.numeric(`%_bases_above_15`) > 80 ~ 'yes',
                                                                  as.numeric(`%_bases_above_15`) < 80 ~ 'no'))
myco_metadata = myco_metadata %>% mutate(amisuccessful= case_when(as.numeric(`%_bases_above_15`) > 80 ~ 'yes',
                                                                  as.numeric(`%_bases_above_15`) < 80 ~ 'no'))

#calculate frac mapped and frac dups for both sample sets
anth_metadata = anth_metadata %>% mutate(frac_mapped = mapped/total.x)
anth_metadata = anth_metadata %>% mutate(frac_duplicates = duplicates/total.x)

myco_metadata = myco_metadata %>% mutate(frac_mapped = mapped/total.x)
myco_metadata = myco_metadata %>% mutate(frac_duplicates = duplicates/total.x)

#remove gatk/flagstat garbage
t= anth_metadata %>% select(-total.y, -secondary, -supplementary, -polyomics, -paired, -read1, -read2, -properly_paired, -with_itself_and_mate_mapped, -mate_on_diff_chr, -mate_on_diff_chrover5, -singletons, -total.y, -organism.y, -organism.x, -filename)
colnames(t)  = c('sample_id', 'carcass_id', 'sequence_project_id', 'library_prep_batch', 'pool', 'date_collected', 'species', 'sample_type', 'max_ct', 'ct_cap', 'ct_lef', 'ct_PL3', 'sample_conc_ng_ul', 'neb_protocol', 'fragmentation_incubation', 'lib_amp_cycles', 'init_lib_conc', 'concentration_of_lib_performed', 'new_lib_conc_ng_ul', 'amplification_cycles_bait_capture', 'cap_lib_conc', 'cap_lib_reamp', 'total_reads', 'duplicates', 'mapped', 'mean', 'granular_third_quartile', 'granular_median', 'granular_first_quartile', 'percent_bases_above_15', '80_percent_over_15X', 'frac_mapped', 'frac_duplicates')
myco_metadata = myco_metadata %>% select(-total.y, -secondary, -supplementary, -paired, -read1, -read2, -properly_paired, -with_itself_and_mate_mapped, -mate_on_diff_chr, -mate_on_diff_chrover5, -singletons, -total.y, -organism.y, -organism.x, -filename)
colnames(myco_metadata) = c('sample_id', 'Number', 'ari_study', 'study_code', 'date', 'extraction', 'udg_PCR_run1', 'udg_PCR_run2', 'udg_PCR_run3_interp', 'udg_PCR_run1_interp', 'extraction.1', 'udg_PCR_run2_interp', 'udg_PCR_run3', 'qubit_HS_ngul', 'dilution_for_qubit', 'init_lib_conc', 'cap_lib_conc', 'Lib', 'max_ct', 'total_reads', 'duplicates', 'mapped', 'X',  'mean', 'granular_third_quartile', 'granular_median', 'granular_first_quartile', 'percent_bases_above_15', 'amisuccessful', 'frac_mapped', 'frac_duplicates')

#wo
write_csv(t, "~/Projects/bactocap/metadata/anthrax_sample_sequencing_coverage_data_st1.csv")
write_csv(myco_metadata, "~/Projects/bactocap/metadata/mycoplasma_sample_sequencing_coverage_data_st2.csv")








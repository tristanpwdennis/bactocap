#ls pkg
pkg = c("tidyverse", "cowplot", "RColorBrewer", "cowplot", "pheatmap" , "data.table")
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#
#setwd to source dir (w/data files too)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#read our shit
metadata <- fread('mlst_metadata.csv') #read sample data
mapstats <- fread('idxstats.txt') #read long=form mapping data (no reads mapped/locus/sample)
targinfo <- fread('targstats.csv')
#add colnames and remove polyomics annoying af sample suffixes
colnames(mapstats) = c('sample_id', 'locus', 'locus_length', 'mapped', 'unmapped')
mapstats$sample_id = gsub("\\_.*","",mapstats$sample_id)

#log mapped reads in case you want to view it in heatmap
mapstats$logmapped = log10(mapstats$mapped)
mapstats[mapstats == -Inf] <- 0 #replace inf with 0 so pheatmap doesn't moan

#join to mdata and sort by sample species composition
x_mapstats = left_join(mapstats, metadata, by = c('sample_id' = 'Sample'))
x_mapstats = x_mapstats[order(x_mapstats$Species),] 

#subset, create wideform table and write - change the selection to 'mapped' if you want raw mapped vs log
mapstatstable = x_mapstats %>% dplyr::select(sample_id, locus, logmapped) %>% pivot_wider(names_from = 'sample_id', values_from = 'logmapped')
write_csv(mapstatstable, "mapstatstable_mlst.csv")

#order rows by target species
mapstatstable$target_organism = targinfo$species
mapstatstable = mapstatstable[order(mapstatstable$target_organism),] 
mapstatstable = select(mapstatstable, -target_organism) #then drop

#make hm matrix
mat = as.matrix(mapstatstable[,2:ncol(mapstatstable)])
rownames(mat) = mapstatstable$locus

#create annotations
sample_annotations = metadata %>% select(Sample, Species) %>% distinct() %>% remove_rownames %>% column_to_rownames(var="Sample")
target_annotations = targinfo %>% select(name, species) %>% distinct() %>% remove_rownames %>% column_to_rownames(var="name")

#plot heapmap
pheatmap(mat,cluster_rows=FALSE, cluster_cols=FALSE, annotation_col = sample_annotations, annotation_row = target_annotations)




#ls pkg
pkg = c("tidyverse", "cowplot", "RColorBrewer", "cowplot", "pheatmap" , "data.table")
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#
#setwd to source dir (w/data files too)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#read our shit
metadata <- fread('../ancillary/metadata/mlst_metadata.csv') #read sample data
mapstats <- fread('../ancillary/metadata/idxstats.txt') #read long=form mapping data (no reads mapped/locus/sample)
zerocovstats <- fread('../ancillary/metadata/mlst_zerocov.txt')
#add colnames and remove polyomics annoying af sample suffixes
colnames(mapstats) = c('sample_id', 'locus', 'locus_length', 'mapped', 'unmapped')
mapstats$sample_id = gsub("\\_.*","",mapstats$sample_id)

#log mapped reads in case you want to view it in heatmap
mapstats$logmapped = log10(mapstats$mapped)
mapstats[mapstats == -Inf] <- 0 #replace inf with 0 so pheatmap doesn't moan

#create wideform table and write
mapstatstable = mapstats %>% dplyr::select(sample_id, locus, logmapped) %>% pivot_wider(names_from = 'sample_id', values_from = 'logmapped')
write_csv(mapstatstable, "../ancillary/metadata/mapstatstable_mlst.csv")

mat = as.matrix(mapstatstable[,2:ncol(mapstatstable)])
rownames(mat) = mapstatstable$locus
pheatmap(mat, kmeans_k = NA)

#read mlst zerocov table
zerocovstats$sample = gsub("_.*","",zerocovstats$id)
zerocovtab = pivot_wider(zerocovstats[,c(3,7,8)], names_from = sample, values_from = V5)
zerocovtab[is.na(zerocovtab)] <- 0
zerocovtab = zerocovtab[zerocovtab$V1 != 'genome',]
write_csv(zerocovtab, "../ancillary/metadata/zerocovtable_mlst.csv")

mat = as.matrix(zerocovtab[,2:ncol(zerocovtab)])
rownames(mat) = zerocovtab$V1
pheatmap(mat, kmeans_k = NA)


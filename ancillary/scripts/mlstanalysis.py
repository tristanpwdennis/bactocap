#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 11:34:37 2020

@author: tristanpwdennis
"""

import os
import pysam 
import io
import pandas as pd

bamdir = "/home/ubuntu/onchogenome/bactocap/datasets/mlst/mlst_only_cox_lep_bart_bcru"

#collect mapping data from bam file indexes
#takes a directory of indexed bam files and outputs samtools index stats (mapping per-target)
#using pysam to collect mapping stats per chromosome
#for each bam file, calculate coverage over each chromosome position - creating a df each time that is appended to 'cov'
#concatenate into one big dataframe with 'sample_id', ''target', 'position', 'depth'

idx = []
for entry in os.scandir(bamdir):
    if(entry.path.endswith(".bam")):
       bamtab = io.StringIO(pysam.idxstats(entry.path))
       bamtab = pd.read_table(bamtab, names = ['locus', 'sequence_length', 'mapped_reads', 'unmapped_reads'])
       bannamestring = os.path.splitext(os.path.basename(entry.path))[0]
       bamtab.insert(0, 'sample_id', bannamestring)
       idx.append(bamtab)

#final big df of per-target mapping stats for all samples 
idxframe = pd.concat(idx)

#export to tsv
idxframe.to_csv('/home/ubuntu/onchogenome/bactocap/datasets/mlst/mlst_only_cox_lep_bart_bcru/idxstats.txt', sep ='\t')

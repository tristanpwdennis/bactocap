#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 11:34:37 2020

@author: tristanpwdennis
"""

import os
import pandas as pd
import pysam 
import pybedtools as pbt
from io import StringIO

bamdir = "/Users/tristanpwdennis/Projects/bactocap/datasets/mlst/results/bam"

#collect mapping data from bam file indexes
#using pysam to collect mapping stats per chromosome
#for each bam file, calculate coverage over each chromosome position - creating a df each time that is appended to 'cov'
#concatenate into one big dataframe with 'sample_id', ''target', 'position', 'depth'
idx = []
for entry in os.scandir(bamdir):
    if(entry.path.endswith(".bam")):
       bamtab = StringIO(pysam.idxstats(entry.path))
       bamtab = pd.read_table(bamtab, names = ['locus', 'sequence_length', 'mapped_reads', 'unmapped_reads'])
       bannamestring = os.path.splitext(os.path.basename(entry.path))[0]
       bamtab.insert(0, 'sample_id', bannamestring)
       idx.append(bamtab)
       
#collect list of dataframes into a single dataframe       
idxstats = pd.concat(idx, axis=0, ignore_index=True)

#collect coverage data from bam files
#for each bam file, calculate coverage over each chromosome position - creating a df each time that is appended to 'cov'
#concatenate into one big dataframe with 'sample_id', ''target', 'position', 'depth'
cov = []
for entry in os.scandir(bamdir):
    if(entry.path.endswith(".bam")):
        bamobj = pbt.BedTool(entry.path)
        covstats = bamobj.genome_coverage(d=True, **{'5': True})
        t = pd.read_table(covstats.fn, names=['chrom', 'pos', 'dp'])
        bannamestring = os.path.splitext(os.path.basename(entry.path))[0]
        t.insert(0, 'sample_id', bannamestring)
        cov.append(t)
        
covs =  pd.concat(cov, axis=0, ignore_index=True)      


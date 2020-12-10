#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 09:52:03 2020

@author: tristanpwdennis
"""

import pybedtools
import pandas as pd
import numpy as np
import glob
import os

org = "mycoplasma"
top = "/home/ubuntu/onchogenome/bactocap/datasets/"
sub = "/results/*/*.bam"

path = top + org + sub
bams= glob.glob(path)

if(org == "mycoplasma"):
    bed = pybedtools.BedTool('/home/ubuntu/onchogenome/bactocap/ancillary/annotation/aln.Myco.4columns.bed')
elif(org == "anthrax"):
    bed = pybedtools.BedTool('/home/ubuntu/onchogenome/bactocap/ancillary/annotation/aln_BaBaits.4columns.bed')
  
    #extract per-position coverage using the bait bed file (wraps bedtools)
def getcovdata(bams, bed):
    testdf = []
    for bam in bams:
        organism = org
        testoutput = bed.coverage(pybedtools.BedTool(bam), d=True)
        df = testoutput.to_dataframe(names=['genome', 'start' 'finish', 'bait', 'pos', 'cov'])
        df['sample'] = os.path.basename(bam)
        df['organism'] = organism
        testdf.append(df)    
    testdf = pd.concat(testdf)
    return(testdf)

#get big covdf - for all samples
output = getcovdata(bams, bed)

#count unique baits where cov  at least in one position is less than the threshold
v=[]
for cov in np.arange(0,200,5):
   df = output[(output['cov'] < cov)].groupby('sample')['bait'].nunique().to_frame().reset_index() 
   df['cov'] = cov
   v.append(df)

#concat to big df
t = pd.concat(v)
#export csv
t.DataFrame.to_csv('/home/ubuntu/onchogenome/bactocap/covdf.csv')

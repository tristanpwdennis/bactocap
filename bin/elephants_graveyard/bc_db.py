#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 09:52:03 2020

@author: tristanpwdennis
"""

import pybedtools
import numpy as np
import glob
import os
import sqlite3
import time
import itertools as it
import os
import cython
import pandas as pd

org = "mycoplasma"
top = "/Users/tristanpwdennis/Projects/bactocap/datasets/"
sub = "/results/*/*.bam"

path = top + org + sub
bams= glob.glob(path)

if(org == "mycoplasma"):
    bed = pybedtools.BedTool('/Users/tristanpwdennis/Projects/bactocap/ancillary/annotation/aln.Myco.4columns.bed')
elif(org == "anthrax"):
    bed = pybedtools.BedTool('/Users/tristanpwdennis/Projects/bactocap/ancillary/annotation/aln_BaBaits.4columns.bed')

bam = "/Users/tristanpwdennis/Projects/bactocap/datasets/mycoplasma/results/09B10514_S83/09B10514_S83.rg.bam"

####is it faster to just return a histogram, or the per-base?
start_time = time.time() 
covhist = bed.coverage(pybedtools.BedTool(bam))
print("--- %s seconds ---" % (time.time() - start_time)) 

start_time = time.time() 
covperbase = bed.coverage(pybedtools.BedTool(bam), d=True)
print("--- %s seconds ---" % (time.time() - start_time))

#let's see roughly how big each is
len(covhist) #24444 rows
#almost 2 million rows
len(covperbase)



    #extract per-position coverage using the bait bed file (wraps bedtools)
def getcovdata(bams, bed):
    testdf = []
    for bam in bams:
        testoutput = bed.coverage(pybedtools.BedTool(bam), d=True)
        df = testoutput.to_dataframe(names=['genome', 'start' 'finish', 'bait', 'pos', 'cov'])
        df['sample'] = os.path.basename(bam)
        testdf.append(df)    
    testdf = pd.concat(testdf)
    return(testdf)

#get big covdf - for all samples
output = getcovdata(bams, bed)


#count unique baits where cov  at least in one position is less than the threshold

def getbaitstocov(bedtoolsoutput):
    v=[]
    for cov in np.arange(0,200,5):
        df = bedtoolsoutput[(bedtoolsoutput['cov'] < cov)].groupby('sample')['bait'].nunique().to_frame().reset_index() 
        df['cov'] = cov
        v.append(df)
    return(v)

#import pandas as pd
t = getbaitstocov(covperbase)

#concat to big df
t = pd.concat(t)
#export big csv
t.to_csv('/Users/tristanpwdennis/Projects/bactocap/covdf.csv')

















covdb = sqlite3.connect('cov.db')
c = covdb.cursor()
#c.execute('''CREATE TABLE covtable
#             (genome text, startfinish text, bait text, pos real, cov real, sample text, organism text)''')

             

#test db import 
df = testoutput.to_dataframe(names=['genome', 'start' 'finish', 'bait', 'pos', 'cov']).reset_index(drop=True) 
df['sample'] = os.path.basename('/Users/tristanpwdennis/Projects/bactocap/datasets/mycoplasma/results/09B10127_S48/09B10127_S48.rg.bam')
df['organism'] = org
df.head()
df.info()
df.to_sql('covtable',covdb, if_exists = 'append')







#let's try speed things up a abit
#create test set
d = []
for s in (it.product(['sample1', 'sample2'], ['bait1', 'bait2', 'bait3', 'bait4'], np.arange(0,80,1))):
    d.append(s)
testoutput = pd.DataFrame(d)
testoutput.columns = ['sample', 'bait', 'pos']
testoutput['cov'] = np.random.choice([1, 9, 20], testoutput.shape[0])
testoutput['genome'] = "genome"
testoutput['startfinish'] = 1

#try the regular#full-fat version
start_time = time.time() 
t = getbaitstocov(output)
print("--- %s seconds ---" % (time.time() - start_time))    


#take a slice of the df containing only columns we care about
#this is a tiny bit faster but not much - maybe it will reduce the amount we put into memory at least
slice = testoutput[['sample', 'bait', 'cov']]
start_time = time.time() 
t = getbaitstocov(slice)
print("--- %s seconds ---" % (time.time() - start_time))  


slice[slice.groupby("bait")["cov"].transform('min') == slice['cov']]

start_time = time.time() 
t=output.loc[output.groupby('bait')['cov'].idxmin()]
print("--- %s seconds ---" % (time.time() - start_time))  






v=[]
for cov in np.arange(0,200,5):
    df = test[(test['cov'] < cov)].groupby('sample')['bait'].nunique().to_frame().reset_index() 
    df['cov'] = cov
    v.append(df)






s = output[['sample', 'bait', 'cov']]
s['baitid'] = s['sample'] + '/' + s['bait']
s = s.drop(['sample', 'bait'], axis=1)
t = s.loc[s.groupby('baitid')['cov'].idxmin()]




def getbaitstocovmod(bedtoolsoutput):
    v=[]
    for cov in np.arange(0,200,5):
        df = bedtoolsoutput[(bedtoolsoutput['cov'] < cov)].groupby('baitid')['baitid'].nunique().to_frame().reset_index() 
        df['cov'] = cov
        v.append(df)
    return(v)

start_time = time.time() 
t  =getbaitstocovmod(s)
print("--- %s seconds ---" % (time.time() - start_time))  




             
testoutput = bed.coverage(pybedtools.BedTool('/Users/tristanpwdennis/Projects/bactocap/datasets/mycoplasma/results/09B10127_S48/09B10127_S48.rg.bam'), d=True)

#convert to df like this
test = testoutput.to_dataframe()

#
d = dict(testoutput)
    


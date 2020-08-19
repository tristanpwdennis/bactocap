import numpy as np
import pandas as pd
import itertools as it
import plotnine as pn

#import data
# turn this into a loop or function sometime
andata = pd.read_csv("/Users/tristanpwdennis/Projects/bactocap/datasets/anthrax/results/covstatstotal.txt", sep=' ') 
mydata = pd.read_csv("/Users/tristanpwdennis/Projects/bactocap/datasets/mycoplasma/covstatstotal.txt", sep=' ')
mymeta = pd.read_csv("/Users/tristanpwdennis/Projects/bactocap-data/mycoplasma/metadata/myco-seq-sample-data.csv")
anmeta = pd.read_csv("/Users/tristanpwdennis/Projects/bactocap-data/anthrax/metadata/anthrax-metadata-1.csv")

#get sampleset identifiers (in this case the number of samples in the set)
#and then get fraction of baits captured to our degree of coverage as a fraction of the total baits for each sampleset
mycolval = pd.Series([len(mymeta.index) for i in range(len(mydata.index))])
mydata.insert(loc = 0, column = 'ninds', value = mycolval)
mydata['frac'] = mydata['num_baits']/24444
#repeat
ancolval = pd.Series([len(anmeta.index) for i in range(len(andata.index))])
andata.insert(loc = 0, column = 'ninds', value = ancolval)
andata['frac'] = andata['num_baits']/148729

df = pd.concat([andata, mydata])

size_df = pd.DataFrame(columns=["frac", "coverage", "ninds", "size"])
for f, c, i in it.product(np.arange(0,1.1,0.1), np.arange(0,410,10),[113, 56]):
    mask = (df.frac > f) & (df.coverage == c) & (df.ninds == i)
    size = mask.sum()
    size_df.loc[len(size_df)] = [f, c, i, size]
    
size_df['frac_baits'] = size_df['size']/size_df['ninds']


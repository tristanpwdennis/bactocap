import numpy as np
import pandas as pd
import itertools as it
import plotnine as p9
#import data
# turn this into a loop or function sometime
andata = pd.read_csv("/Users/tristanpwdennis/Projects/bactocap/datasets/anthrax/results/covstatstotal.txt", sep=' ', names=["num_baits", "sample", "depth"]) 
mydata = pd.read_csv("/Users/tristanpwdennis/Projects/bactocap/datasets/mycoplasma/covstatstotal.txt", sep=' ', names=["num_baits", "sample", "depth"])
mymeta = pd.read_csv("/Users/tristanpwdennis/Projects/bactocap-data/mycoplasma/metadata/myco-seq-sample-data.csv")
anmeta = pd.read_csv("/Users/tristanpwdennis/Projects/bactocap-data/anthrax/metadata/anthrax-metadata-1.csv")

#get sampleset identifiers (in this case the number of samples in the set)
#and then get fraction of baits captured to our degree of coverage as a fraction of the total baits for each sampleset
def wrangle(datafile, metadatafile, numbaits): 
    colval = pd.Series([len(metadatafile.index) for i in range(len(datafile.index))])
    datafile.insert(loc = 0, column = 'ninds', value = colval)
    datafile['frac'] = datafile['num_baits']/numbaits
    return(datafile)
    
#bind together
df = pd.concat([wrangle(mydata, mymeta, 24444), wrangle(andata, anmeta, 148729)])

df['organism'] = np.where(df['ninds']==113, 'anthrax', 'mycoplasma')

#intialise empty dataframe
size_df = pd.DataFrame(columns=["frac", "coverage", "ninds", "size"])

#we want to calculate the number of individuals where we covered the genome to a depth of 0X to 400X, and a breadth of 0-1 (as fraction)
#so we create the cartesian product of these ranges, expand them into a triplet in a for loop and iterate over them all (3-way nested loop)
#the count the number of rows matching each combination of conditions, populating a dataframe with the result
for f, c, i in it.product(np.arange(0.1,1.0,0.1), np.arange(0,410,10),[113, 56]):
    mask = (df.frac > f) & (df.depth == c) & (df.ninds == i)
    size = mask.sum()
    size_df.loc[len(size_df)] = [f, c, i, size]
    
#get the fraction of individuals per sample set by dividing the number of individuals matching each condition
size_df['frac_inds'] = size_df['size']/size_df['ninds']
#add organism name for ninds as above (boring and redundant)
size_df['organism'] = np.where(size_df['ninds']==113, 'anthrax', 'mycoplasma')

#plot depth of coverage to which an increasing fraction of baits are totally covered in both
#sample sets

(
 p9.ggplot(size_df, p9.aes(x = 'coverage', y = 'frac_inds', colour = 'organism'))+
  p9.facet_wrap('frac') +
  p9.geom_line() +
  p9.xlab("Depth of Coverage") +
  p9.ylab("Fraction of Individuals") +
  p9.theme_minimal()
 )

#plot fraction of baits covered to x depth in all individuals for both anth and myco
(
 p9.ggplot(df, p9.aes(x='depth', y ='frac', colour = 'sample')) +
     p9.geom_line() +
     p9.facet_wrap('organism') +
     p9.scale_color_discrete(guide=False) +
     p9.theme_minimal() 
 )


import numpy as np
import pandas as pd
import itertools as it
from plotnine import *

#import data
# turn this into a loop or function sometime
andata = pd.read_csv("/Users/tristanpwdennis/Projects/bactocap/datasets/anthrax/results/covstatstotal.txt", sep=' ', names=["num_baits", "sample", "depth"]) 
mydata = pd.read_csv("/Users/tristanpwdennis/Projects/bactocap/datasets/mycoplasma/covstatstotal.txt", sep=' ', names=["num_baits", "sample", "depth"])
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
#bind together
df = pd.concat([andata, mydata])
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


test = size_df[(size_df.ninds == 56)]

size_df['organism'] = np.where(size_df['ninds']==113, 'anthrax', 'mycoplasma')



#plot depth of coverage to which an increasing fraction of baits are totally covered in both
#sample sets

(
 ggplot(size_df, aes(x = 'coverage', y = 'frac_inds', colour = 'organism'))+
  facet_wrap('frac') +
  geom_line() +
  xlab("Depth of Coverage") +
  ylab("Fraction of Individuals") +
  theme_minimal()
 )

#plot fraction of baits covered to x depth in all individuals for both anth and myco
(
 ggplot(df, aes(x='depth', y ='frac', colour = 'sample')) +
     geom_line() +
     facet_wrap('organism') +
     scale_color_discrete(guide=False) +
     theme_minimal() 
 )


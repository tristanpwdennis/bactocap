#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 11:34:37 2020

@author: tristanpwdennis
"""

import os
import numpy as np
import pandas as pd
import pysam 
from io import StringIO

bamdir = "/Users/tristanpwdennis/Projects/bactocap/datasets/mlst/results/bam"

li = []

for entry in os.scandir(bamdir):
    if(entry.path.endswith(".bam")):
       bamtab = StringIO(pysam.idxstats(entry.path))
       bamtab = pd.read_table(bamtab, names = ['locus', 'sequence_length', 'mapped_reads', 'unmapped_reads'])
       bannamestring = os.path.splitext(os.path.basename(entry.path))[0]
       bamtab.insert(0, 'bamname', bannamestring)
       li.append(bamtab)
       
idxstats = pd.concat(li, axis=0, ignore_index=True)

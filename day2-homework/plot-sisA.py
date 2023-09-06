#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Get dataset to recreate Fig 3B from Lott et al 2011 PLoS Biology https://pubmed.gov/21346796
# wget https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/bulk_RNA-seq/extra_data/all_annotated.csv

transcripts = np.loadtxt( "all_annotated.csv", delimiter=",", usecols=0, dtype="<U30", skiprows=1 )
print( "transcripts: ", transcripts[0:5] )

samples = np.loadtxt( "all_annotated.csv", delimiter=",", max_rows=1, dtype="<U30" )[2:]
print( "samples: ", samples[0:5] )

data = np.loadtxt( "all_annotated.csv", delimiter=",", dtype=np.float32, skiprows=1, usecols=range(2, len(samples) + 2) )
print( "data: ", data[0:5, 0:5] )

# Find row with transcript of interest
for i in range(len(transcripts)):
    if transcripts[i] == 'FBtr0073461':
        row = i

# Find columns with samples of interest
cols = []
for i in range(len(samples)):
    if "female" in samples[i]:
        cols.append(i)
#male columns
cols_m = []
for i in range(len(samples)):
    if "male" in samples[i] and "female" not in samples[i]:
        cols_m.append(i)

# Subset data of interest
expression = data[row, cols]
expression_m = data[row, cols_m]
# Prepare data
x = samples[cols]
y = expression
x_m = samples[cols_m]
y_m = expression_m
# Plot data
fig, ax = plt.subplots()
ax.set_title( "FBtr0073461" )
ax.plot( x, y )
ax.plot(x_m, y_m )
fig.savefig( "FBtr0073461.png" )
plt.close( fig )

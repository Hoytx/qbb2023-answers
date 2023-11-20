#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sulfite = pd.read_table('ONT_Data/bisulfite.cpg.chr2.bedgraph', sep = '\t', names=["CHR", "POS1", "POS2", "PERCENT_sulfite", "COVERAGE_sulfite"])
onp = pd.read_table('ONT_Data/ONT.cpg.chr2.bedgraph', sep = '\t', names=["CHR", "POS1", "POS2", "PERCENT_onp", "COVERAGE_onp"])

sulfite_pos = sulfite["POS1"]
onp_pos = onp["POS1"]

overlap_pos = set(sulfite_pos).intersection(onp_pos)
only_sulfite = np.setdiff1d(onp_pos,sulfite_pos)
only_onp = np.setdiff1d(sulfite_pos, onp_pos)
overlap_count = len(overlap_pos)
tot_pos = len(overlap_pos) + len(only_sulfite) + len(only_onp)
unique_pos = len(only_sulfite) + len(only_onp)


sx = overlap_count/tot_pos
nx = len(only_onp)/tot_pos
bx = len(only_sulfite)/tot_pos
print(f"Percentage of shared sites: {sx}")
print(f"Percentage of nanopore unique sites: {nx}")
print(f"Percentage of bisulfite unique sites: {bx}")

#make list of coverages in samples and graph them

sulfite_coverage = sulfite["COVERAGE_sulfite"]
onp_coverage = onp["COVERAGE_onp"]

fig, (ax, ay, az) = plt.subplots(3, 1)

ax.hist([sulfite_coverage, onp_coverage], alpha = 0.4, bins = 100, range = (0, 100))
ax.legend(["Bismark", "Oxford Nanopore"])
ax.set_xlabel("Coverage Depth")
ax.set_ylabel("# of sites")

#make new dataframes with overlapping position's data
sulfite_shared = sulfite.loc[sulfite["POS1"].isin(overlap_pos)]
onp_shared = onp.loc[onp["POS1"].isin(overlap_pos)]
sulfite_and_onp = pd.merge(sulfite_shared,onp_shared[['POS1', 'PERCENT_onp', 'COVERAGE_onp']],on='POS1', how='left')
sulfite_shared.reset_index(drop=True, inplace=True)
onp_shared.reset_index(drop=True, inplace=True)
sulfite_and_onp.reset_index(drop=True, inplace=True)

#select data and contruct 2d histogram
x_hist = sulfite_and_onp["PERCENT_sulfite"]
y_hist = sulfite_and_onp["PERCENT_onp"]

relationship, xedges, yedges = np.histogram2d(x_hist,y_hist, bins=100)
transformed = np.log10(relationship + 1)
print(np.corrcoef(x_hist,y_hist))
ay.imshow(transformed)
ay.set_title("Pearson R coeff = 0.936", fontsize = 8)
ay.set_xlabel("Bismark")
ay.set_ylabel("Oxford Nanopore")

#load matched samples
normal_onp = pd.read_table('ONT_Data/normal.ONT.chr2.bedgraph', sep = '\t', names=["CHR", "POS1", "POS2", "PERCENT_normal", "COVERAGE_normal"])
tumor_onp = pd.read_table('ONT_Data/tumor.ONT.chr2.bedgraph', sep = '\t', names=["CHR", "POS1", "POS2", "PERCENT_tumor", "COVERAGE_tumor"])
normal_sulfite = pd.read_table('ONT_Data/normal.bisulfite.chr2.bedgraph', sep = '\t', names=["CHR", "POS1", "POS2", "PERCENT_normal", "COVERAGE_normal"])
tumor_sulfite = pd.read_table('ONT_Data/tumor.bisulfite.chr2.bedgraph', sep = '\t', names=["CHR", "POS1", "POS2", "PERCENT_tumor", "COVERAGE_tumor"])

#find overlapping sites for each method
onp_overlap = set(normal_onp["POS1"]).intersection(tumor_onp["POS1"])
normal_onp_shared = normal_onp.loc[normal_onp["POS1"].isin(onp_overlap)]
tumor_onp_shared = tumor_onp.loc[tumor_onp["POS1"].isin(onp_overlap)]
sulfite_overlap = set(normal_sulfite["POS1"]).intersection(tumor_sulfite["POS1"])
normal_sulfite_shared = normal_sulfite.loc[normal_sulfite["POS1"].isin(sulfite_overlap)]
tumor_sulfite_shared = tumor_sulfite.loc[tumor_sulfite["POS1"].isin(sulfite_overlap)]

onp_overlap_data = pd.merge(normal_onp_shared, tumor_onp_shared[['POS1', 'PERCENT_tumor', 'COVERAGE_tumor']],on='POS1', how='left')
sulfite_overlap_data = pd.merge(normal_sulfite_shared,tumor_sulfite_shared[['POS1', 'PERCENT_tumor', 'COVERAGE_tumor']],on='POS1', how='left')

#calculate difference in methylation only for sites shared between methods
method_overlap_positions = set(onp_overlap_data["POS1"]).intersection(sulfite_overlap_data["POS1"])
shared_site_onp = onp_overlap_data.loc[onp_overlap_data["POS1"].isin(method_overlap_positions)]
shared_site_sulfite = sulfite_overlap_data.loc[sulfite_overlap_data["POS1"].isin(method_overlap_positions)]
shared_site_difference_onp = np.subtract(shared_site_onp["PERCENT_tumor"], shared_site_onp["PERCENT_normal"])
shared_site_difference_sulfite = np.subtract(shared_site_sulfite["PERCENT_tumor"], shared_site_sulfite["PERCENT_normal"])
#calculate pearson R coeff and print the result
print(np.corrcoef(shared_site_difference_onp,shared_site_difference_sulfite))

#resetting index to allow for subtracting each item correctly
normal_onp_shared.reset_index(drop=True, inplace=True)
tumor_onp_shared.reset_index(drop=True, inplace=True)
normal_sulfite_shared.reset_index(drop=True, inplace=True)
tumor_sulfite_shared.reset_index(drop=True, inplace=True)

onp_difference = np.subtract(tumor_onp_shared["PERCENT_tumor"], normal_onp_shared["PERCENT_normal"])
onp_difference = onp_difference[onp_difference != 0]
onp_difference.reset_index(drop=True, inplace=True)

sulfite_difference = np.subtract(tumor_sulfite_shared["PERCENT_tumor"], normal_sulfite_shared["PERCENT_normal"])
sulfite_difference = sulfite_difference[sulfite_difference != 0]
sulfite_difference.reset_index(drop=True, inplace=True)


az.violinplot([onp_difference, sulfite_difference], showmeans = True)
az.set_xticks([1,2])
az.set_xticklabels(["Oxford Nanopore", "Bismark"])
az.set_ylabel("Methylation Difference")
az.set_title("Pearson R coeff = 1")

fig.tight_layout()
fig.savefig("Plot.png", dpi=1000)
plt.show()
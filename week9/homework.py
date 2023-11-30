#!/usr/bin/env python

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats import multitest
from pydeseq2 import preprocessing
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt

# read in data
counts_df = pd.read_csv("gtex_whole_blood_counts_formatted.txt", index_col = 0)

# read in metadata
metadata = pd.read_csv("gtex_metadata.txt", index_col = 0)

#normalization
counts_df_normed = preprocessing.deseq2_norm(counts_df)[0]
#log2 of normalization
counts_df_normed = np.log2(counts_df_normed + 1)
#Add the metadata
full_design_df = pd.concat([counts_df_normed, metadata], axis=1)

#run linear regression for a gene
model = smf.ols(formula = 'Q("DDX11L1") ~ SEX', data=full_design_df)
results = model.fit()

slope = results.params[1]
pval = results.pvalues[1]
shape = counts_df_normed.shape
sex_difference = pd.DataFrame(columns = ["Gene", "Slope", "Pvalue"])
gene = ""
#run through the dataframe and compare expression between sexes for every gene
'''
for i in range(shape[1]):
    gene = full_design_df.columns[i]
    if gene == "SUBJECT_ID":
        continue
    else:
        model = smf.ols(formula = f'Q("{gene}") ~ SEX', data=full_design_df)
        results = model.fit()
        slope = results.params[1]
        pval = results.pvalues[1]
        sex_difference.loc[len(sex_difference.index)] = [gene, slope, pval]
#filling in blank pvalues
sex_difference['Pvalue'] = sex_difference['Pvalue'].fillna(1.0)
#write to file
sex_difference.to_csv("Sex_Difference.csv")
'''
load = pd.read_csv("Sex_Difference.csv")
load["rejected"], load["pvalue-corrected"] = multitest.fdrcorrection(load["Pvalue"], alpha = 0.1)
load.to_csv("Sex_Difference.csv")

FDR_hits = load.loc[load["pvalue-corrected"] <= 0.1]
FDR_hits.to_csv("FDR_corrected_step1-5.csv")

manual_hits = FDR_hits["Gene"]

file1 = open("FDR_corrected_step1-5.txt", "r+")
file1.writelines(manual_hits)

dds = DeseqDataSet(counts=counts_df, metadata=metadata, design_factors="SEX", n_cpus=4)

dds.deseq2()
stat_res = DeseqStats(dds)
stat_res.summary()
results = stat_res.results_df
#saving extra copy for use in volcano plot
unfiltered_results = stat_res.results_df
results = results.loc[results["padj"] <= 0.1]
results["Gene"] = results.index

dds_genes = results["Gene"]

file2 = open("PyDESeq2_result_genes.txt", "r+")
file2.writelines(dds_genes)


results.to_csv("PyDESeq2_result_genes.csv")

manual_10FDR = manual_hits
dds_10FDR = dds_genes.tolist()

#print(manual_10FDR)
#print(dds_10FDR)

in_both = list(set(manual_10FDR).intersection(dds_10FDR))
in_either = list(manual_10FDR) + list(dds_10FDR)
in_either = set(in_either)
#print(len(in_both))
#print(len(in_either))
overlap_percentage = (len(in_both) / len(in_either)) * 100
print(f"Overlap is {overlap_percentage} %")

file1.close()
file2.close()

#volcano plot of results

fig, ax = plt.subplots()
col = []
for i in range(results.shape[0]):
    if results["log2FoldChange"][i] > 1 and results["padj"][i] <= 0.1:
        col.append("red")
    else:
        col.append("gray")

y = -(np.log10(results["padj"] + 1e-250))

ax.scatter(results["log2FoldChange"], y, c=col, s=2 )
ax.set_xlabel("log2FoldChange")
ax.set_ylabel("-log10(q-value)")

fig.savefig("Volcano_Plot.png", dpi=1000)
plt.show()

#!/usr/bin/env python

import numpy as np
import pandas as pd
from pydeseq2 import preprocessing
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# read in data
counts_df = pd.read_csv("gtex_whole_blood_counts_formatted.txt", index_col = 0)

# read in metadata
metadata = pd.read_csv("gtex_metadata.txt", index_col = 0)

# normalize
counts_df_normed = preprocessing.deseq2_norm(counts_df)[0]

# log
counts_df_logged = np.log2(counts_df_normed + 1)

#1 male 2 female in metadata
# merge with metadata
full_design_df = pd.concat([counts_df_logged, metadata], axis=1)

subject_exp = counts_df_logged.loc["GTEX-113JC"]


fig, ax = plt.subplots()

ax.hist(subject_exp, bins=40, range=[0, 20], color="green")
ax.set_xlabel("log2 Normed Expression")
ax.set_ylabel("# of Genes")
ax.set_xticks(ticks=range(0,21))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
fig.savefig("GTEX-113JC_Gene_Expression.png")

fig2, ay = plt.subplots()
single_gene = full_design_df[["MXD4", "SEX"]]
male_exp = single_gene.loc[single_gene['SEX'] == 1]
female_exp = single_gene.loc[single_gene['SEX'] == 1]
ay.boxplot([male_exp["MXD4"], female_exp["MXD4"]])
ay.set_xticks(ticks=[1,2], labels=["Male","Female"])
ay.set_ylabel("log2 Normed MXD4 Expression")

fig2.savefig("MXD4_Expression_By_Sex.png")

fig3, az = plt.subplots()
ages = [full_design_df["AGE"].value_counts()["20-29"], full_design_df["AGE"].value_counts()["30-39"], full_design_df["AGE"].value_counts()["40-49"], full_design_df["AGE"].value_counts()["50-59"], full_design_df["AGE"].value_counts()["60-69"], full_design_df["AGE"].value_counts()["70-79"]]
#print(ages)
az.bar(["20-29","30-39","40-49","50-59","60-69","70-79"], height = ages, color ="green")
az.set_ylabel("# of Subjects")
az.set_xlabel("Ages")
fig3.savefig("Subjects_by_Age_Catagory.png")

fig4, av =plt.subplots()
LPXN_data = full_design_df[["LPXN", "SEX", "AGE"]]

male_LPXN = []
male_LPXN_std = []
female_LPXN = []
female_LPXN_std = []

#filter the log normalized counts based on sex and age catagory and save results as lists within a master list
male_LPXN.append(LPXN_data["LPXN"][(LPXN_data['SEX'] == 1) & (LPXN_data["AGE"] == "20-29")].values.tolist())
male_LPXN.append(LPXN_data["LPXN"][(LPXN_data['SEX'] == 1) & (LPXN_data["AGE"] == "30-39")].values.tolist())
male_LPXN.append(LPXN_data["LPXN"][(LPXN_data['SEX'] == 1) & (LPXN_data["AGE"] == "40-49")].values.tolist())
male_LPXN.append(LPXN_data["LPXN"][(LPXN_data['SEX'] == 1) & (LPXN_data["AGE"] == "50-59")].values.tolist())
male_LPXN.append(LPXN_data["LPXN"][(LPXN_data['SEX'] == 1) & (LPXN_data["AGE"] == "60-69")].values.tolist())
male_LPXN.append(LPXN_data["LPXN"][(LPXN_data['SEX'] == 1) & (LPXN_data["AGE"] == "70-79")].values.tolist())


#find the average for each age catagory and store the results as a new list
male_LPXN_avg = [(sum(male_LPXN[0])/len(male_LPXN[0])), (sum(male_LPXN[1])/len(male_LPXN[1])), (sum(male_LPXN[2])/len(male_LPXN[2])), (sum(male_LPXN[3])/len(male_LPXN[3])), (sum(male_LPXN[4])/len(male_LPXN[4])), (sum(male_LPXN[5])/len(male_LPXN[5]))]

#find standard deviation for each age group

male_LPXN_std.append(np.std(male_LPXN[0]))
male_LPXN_std.append(np.std(male_LPXN[1]))
male_LPXN_std.append(np.std(male_LPXN[2]))
male_LPXN_std.append(np.std(male_LPXN[3]))
male_LPXN_std.append(np.std(male_LPXN[4]))
male_LPXN_std.append(np.std(male_LPXN[5]))

#repeat above but for female data
female_LPXN.append(LPXN_data["LPXN"][(LPXN_data['SEX'] == 2) & (LPXN_data["AGE"] == "20-29")].values.tolist())
female_LPXN.append(LPXN_data["LPXN"][(LPXN_data['SEX'] == 2) & (LPXN_data["AGE"] == "30-39")].values.tolist())
female_LPXN.append(LPXN_data["LPXN"][(LPXN_data['SEX'] == 2) & (LPXN_data["AGE"] == "40-49")].values.tolist())
female_LPXN.append(LPXN_data["LPXN"][(LPXN_data['SEX'] == 2) & (LPXN_data["AGE"] == "50-59")].values.tolist())
female_LPXN.append(LPXN_data["LPXN"][(LPXN_data['SEX'] == 2) & (LPXN_data["AGE"] == "60-69")].values.tolist())
female_LPXN.append(LPXN_data["LPXN"][(LPXN_data['SEX'] == 2) & (LPXN_data["AGE"] == "70-79")].values.tolist())

female_LPXN_avg = [(sum(female_LPXN[0])/len(female_LPXN[0])), (sum(female_LPXN[1])/len(female_LPXN[1])), (sum(female_LPXN[2])/len(female_LPXN[2])), (sum(female_LPXN[3])/len(female_LPXN[3])), (sum(female_LPXN[4])/len(female_LPXN[4])), (sum(female_LPXN[5])/len(female_LPXN[5]))]

female_LPXN_std.append(np.std(female_LPXN[0]))
female_LPXN_std.append(np.std(female_LPXN[1]))
female_LPXN_std.append(np.std(female_LPXN[2]))
female_LPXN_std.append(np.std(female_LPXN[3]))
female_LPXN_std.append(np.std(female_LPXN[4]))
female_LPXN_std.append(np.std(female_LPXN[5]))



av.errorbar([1.1,2.1,3.1,4.1,5.1,6.1], male_LPXN_avg, yerr = male_LPXN_std, alpha = 0.4, capsize = 5, fmt = "o")
av.errorbar([0.9,1.9,2.9,3.9,4.9,5.9], female_LPXN_avg, yerr = female_LPXN_std, alpha = 0.4, capsize = 5, fmt = "s")
av.legend(labels=["Male", "Female"], loc = "lower right")
av.set_xticklabels(["","20-29","30-39","40-49","50-59","60-69","70-79"])
av.set_xlabel("Age")
av.set_ylabel("Median log2 Normed Expression")
fig4.savefig("LPXN_Expression_By_Age_&_Sex.png")

plt.show()
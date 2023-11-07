#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

pc1 = []
pc2 = []

for line in open('genotypes.eigenvec'):
    if line.startswith('FID'):
        continue
    fields = line.rstrip('\n').split(' ')

    pc1.append(float(fields[2]))
    pc2.append(float(fields[3]))

fig, ax = plt.subplots()

ax.scatter(pc1,pc2, s=1)
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
fig.savefig("PCA_plot.png", dpi=1000)
plt.close()

allele_frequencies = []

for line in open('plink.frq'):
    if line.find('CHR') != -1:
        continue
    fields = line.rstrip('\n').split()
    allele_frequencies.append(float(fields[4]))

#convert frequency to percent scale for graphing
allele_percents = []
for i in range(len(allele_frequencies)):
    x = float(allele_frequencies[i]) * 100
    allele_percents.append(x)

fig, ax = plt.subplots()
plt.xticks(fontsize = 6,rotation=60)
plt.yticks(fontsize = 6)

ax.hist(allele_percents, bins = range(0, 100, 5))
ax.set_xlabel("Allelic Frequency %")
ax.set_ylabel("# of Variants")
fig.savefig("AFS.png", dpi=1000)
plt.close()

######
#Adapted from code I found here: https://piperwrites95180714.wordpress.com/2018/04/04/genome-wide-association-study-manhattan-plot-tutorial/
######
#Commented out to avoid running over and over as I code


pd.set_option('expand_frame_repr', False)
df = pd.read_table('GS451_IC50_gwas_results.assoc.linear', sep = '\s+')
#print(df.head(10))

df['p_adj'] = -np.log10(df['P'])
df['CHR'] = df['CHR'].astype('category')

#print(df.head(10))

df['ind'] = range(len(df))
df_filtered = df[df['TEST'] == 'ADD']
df_grouped = df_filtered.groupby(('CHR'))
#print(df_grouped.head(10))

fig, (ax, ay) = plt.subplots(2, 1)
cmap, norm = mcolors.from_levels_and_colors([0,5], ['gray', 'red'], extend = 'max')
x_labels = []
x_labels_pos = []
chromosome = 0
for num, (name, group) in enumerate(df_grouped):
    chromosome += 1
    vals = df_filtered.loc[df_filtered['CHR'] == chromosome, 'p_adj'].to_list()
    group.plot(kind='scatter', x='ind', y='p_adj', c = vals, cmap = cmap, norm = norm, ax = ax, s=4)
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))

ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.tick_params(labelsize = 6, rotation = 60)
ax.set_xlim([0, len(df)])
ax.set_ylim([0, 15])
ax.set_xlabel('Chromosome')
ax.set_ylabel('-log10(p_value)')
ax.set_title('GS451_IC50_gwas_results')

df = pd.read_table('CB1908_IC50_gwas_results.assoc.linear', sep = '\s+')
#print(df.head(10))

df['p_adj'] = -np.log10(df['P'])
df['CHR'] = df['CHR'].astype('category')

#print(df.head(10))

df['ind'] = range(len(df))
df_filtered = df[df['TEST'] == 'ADD']
df_grouped = df_filtered.groupby(('CHR'))
#print(df_grouped.head(10))

x_labels = []
x_labels_pos = []
chromosome = 0
for num, (name, group) in enumerate(df_grouped):
    chromosome += 1
    vals = df_filtered.loc[df_filtered['CHR'] == chromosome, 'p_adj'].to_list()
    group.plot(kind='scatter', x='ind', y='p_adj', c = vals, cmap = cmap, norm = norm, ax = ay, s=4)
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))

ay.set_xticks(x_labels_pos)
ay.set_xticklabels(x_labels)
ay.tick_params(labelsize = 6, rotation = 60)
ay.set_xlim([0, len(df)])
ay.set_ylim([0, 15])
ay.set_xlabel('Chromosome')
ay.set_ylabel('-log10(p_value)')
ay.set_title('CB1908_IC50_gwas_results')

plt.tight_layout()
fig.savefig('Manhattan_plots.png', dpi=1000)
plt.close()


df = pd.read_table('CB1908_IC50_gwas_results.assoc.linear', sep = '\s+')

#Find the smallest p-value and save the SNP id
top_SNP = df[df['P']==df['P'].min()]
top_ID = top_SNP.iloc[0]['SNP']
#print(top_SNP)
#9 and beyond is where the genotypes are in the vcf

genotypes = []

for line in open('gwas_data/genotypes.vcf'):
    if line.startswith('#'):
        continue
    if line.find(top_ID) == -1:
        continue
    target_line = line.rstrip().split('\t')
    genotypes = target_line[9:]

phenotypes = pd.read_table('gwas_data/CB1908_IC50.txt', sep = '\t')

phenotype_vals = phenotypes['CB1908_IC50'].tolist()

het_IC50 = []
homo_IC50 = []
wt_IC50 = []

#Determine if each sample is mutant homozygous, mutant heterozygous, or wt homozygous.
#Then place IC50s into seperate lists
#additionally filter out missing IC50 samples 
for i in range(len(genotypes)):
    if genotypes[i] == '0/0':
        genotypes[i] = 'AA'
        if np.isnan(phenotype_vals[i]) == False:
            wt_IC50.append(phenotype_vals[i])
        else:
            continue
    elif genotypes[i] in ['0/1', '1/0']:
        genotypes[i] = 'AG'
        if np.isnan(phenotype_vals[i]) == False:
            het_IC50.append(phenotype_vals[i])
        else:
            continue
    elif genotypes[i] == '1/1':
        genotypes[i] = 'GG'
        if np.isnan(phenotype_vals[i]) == False:
            homo_IC50.append(phenotype_vals[i])
        else:
            continue
    else:
        genotypes[i] = 'No Data'


fig, ax = plt.subplots()

ax.boxplot([wt_IC50, het_IC50, homo_IC50])
ax.set_ylabel("IC50")
ax.set_xlabel(f'{top_ID} Genotype')
ax.set_title(f'IC50 by {top_ID} Genotype')
ax.set_xticklabels(['AA', 'AG', 'GG'])
fig.savefig('Boxplot.png', dpi = 1000)




'''
df = pd.read_table('GS451_IC50_gwas_results.assoc.linear', sep = '\s+')

#Find the smallest p-value and save the SNP id
top_SNP = df[df['P']==df['P'].min()]
top_ID = top_SNP.iloc[0]['SNP']
print(top_SNP)
'''
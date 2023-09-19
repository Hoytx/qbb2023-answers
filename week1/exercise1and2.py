#!/usr/bin/env python

import pandas as pd
import numpy as np

df = pd.read_csv('aau1043_dnm.csv')

converted = df[['Proband_id', "Phase_combined"]].value_counts()

Valuecount_dataframe = df[['Proband_id', "Phase_combined"]].value_counts().reset_index().rename(columns = {"index": "Proband_id", 0: "Phase_combined", 1: "Count"})
#above code converts value count output to a dataframe
#need to create a dictionary with key values equal to all of the proband ids and the value for each being equal to [0,0]
#then use for loop to grab the actual data and insert the results to the dictionary
ids = list(df["Proband_id"])

deNovoCount = { }


Valuecount_dataframe.sort_values(by=['Phase_combined', 'Proband_id'], ascending=True, inplace=True) 
#sort dataframe so all mother and father data points are grouped with father above mother

for i in range(len(ids)):
    deNovoCount[ids[i]] = [0,0]
loopdex = 0
# maternal_key = 0
# maternal_count = 0
# paternal_key = 0
# paternal_count = 0

for i in range(Valuecount_dataframe.shape[0]):
    #grab maternal or paternal data
    if Valuecount_dataframe.iat[loopdex, 1] == 'mother':
        maternal_count = 0
        maternal_key = 0
        maternal_count = Valuecount_dataframe.iat[loopdex, 2]
        maternal_key = Valuecount_dataframe.iat[loopdex, 0]
        deNovoCount[maternal_key] = [maternal_count, deNovoCount.pop(maternal_key)]
        #because data for mother will be read after father data I need to keep the old data
        #I used pop to grab it and stick it in with the new data because values associated with keys are only overwritten
    elif Valuecount_dataframe.iat[loopdex, 1] == 'father':    
        paternal_count = 0
        paternal_key = 0
        paternal_count = Valuecount_dataframe.iat[loopdex, 2]
        paternal_key = Valuecount_dataframe.iat[loopdex, 0]
        deNovoCount[paternal_key] = paternal_count
    #update proband_id key with values
    loopdex = loopdex + 1


deNovoCountDF = pd.DataFrame.from_dict(deNovoCount, orient = 'index', columns = ['maternal_dnm', 'paternal_dnm'])
#converts my dictionary into a dataframe 

df2 = pd.read_csv('aau1043_parental_age.csv', index_col = 'Proband_id') #needed to add index_col so it would treat the proband id as the index numbers

combinedDF = pd.concat([df2, deNovoCountDF], axis = 1, join = 'outer')

import statsmodels.formula.api as smf
import matplotlib.pyplot as plt

maternal_age_list =[]
paternal_age_list = []
maternal_mut = []
paternal_mut = []

for i in range(combinedDF.shape[0]):
    paternal_age_list.append(combinedDF.iat[i, 0])
    maternal_age_list.append(combinedDF.iat[i, 1])
    maternal_mut.append(combinedDF.iat[i, 2])
    paternal_mut.append(combinedDF.iat[i, 3])
fig, ax = plt.subplots()

ax.scatter(maternal_age_list, maternal_mut)
ax.set_xlabel("Maternal age")
ax.set_ylabel("Maternal de novo mutations")
plt.savefig('ex2_a.png')
plt.close()

fig, ax = plt.subplots()

ax.scatter(paternal_age_list, paternal_mut)
ax.set_xlabel("Paternal age")
ax.set_ylabel("Paternal de novo mutations")
plt.savefig('ex2_b.png')
plt.close()

model = smf.ols(formula = "maternal_dnm ~ 1 + Mother_age", data = combinedDF)
results = model.fit()
print(results.summary())

model2 = smf.ols(formula = "paternal_dnm ~ 1 + Father_age", data = combinedDF)
results2 = model2.fit()
print(results2.summary())

fig, ax = plt.subplots()

ax.hist(paternal_mut, alpha = 0.5)
ax.hist(maternal_mut, alpha = 0.5)
ax.set_xlabel('De novo mutations')
ax.set_ylabel('# of probands')
ax.legend(['Parental', 'Maternal'])
plt.savefig('ex2_c.png')
plt.close()

import scipy.stats as sts

stat_test = sts.ttest_ind(maternal_mut, paternal_mut)
print(stat_test)

#size is slope (coeffiecnt)
#probability of f-statistic tells you if it is significant

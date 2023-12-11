#!/usr/bin/env python

import numpy as np
import pandas as pd
from pydeseq2 import preprocessing
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#using 2018 week 35 data https://github.com/rfordatascience/tidytuesday/tree/master/data/2018/2018-11-27

#read in data

bmore_bridges = pd.read_csv("baltimore_bridges.csv")

print(set(bmore_bridges["owner"]))

#bar plot of disrepair status stratified by ownership
park_owned = bmore_bridges.loc[bmore_bridges["owner"] == "National Park Service"] 
state_owned = bmore_bridges.loc[bmore_bridges["owner"] == "State Highway Agency"]
other_state = bmore_bridges.loc[bmore_bridges["owner"] == "Other State Agencies"]
town_owned = bmore_bridges.loc[bmore_bridges["owner"] == "Town or Township Highway Agency"]  
state_toll = bmore_bridges.loc[bmore_bridges["owner"] == "State Toll Authority"] 
navy_marine_owned = bmore_bridges.loc[bmore_bridges["owner"] == "Navy/Marines"]
army_owned = bmore_bridges.loc[bmore_bridges["owner"] == "Army"]
county_owned = bmore_bridges.loc[bmore_bridges["owner"] == "County Highway Agency"]
city_owned = bmore_bridges.loc[bmore_bridges["owner"] == "City or Municipal Highway Agency"] 
private_owned = bmore_bridges.loc[bmore_bridges["owner"] == "Private (other than railroad)"] 
other_local = bmore_bridges.loc[bmore_bridges["owner"] == "Other Local Agencies"] 

owners = [county_owned, state_toll, city_owned, state_owned]
names = ['County Highway Agency', 'State Toll Authority', 'City or Municipal Highway Agency', 'State Highway Agency']



def find_counts(owner):
    good_cond = owner.loc[owner["bridge_condition"] == "Good"].shape[0]
    fair_cond = owner.loc[owner["bridge_condition"] == "Fair"].shape[0]
    poor_cond = owner.loc[owner["bridge_condition"] == "Poor"].shape[0]
    return [good_cond, fair_cond, poor_cond]

counts = {}

good = []
fair = []
poor = []
for i in range(len(owners)):
    good.append(find_counts(owners[i])[0])
    fair.append(find_counts(owners[i])[1])
    poor.append(find_counts(owners[i])[2])

counts["Good"] = good
counts["Fair"] = fair
counts["Poor"] = poor


print(counts.items())

x = np.arange(len(names))  # the label locations
width = 0.25  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for conditions, count in counts.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, count, width, label=conditions)
    multiplier += 1

ax.set_xticks(x + width, names, size = 'x-small', rotation = 45, horizontalalignment = 'right')
ax.legend(loc='upper left', ncols=3)
fig.savefig("Bridge_Condition_By_Ownership.png")


#box plot of average traffic for each disrepair status
good = bmore_bridges.loc[bmore_bridges["bridge_condition"] == "Good"]
fair = bmore_bridges.loc[bmore_bridges["bridge_condition"] == "Fair"]
poor = bmore_bridges.loc[bmore_bridges["bridge_condition"] == "Poor"]

fig2, ay = plt.subplots()

ay.boxplot([np.log(good["avg_daily_traffic"]), np.log(fair["avg_daily_traffic"]), np.log(poor["avg_daily_traffic"])], labels = ["Good Condition", "Fair Condition", "Poor Condition"])
ay.set_ylabel("log(Average Daily Traffic)")

fig2.savefig("Condition_vs_Traffic_Volume.png")

#scatterplot of bridges condition within city limits based on logitude + latitude
city_bridges = bmore_bridges.loc[bmore_bridges["county"] == "Baltimore City"]

city_good = good.loc[good["county"] == "Baltimore City"]
city_fair = fair.loc[fair["county"] == "Baltimore City"]
city_poor = poor.loc[poor["county"] == "Baltimore City"]

fig4, az = plt.subplots()

az.scatter(city_poor["long"], city_poor["lat"], s=4, color = "red", alpha = 0.6)
az.scatter(city_fair["long"], city_fair["lat"], s=4, color = "orange", alpha = 0.6)
az.scatter(city_good["long"], city_good["lat"], s=4, color = "blue", alpha = 0.6)
az.legend(["Poor", "Fair", "Good"], loc = "lower left")
az.set_xlabel("Longitude")
az.set_ylabel("Latitude")

fig4.savefig("Baltimore_City_Bridges.png")

plt.show()
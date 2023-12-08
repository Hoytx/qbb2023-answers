#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd

#loading needed files
baitmap_file = sys.argv[1]
washu_file = sys.argv[2]
output_file = sys.argv[3]

#read baitmap and output washu file as dataframes
baitmap = pd.read_table(baitmap_file, names = ["chrom", "start", "end", "id", "gene"])
washu = pd.read_table(washu_file, names = ["frag1", "frag2", "strength"])
#seperate the fragment data into seperate columns
washu[["frag1chrom", "frag1Start", "frag1End"]] = washu["frag1"].str.split(',', expand=True)
washu[["frag2chrom", "frag2Start", "frag2End"]] = washu["frag2"].str.split(',', expand=True)

#find max strength
max_strength = washu["strength"].max()

#calculate the score for each interaction
strengths_list = washu["strength"].to_list()
score_list = []
for i in range(len(strengths_list)):
    score_list.append(int((strengths_list[i]/max_strength) * 1000))
washu["score"] = score_list


#reorganize the interactions based on chromosome position
for i in range(washu.shape[0]):
    if int(washu.iloc[i]["frag1Start"]) > int(washu.iloc[i]["frag2End"]):
        #temporary storage of values
        frag1s = washu.iloc[i]["frag1Start"]
        frag1e = washu.iloc[i]["frag1End"]
        frag2s = washu.iloc[i]["frag2Start"]
        frag2e = washu.iloc[i]["frag2End"]
        #reordering values
        washu.at[i, "frag1Start"] = frag2s
        washu.at[i, "frag1End"] = frag2e
        washu.at[i, "frag2Start"] = frag1s
        washu.at[i, "frag2End"] = frag1e
    else:
        continue


#for each interaction I need to find the bait fragment in the baitmap and add the info to the washu dataframe
#create new columns to identify the source and the target
#if both are bait lower fragment will be treated as source, higher as target
#set targetStrand to '+' if both fragments were bait and '-' if it was not


for i in range(washu.shape[0]):
    if (int(washu.iloc[i]["frag1Start"]) in baitmap["start"].tolist()) and (int(washu.iloc[i]["frag2Start"]) in baitmap["start"].tolist()):
        washu.at[i, "sourceChrom"] = washu.iloc[i]["frag1chrom"]
        washu.at[i, "sourceStart"] = washu.iloc[i]["frag1Start"]
        washu.at[i, "sourceEnd"] = washu.iloc[i]["frag1End"]
        washu.at[i, "targetChrom"] = washu.iloc[i]["frag2chrom"]
        washu.at[i, "targetStart"] = washu.iloc[i]["frag2Start"]
        washu.at[i, "targetEnd"] = washu.iloc[i]["frag2End"]
        washu.at[i, "targetStrand"] = "+"
    elif (int(washu.iloc[i]["frag1Start"]) in baitmap["start"].tolist()) and (int(washu.iloc[i]["frag2Start"]) not in baitmap["start"].tolist()):
        washu.at[i, "sourceChrom"] = washu.iloc[i]["frag1chrom"]
        washu.at[i, "sourceStart"] = washu.iloc[i]["frag1Start"]
        washu.at[i, "sourceEnd"] = washu.iloc[i]["frag1End"]
        washu.at[i, "targetChrom"] = washu.iloc[i]["frag2chrom"]
        washu.at[i, "targetStart"] = washu.iloc[i]["frag2Start"]
        washu.at[i, "targetEnd"] = washu.iloc[i]["frag2End"]
        washu.at[i, "targetStrand"] = "-"
    elif (int(washu.iloc[i]["frag1Start"]) not in baitmap["start"].tolist()) and (int(washu.iloc[i]["frag2Start"]) in baitmap["start"].tolist()):
        washu.at[i, "sourceChrom"] = washu.iloc[i]["frag2chrom"]
        washu.at[i, "sourceStart"] = washu.iloc[i]["frag2Start"]
        washu.at[i, "sourceEnd"] = washu.iloc[i]["frag2End"]
        washu.at[i, "targetChrom"] = washu.iloc[i]["frag1chrom"]
        washu.at[i, "targetStart"] = washu.iloc[i]["frag1Start"]
        washu.at[i, "targetEnd"] = washu.iloc[i]["frag1End"]
        washu.at[i, "targetStrand"] = "-"
    else:
        print("error")
        break

print(washu)
#create empty dataframe with format of final bed file
UCSC_bed_file = pd.DataFrame(columns = ["chrom", "chromStart", "chromEnd", "name", "score", "value", "exp", "color", "sourceChrom", "sourceStart", "sourceEnd", "sourceName", "sourceStrand", "targetChrom", "targetStart", "targetEnd", "targetName", "targetStrand"])

#add all of the interactions to the bed file
UCSC_bed_file["chrom"],UCSC_bed_file["chromStart"],UCSC_bed_file["chromEnd"] = washu["frag1chrom"],washu["frag1Start"],washu["frag2End"]
UCSC_bed_file["score"],UCSC_bed_file["value"] = washu["score"],washu["strength"]
UCSC_bed_file["sourceChrom"], UCSC_bed_file["sourceStart"], UCSC_bed_file["sourceEnd"] = washu["sourceChrom"], washu["sourceStart"], washu["sourceEnd"]
UCSC_bed_file["targetChrom"], UCSC_bed_file["targetStart"], UCSC_bed_file["targetEnd"] = washu["targetChrom"], washu["targetStart"], washu["targetEnd"]
UCSC_bed_file["targetStrand"] = washu["targetStrand"]

#set name, exp, color, and sourceStrand data for all rows since they are the same for all interactions
UCSC_bed_file["name"] = "."
UCSC_bed_file["exp"] = "."
UCSC_bed_file["color"] = "0"
UCSC_bed_file["sourceStrand"] = "+"

#for each interaction find the name of the bait fragment genes and add them to the dataframe, if only one fragment is a bait add "." to the target name column
for i in range(UCSC_bed_file.shape[0]):
    if UCSC_bed_file.iloc[i]["targetStrand"] == "+":
        UCSC_bed_file.at[i, "sourceName"] = baitmap.loc[baitmap["start"] == int(UCSC_bed_file.iloc[i]["sourceStart"]), "gene"].iloc[0]
        UCSC_bed_file.at[i, "targetName"] = baitmap.loc[baitmap["start"] == int(UCSC_bed_file.iloc[i]["targetStart"]), "gene"].iloc[0]
    elif UCSC_bed_file.iloc[i]["targetStrand"] == "-":
        UCSC_bed_file.at[i, "sourceName"] = baitmap.loc[baitmap["start"] == int(UCSC_bed_file.iloc[i]["sourceStart"]), "gene"].iloc[0]
        UCSC_bed_file.at[i, "targetName"] = "."
    else:
        print("error")
        break



output = open(output_file, 'w+')
output.write('track type=interact name="pCHIC" description="Chromatin interactions" useScore=on maxHeightPixels=200:100:50 visibility=full \n')
output.write(UCSC_bed_file.to_string(index=False, header=False))



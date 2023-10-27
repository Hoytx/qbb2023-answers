#!/usr/bin/env python

import matplotlib.pyplot as plt

read_depths = []
genotype_qualities = []
allele_frequencies = []
#variant_type_counts = {"snp": 0, "mnp": 0, "ins": 0, "del": 0, "complex": 0}
variant_effects = []
variant_count = []
for line in open('annotated_variants.vcf'):
    if line.startswith('#'):
        continue
    fields = line.rstrip('\n').split('\t')
    
    samples = fields[9:]
    info = fields[7]
    for i in range(len(samples)):
        reads_split = samples[i].split(':')
        if reads_split[2] != ".":
            read_depths.append(int(reads_split[2]))
        elif reads_split[2] == ".":
            read_depths.append(0)
        else:
            print("edge case")
    
    for i in range(len(samples)):
        reads_split = samples[i].split(':')
        if reads_split[1] != ".":
            genotype_qualities.append(float(reads_split[1]))
        elif reads_split[1] == ".":
            genotype_qualities.append(0)
        else:
            print("edge case")
    split_info = info.split(";")
    
    #seperate out allele frequency data to deal with multiple alleles at the same site
    allele_frequency = split_info[3][3:]
    split_frequencies = allele_frequency.split(',')
    allele_frequencies.extend(split_frequencies)
    
    for i in range(len(split_info)):
        if "ANN=" in split_info[i]:
            annotations = split_info[i] 
    annotation_split = annotations.split("|")
    for i in annotation_split:
             x = annotation_split[1]
             variant_effects.append(x)

for i in set(variant_effects):
    variant_count.append(variant_effects.count(i))

types_of_effects = list(set(variant_effects))

#convert frequency to percent scale for graphing
allele_percents = []
for i in range(len(allele_frequencies)):
    x = float(allele_frequencies[i]) * 100
    allele_percents.append(x)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

ax1.hist(read_depths, bins = range(min(read_depths), max(read_depths) + 1, 1))
ax1.set_xlim(0,50)
ax1.set_xlabel("Coverage Depth")
ax1.set_ylabel("# of Samples")

ax2.hist(genotype_qualities, bins = range(0, 200, 5))
ax2.set_xlabel("Genotype Quality (Phred-scaled)")
ax2.set_ylabel("# of Samples")

ax3.hist(allele_percents)
ax3.set_xlabel("Allelic Frequency %")
ax3.set_ylabel("# of Variants")

ax4.bar(types_of_effects, variant_count)
ax4.set_xticklabels(types_of_effects, rotation = 45, ha= "right", fontsize = 5.5)
ax4.set_ylabel("# of Variants")

plt.tight_layout()
plt.savefig("Exploratory_Analysis.png", dpi=1200)
plt.show()




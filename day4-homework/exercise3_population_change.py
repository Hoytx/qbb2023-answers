#!/usr/bin/env python

import numpy as np 
import matplotlib.pyplot as plt


def simulate_wf(population_size, frequency):

    tot_alleles = 2 * population_size
    allele_frequencies = [frequency]
    trial = 0

    while 0 < frequency < 1:
        trial  = np.random.binomial(tot_alleles, frequency)

        frequency = trial / tot_alleles

        allele_frequencies.append(frequency)

    generations = len(allele_frequencies)    


    return allele_frequencies, generations


fig, ax = plt.subplots()

r_1 = 0
r_2 = 0
r_3 = 0
r_4 = 0
r_5 = 0
r_6 = 0
r_7 = 0
r_8 = 0
r_9 = 0
r_10 = 0

for i in range(50):
    output = simulate_wf(100, 0.5)
    r_1 = r_1 + output[1]
for i in range(50):
    output = simulate_wf(200, 0.5)
    r_2 = r_2 + output[1]
for i in range(50):
    output = simulate_wf(250, 0.5)
    r_3 = r_3 + output[1]
for i in range(50):
    output = simulate_wf(500, 0.5)
    r_4 = r_4 + output[1]
for i in range(50):
    output = simulate_wf(1000, 0.5)
    r_5 = r_5 + output[1]


av_1 = r_1 / 50
av_2 = r_2 / 50
av_3 = r_3 / 50
av_4 = r_4 / 50
av_5 = r_5 / 50    


change_pop = [av_1, av_2, av_3, av_4, av_5]
population_trial = [100, 200, 250, 500, 1000]
ax.plot(population_trial, change_pop)
plt.xlabel("population")
plt.ylabel("generations to fixation")


#    for i in range(50):
#        output = simulate_wf(250, 0.5)
#        ax.plot(range(len(output[0])), output[0])
#dont need now but don't want to delete
plt.savefig( "Population_Change_Fixation.jpg" )
plt.show()



#print(simulate_wf(250, 0.5))

# Get starting frequency and a population size
# Get imput parameters for function

# Make list to store our allele frequencies

# While our allele frequency is between 0 and 1:
#   Get the new allele frequency for next generation
#   by drawing from the binomial distribution
#   (convert number of successes into a frequency)
#
#   Store our allele frequency in the AF list


# Return list of allele frequency at each time 
# Number of gererations to fixation
# is the length of your list

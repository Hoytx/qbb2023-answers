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

    fig, ax = plt.subplots()
    ax.plot(range(len(allele_frequencies)), allele_frequencies)


    return allele_frequencies, generations




print(simulate_wf(250, 0.5))
plt.xlabel("Generations")
plt.ylabel("Allele Frequency")
plt.savefig( "exercise_1_graph.jpg" )
plt.show()

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

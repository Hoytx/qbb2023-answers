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

#adjusted to just give generations to fixation
    return generations


fig, ax = plt.subplots()

average_result = [] #created empty list to hold results of avarage time to fixation for each allele frequency

trial_frequencies = [0.2, 0.5, 0.6, 0.7, 0.9] #allele frequencies to try

#make a loop that runs 5 times
#in that loop make a loop that runs 20 times
#the inner loop takes the trial frequency at the index number equal to how many times we have run the outer loop (starting with zero)
#and adds the resulting time to fixation to the output value
#the outer loop should select the next trial frequency, run the inner loop, convert the output into the average time to fixation
#and append the result to my list of results
#outer loop should also reset the output value to zero right before restarting the loop


output = 0

for i in range(5):
    frequency_index = i
    for j in range(40):
        output = output + simulate_wf(1000, trial_frequencies[frequency_index])
    output = output / 40
    average_result.append(output)
    print("starting frequency: " + str(trial_frequencies[frequency_index]) + " average generations to fixation: "+ str(output))
    #this print statement shows what the starting allele frequency is for each time through the loop and the average time to fixation
    #this mainly served as a sanity check for the graphed results
    output = 0 


ax.plot(trial_frequencies, average_result)
plt.xlabel("allele frequency")
plt.ylabel("average generations to fixation")


#    for i in range(50):
#        output = simulate_wf(250, 0.5)
#        ax.plot(range(len(output[0])), output[0])
#dont need now but don't want to delete
plt.savefig( "Allele_Frequency_Change_Fixation.jpg" )
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

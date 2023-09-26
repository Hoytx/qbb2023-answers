#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

#Exercise 1

def simulate_coverage(coverage, genome_len, read_len, figname):

    coverage_arr = np.zeros(genome_len)
    
    num_reads = int(coverage * genome_len / read_len)

    low = 0
    high = genome_len - read_len 

    start_positions = np.random.randint(low=0, high= high + 1, size = num_reads) # high is exclusive

    for start in start_positions:
        coverage_arr[start: start+read_len] += 1

    x = np.arange(0, max(coverage_arr)+1)

    sim_0cov = genome_len - np.count_nonzero(coverage_arr)
    sim_0cov_pct = 100* sim_0cov / genome_len

    print(f'In the simulation, there are {sim_0cov} bases with 0 coverage')
    print(f'This is {sim_0cov_pct}% of the genome')

    # Get poisson distribution
    y_poisson = stats.poisson.pmf(x ,mu=coverage) * genome_len
   
    #Get normal_distribution
    y_normal = stats.norm.pdf(x, loc = coverage, scale = np.sqrt(coverage)) * genome_len
    

    fig, ax = plt.subplots()
    ax.hist(coverage_arr, bins=x, label = 'Simulation', alpha = 0.75)
    ax.plot(x, y_poisson, label = 'Poisson', alpha = 0.75)
    ax.plot(x, y_normal, label = 'Normal', alpha = 0.75)
    ax.set_xlabel('Coverage')
    ax.set_ylabel('Frequency (bp)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(figname)


simulate_coverage(3, 1_000_000, 100, 'ex1_3x_cov.png')

simulate_coverage(10, 1_000_000, 100, 'ex1_10x_cov.png')

simulate_coverage(30, 1_000_000, 100, 'ex1_30x_cov.png')

#----------
# Exercise 2
#----------

reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']

k = 3
kmers = []

for i in range(len(reads)):
    read = reads[i]
    for j in range(len(read) - k):
        kmer1 = read[j:j+k]
        kmer2 = read[j+1:j+k+1]
        kmers.append(kmer1 + " -> " + kmer2)


edges = list(set(kmers))

formatted_edges = ['digraph sample {\n']

file_for_dots = open("dots.dot", "w")

for i in range(len(edges)):
    item = edges[i] + ';' + '\n'
    formatted_edges.append(item)       
formatted_edges.append('}')    
file_for_dots.writelines(formatted_edges)
file_for_dots.close()

print(edges)
print(kmers)


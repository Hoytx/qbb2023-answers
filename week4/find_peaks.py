#!/usr/bin/env python

import sys

from model_peaks import load_bedgraph, bin_array
import numpy
import scipy.stats
import matplotlib.pyplot as plt


def main():
    # Load file names and fragment width
    forward = sys.argv[1]
    reverse = sys.argv[2]
    controlfwd = sys.argv[3]
    controlrev = sys.argv[4]  
    fragment_width = int(sys.argv[5])
    output_wiggle = sys.argv[6]
    output_bed = sys.argv[7]
    # Define what genomic region we want to analyze
    chrom = "chr2R"
    chromstart = 10000000
    chromend =  12000000
    chromlen = chromend - chromstart
    offset = fragment_width 

    # Load the sample bedgraph data, reusing the function we already wrote
    forward_data = load_bedgraph(forward, chrom, chromstart, chromend)
    reverse_data = load_bedgraph(reverse, chrom, chromstart, chromend)
    # Combine tag densities, shifting by our previously found fragment width
    data_combined = forward_data[offset:] + reverse_data[:-offset]
    # Load the control bedgraph data, reusing the function we already wrote
    control_forward = load_bedgraph(controlfwd, chrom, chromstart, chromend)    
    control_reverse = load_bedgraph(controlrev, chrom, chromstart, chromend) 
    # Combine tag densities
    control_combined = control_forward + control_reverse
    # Adjust the control to have the same coverage as our sample
    print(data_combined.shape)
    #adj_control = control_combined[offset//2:len(control_combined)]
    adj_control = control_combined
    adj_control = adj_control / numpy.sum(adj_control) * numpy.sum(data_combined)

    # Create a background mean using our previous binning function and a 1K window
    # Make sure to adjust to be the mean expected per base
    background_mean = bin_array(adj_control, 1000)/1000
    # Find the mean tags/bp and make each background position the higher of the
    # the binned score and global background score
    global_mean = numpy.mean(adj_control)

    local_background = numpy.zeros(data_combined.shape[0], dtype = float)

    for i in range(len(local_background)):
        if background_mean[i] > global_mean:
            local_background[i] = background_mean[i]
        else:
            local_background[i] = global_mean

    # Score the sample using a binsize that is twice our fragment size
    # We can reuse the binning function we already wrote
    sample_score = bin_array(data_combined, fragment_width * 2)

    # Find the p-value for each position (you can pass a whole array of values
    # and and array of means). Use scipy.stats.poisson for the distribution.
    # Remeber that we're looking for the probability of seeing a value this large
    # or larger
    # Also, don't forget that your background is per base, while your sample is
    # per 2 * width bases. You'll need to adjust your background
    
    #cdf function of poisson 

    values = scipy.stats.poisson.cdf(sample_score, (local_background * 2 * fragment_width))


    # Transform the p-values into -log10
    # You will also need to set a minimum pvalue so you doen't get a divide by
    # zero error. I suggest using 1e-250
    pvalues = 1e-250
    pvalues = -numpy.log10(values)
    # Write p-values to a wiggle file
    # The file should start with the line
    # "fixedStep chrom=CHROM start=CHROMSTART step=1 span=1" where CHROM and
    # CHROMSTART are filled in from your target genomic region. Then you have
    # one value per line (in this case, representing a value for each basepair).
    # Note that wiggle files start coordinates at 1, not zero, so add 1 to your
    # chromstart. Also, the file should end in the suffix ".wig"
    write_wiggle(pvalues, chrom, chromstart + 1, output_wiggle)
    # Write bed file with non-overlapping peaks defined by high-scoring regions 
    write_bed(pvalues, chrom, chromstart, chromend, fragment_width, output_bed)








def write_wiggle(pvalues, chrom, chromstart, fname):
    output = open(fname, 'w')
    print(f"fixedStep chrom={chrom} start={chromstart + 1} step=1 span=1",
          file=output)
    for i in pvalues:
        print(i, file=output)
    output.close()

def write_bed(scores, chrom, chromstart, chromend, width, fname):
    chromlen = chromend - chromstart
    output = open(fname, 'w')
    while numpy.amax(scores) >= 10:
        pos = numpy.argmax(scores)
        start = pos
        while start > 0 and scores[start - 1] >= 10:
            start -= 1
        end = pos
        while end < chromlen - 1 and scores[end + 1] >= 10:
            end += 1
        end = min(chromlen, end + width - 1)
        print(f"{chrom}\t{start + chromstart}\t{end + chromstart}", file=output)
        scores[start:end] = 0
    output.close()


if __name__ == "__main__":
    main()
#!/usr/bin/env python

import numpy as np
from fasta import readFASTA
import pandas as pd
import sys


fasta_file = sys.argv[1]
score_matrix_file = sys.argv[2]
gap_penalty = float(sys.argv[3])
output_file = sys.argv[4]

input_sequences = readFASTA(open(fasta_file))

seq1_id, sequence1 = input_sequences[0]
seq2_id, sequence2 = input_sequences[1]

score_matrix = pd.read_csv(score_matrix_file, sep = '\s+')

#Make 2 empty matrices, one for optimal local alightments, one for the optimal global path

F_matrix = np.zeros((len(sequence1)+1, len(sequence2)+1))

traceback_matrix = np.empty((len(sequence1)+1, len(sequence2)+1), dtype = str)

#Fill in first column
for i in range(len(sequence1)+1):
    F_matrix[i,0] = i*gap_penalty
    traceback_matrix[i,0] = "U"
#Fill in first row
for j in range(len(sequence2)+1):
    F_matrix[0,j] = j*gap_penalty
    traceback_matrix[0,j] = "L"    
#reset top left of traceback matrix to indicate the end of the sequence
traceback_matrix[0,0] = "E"


for i in range(1, len(sequence1)+1): # loop through rows
    for j in range(1, len(sequence2)+1): # loop through columns
        x = sequence2[j-1]
        y = sequence1[i-1]
        if sequence1[i-1] == sequence2[j-1]: # if sequence1 and sequence2 match at positions i and j, respectively...
            d = F_matrix[i-1, j-1] + float(score_matrix.loc[x, y])            
        else: # if sequence1 and sequence2 don't match at those positions...the score matrix will give the mismatch penalty
            d = F_matrix[i-1, j-1] + float(score_matrix.loc[x, y])
        h = F_matrix[i,j-1] + gap_penalty
        v = F_matrix[i-1,j] + gap_penalty
        if max(d,h,v) == d:
            F_matrix[i,j] = d
            traceback_matrix[i,j] = 'D'
        elif max(d,h,v) == h:
            F_matrix[i,j] = h
            traceback_matrix[i,j] = 'L'
        elif max(d,h,v) == v:
            F_matrix[i,j] = v
            traceback_matrix[i,j] = 'U'
#Setting myself up to start in the bottom right corner
i = len(sequence1)
j = len(sequence2)
alignment_score = F_matrix[i,j]
seq1_align = ""
seq2_align = ""
path = ''
gaps1 = 0
gaps2 = 0
while i != 0 or j != 0:
    if traceback_matrix[i,j] == 'D':
        seq1_align = sequence1[i-1] + seq1_align
        seq2_align = sequence2[j-1] + seq2_align
        i -= 1
        j -= 1
        path = path + 'D'
    elif traceback_matrix[i,j] == 'U':
        seq1_align = sequence1[i-1] + seq1_align
        seq2_align = '-' + seq2_align
        i -= 1
        path = path + 'U'
        gaps2 += 1
    elif traceback_matrix[i,j] == 'L':
        seq1_align = '-' + seq1_align
        seq2_align = sequence2[j-1] + seq2_align
        j -= 1
        path = path + 'L'
        gaps1 += 1       

print(F_matrix)
print(traceback_matrix)
print("The alignment score was: " + str(alignment_score))
#print(seq1_align)
#print(seq2_align)
print(f"Sequence one had", gaps1, "gaps.")
print(f"Sequence two had", gaps2, "gaps.")

alignment_output = open(output_file, "w")
alignment_output.write("Sequence 1:\n " + seq1_align + "\nSequence 2:\n " + seq2_align)
alignment_output.close()
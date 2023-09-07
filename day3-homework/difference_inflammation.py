#!/usr/bin/env python

def inflam_mean_dif(patient1_id, patient2_id, filename = ""):
    
    index = 0
    f = open(filename, "r")
    patients = f.readlines()

    for patient in patients:
        patient = patient.rstrip() #strip off linebreak
        patient = patient.split(",") # split data points to parts of a list
        patients[index] = patient #creating list of lists with the index number being patient number
        index = index + 1  

    patient1 = []        
    patient2 = []

    for i in patients[patient1_id]:
        total = int(i)
        patient1.append(total) 

    for i2 in patients[patient2_id]:
        total = int(i2)
        patient2.append(total)

    differences = []    

    index_1 = 0

    for diff in range(40):
        differences.append(patient1[index_1] - patient2[index_1])
        index_1 = index_1 + 1
    
    return differences
    

print(inflam_mean_dif(1,2,"inflammation-01.csv"))
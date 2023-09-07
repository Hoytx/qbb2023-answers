#!/usr/bin/env python

def inflam_mean(patient_id, filename = ""):
    
    index = 0
    f = open(filename, "r")
    patients = f.readlines()

    for patient in patients:
        patient = patient.rstrip() #strip off linebreak
        patient = patient.split(",") # split data points to parts of a list
        patients[index] = patient #creating list of lists with the index number being patient number
        index = index + 1   

    total = 0
    for number in patients[patient_id]:
        total = int(total)
        total = total + int(number)
    return(total / 40)        

print(inflam_mean(0, "inflammation-01.csv"))            
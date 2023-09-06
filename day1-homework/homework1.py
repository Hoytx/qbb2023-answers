#!/usr/bin/env python
import numpy

f = open("inflammation-01.csv", "r")


patients = f.readlines()

print("Excercise 1")

index = 0

for patient in patients:
	patient = patient.rstrip() #strip off linebreak
	patient = patient.split(",") # split data points to parts of a list
	patients[index] = patient #creating list of lists with the index number being patient number
	index = index + 1

#--excercise 1----------
print(patients[4][0]) #flare ups of 5th patient on 1st day
print(patients[4][9]) #flare ups of 5th patient on 10th day
print(patients[4][-1]) #flare ups of 5th patient on last day
#----------------
#Results 0, 4, 1
#----------------

print("Excercise 2")

#---exercise 2---------
patients_int = []
for patient in patients:
	patient_int = []
	for day in patient:
		day = int(day)
		patient_int.append(day)
	patients_int.append(patient_int)

#calculate mean for first 10 patients

for flareups in patients_int[0:10]:
	totalflareups = numpy.mean(flareups)
	print(totalflareups)
#-------------
#Results 5.45, 5.425, 6.1, 5.9, 5.55, 6.225, 5.975, 6.65, 6.625, 6.525
#--------------

print("Exercise 3")


#---excercise 3---------
flareup_means = []
for flareups in patients_int[0:10]: #for first 10 patients
	totalflareups = numpy.mean(flareups) #find mean
	flareup_means.append(totalflareups) #add each mean as item in list
print(numpy.max(flareup_means)) #find max of means
print(numpy.min(flareup_means)) #find min of means
#-------------
#Results Max = 6.65 Min = 5.425
#-----------------------

print("Excercise 4")

#---excercise 4------
day_index = 0 #create new index for adding to list
differences = [] #empty list
for i in patients_int[0:40]: #do loop for all 40 days
	difference = (patients_int[0][day_index] - patients_int[4][day_index]) #find difference between patient 1 and 5 for number day
	day_index = day_index + 1 
	differences.append(difference) #append result to list
print(differences) #print list of 40 differences
print(len(differences)) #checking that I had the right number of values

#--------------
#Results [0, -1, 0, 0, -2, 1, 1, 2, 6, -1, -1, -4, 4, 0, 4, -6, -1, -3, 6, 1, -3, -1, 2, 4, -6, -2, -8, 0, 1, 1, -5, -2, 2, 5, 1, 0, 0, 3, -1, -1]
#-----------------

#print(numpy.mean(start_coords))
#print(numpy.max(start_coords))
#print(numpy.min(start_coords))

f.close()

#!/usr/bin/env python

def my_mean(a = []):
    total = 0
    for number in a:
        total = total + number
    return(total / len(a))

f = open("my_integers.txt", "r")
lines = f.readlines()
number_strip = []
for numbers in lines:
    striped = numbers.rstrip()
    striped = int(striped)
    number_strip.append(striped)

print(number_strip)

print(my_mean(number_strip))


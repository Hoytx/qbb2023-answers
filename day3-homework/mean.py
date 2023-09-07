#!/usr/bin/env python

def my_mean(a = []):
    total = 0
    for number in a:
        total = total + number
    return(total / len(a))

my_list = [1,2,1,2]
print(my_mean(my_list))   
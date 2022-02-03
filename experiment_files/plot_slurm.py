#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import sys

with open(sys.argv[1]) as f:
    reading=False
    i=0

    energies=[]
    hit_rates=[]

    for line in f:
        if reading == False:
            if line.endswith("good ones?\n"):
                reading = True
            continue

        if i == 0:
            energy, hit_rate = line.split(' ')
            energies.append(float(energy))
            hit_rates.append(float(hit_rate))

        if i == 1:
            #print("b",line)
            pass
        
        i+=1
        i=i%2
    print(len(energies)) 
    plt.plot(np.arange(len(energies)), energies)
    plt.plot(hit_rates)
    plt.show()

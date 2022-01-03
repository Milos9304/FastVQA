import csv
import numpy as np
import sys
import os
import glob
import subprocess
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from glob import glob

dims=[16,17,18,19,20,21,22,23] #[15
means = []
variances = []
stds = []
solved_mean = []

lbs=[]
ubs=[]

n_samples = 1000

for dim in dims:

    hit_rates = []
    solved=[]

    for i in range(127+1):

        filename = 'rank_'+str(dim)+"/statsfile_"+str(i)+"_"+str(dim)+".txt"

        with open(filename) as f:

            line = subprocess.check_output(['tail', '-1', filename])[:-1].decode("utf-8")
            splitted=line.split(' ')
            opt_energy=splitted[0]
            hit_rate=splitted[1]
            x_vect=splitted[2:]

            hit_rates.append(float(hit_rate))
            if float(hit_rate) > 0 and 1/float(hit_rate) < (n_samples+(dim-dims[0])*0): 
                solved.append(1)
            else:
                solved.append(0)

    means.append(np.mean(hit_rates))
    stds.append(np.std(hit_rates))
    solved_mean.append(np.mean(solved))
    lbs.append(np.min(hit_rates))
    ubs.append(np.max(hit_rates))

means=np.array(means)
stds=np.array(stds)
solved_mean=np.array(solved_mean)

plt.errorbar(dims, means, stds, ls='--', capsize=2)
#plt.fill_between(dims, lbs, ubs, alpha=.1)
plt.plot(dims, solved_mean)

plt.ylim((0,1))
plt.xlabel("Rank")
plt.ylabel("Probability")
plt.legend(["Probability of finding ground states with " + str(n_samples) + " samples", "Averaged overlap (+stdev) of final ansatz state with ground state"])
plt.show()

import csv
import numpy as np
import sys
import os
import glob
import subprocess
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from glob import glob

dims=[15,16,17,18,19,20,21,22,23]
means = []
variances = []
stds = []

for dim in dims:

    hit_rates = []
    for i in range(127+1):

        filename = 'rank_'+str(dim)+"/statsfile_"+str(i)+"_"+str(dim)+".txt"

        with open(filename) as f:

            line = subprocess.check_output(['tail', '-1', filename])[:-1].decode("utf-8")
            splitted=line.split(' ')
            opt_energy=splitted[0]
            hit_rate=splitted[1]
            x_vect=splitted[2:]

            hit_rates.append(float(hit_rate))

    means.append(np.mean(hit_rates))
    stds.append(np.std(hit_rates))

means=np.array(means)
stds=np.array(stds)

plt.errorbar(dims, means, stds, ls='--', capsize=2)
plt.fill_between(dims, means-stds, means+stds, alpha=.1)
plt.ylim((0,1))
plt.show()

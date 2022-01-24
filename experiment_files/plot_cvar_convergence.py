import csv
import numpy as np
import sys
import os
import glob
import subprocess
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from glob import glob

dims=[16]
means = []
variances = []
stds = []

cvars=[0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325]

low_bounds = []
high_bounds = []

w_size=50
def window_convergence(x):
    n=len(x)
    mu=np.mean(x)

    return (np.sum(np.power(x-mu, 2))/(n-1))**(1/2)

for dim in dims:
    ctimes = []

    for cvar in cvars:

        print(cvar)

        for i in range(127+1):

            filename = 'rank_'+str(dim)+"/cvar"+str(cvar)+"_statsfile_"+str(i)+"_"+str(dim)+".txt"

            hit_rates = []
            energies = []

            with open(filename) as f:

                #line = subprocess.check_output(['tail', '-1', filename])[:-1].decode("utf-8")
                hit_rates = []
                energies = []
                for line in f:

                    splitted=line.split(' ')

                    if len(splitted) > 2:
                        break

                    energy = float(splitted[0])
                    hit_rate= float(splitted[1])
                    
                    energies.append(energy)
                    hit_rates.append(hit_rate)
                    
                energies = np.array(energies)
                hit_rates = np.array(hit_rates)
                
                rmsds1=[]
                rmsds2=[]
                for i in range(w_size, len(energies)):
                    window1 = energies[i-w_size:i]
                    window2 = hit_rates[i-w_size:i]
                    rmsds1.append(window_convergence(window1))
                    rmsds2.append(window_convergence(window2))

                r=0.05
                t1=np.max(rmsds1)*r
                t2=np.max(rmsds2)*r
                for i in range(len(rmsds1)):
                    if rmsds1[i] < t1 and rmsds2[i] < t2:
                        ctimes.append(i)
                        break
        low_bounds.append(np.min(ctimes))
        high_bounds.append(np.max(ctimes))

        means.append(np.mean(ctimes))
        stds.append(np.std(ctimes))

means=np.array(means)
stds=np.array(stds)

plt.errorbar(cvars, means, stds, ls='--', capsize=2)
#plt.fill_between(dims, means-stds, means+stds, alpha=.1)
plt.fill_between(cvars, low_bounds, high_bounds, alpha=.1)
#plt.ylim((0,1))

print(means)
print(stds)

plt.title("Convergence time")
plt.xlabel("Rank")
plt.ylabel("Time = number of iterations until convergence")
plt.show()

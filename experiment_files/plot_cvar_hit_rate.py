import csv
import numpy as np
import sys
import os
import glob
import subprocess
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from glob import glob

dims=[15]
cvars=[0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1]

means = []
medians=[]
variances = []
stds = []
solved_mean = []
cvars_total=[]

hit_rates_total=[]

lbs=[]
ubs=[]

n_samples = 5000

for dim in dims:

    for cvar in cvars:

        hit_rates = []
        solved=[]

        for i in range(1023+1):

            filename = 'rank_'+str(dim)+"/cvar"+str(cvar)+"_statsfile_"+str(i)+"_"+str(dim)+".txt"

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
            cvars_total.append(cvar)
        means.append(np.mean(hit_rates))
        medians.append(np.median(hit_rates))
        stds.append(np.std(hit_rates))
        solved_mean.append(np.mean(solved))
        lbs.append(np.min(hit_rates))
        ubs.append(np.max(hit_rates))
        hit_rates_total+=hit_rates
        
means=np.array(means)
stds=np.array(stds)
solved_mean=np.array(solved_mean)
hit_rates_total=np.array(hit_rates_total)

#plt.errorbar(cvars, means, stds, ls='--', capsize=2)
#plt.fill_between(dims, lbs, ubs, alpha=.1)
#plt.plot(cvars, solved_mean)
#plt.plot(cvars, medians)
#print(len(cvars_total))
#print(len(hit_rates_total))
#plt.scatter(cvars_total+np.random.normal(0, 0.01, len(cvars_total)), hit_rates_total, alpha=0.03)

print("x y")
for i in range(len(means)):
    print(cvars[i], ' ', solved_mean[i])#medians[i])#means[i], ' ', #stds[i])

#plt.ylim((0,1))
plt.xlabel("CVaR")
plt.ylabel("Probability")
plt.legend(["Probability of finding ground states with " + str(n_samples) + " samples", "Averaged overlap (+stdev) of final ansatz state with ground state"])
plt.show()

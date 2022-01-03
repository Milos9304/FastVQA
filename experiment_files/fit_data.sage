import subprocess
import numpy as np
import os
import sys
import csv
from glob import glob

dims=[16,17]#,18,19,20,21,22,23]

n_samples = 1000

w_size=50
def window_convergence(x):
    n=len(x)
    mu=np.mean(x)

    return (np.sum(np.power(x-mu, 2))/(n-1))**(1/2)

times = []
probabilities = []
for dim in dims:
    
    print(dim)

    solved=[]
    ctimes=[]

    for i in range(127+1):

        filename = 'rank_'+str(dim)+"/statsfile_"+str(i)+"_"+str(dim)+".txt"
        with open(filename) as f:
            line = subprocess.check_output(['tail', '-1', filename])[:-1].decode("utf-8")

            splitted=line.split(' ')
            hit_rate=splitted[1]

            #hit_rates.append(float(hit_rate))
            
            if float(hit_rate) > 0 and 1/float(hit_rate) < n_samples: 
                solved.append(1)
            else:
                solved.append(0)

            energies = []
            hit_rates = []
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

    times.append(np.mean(ctimes))
    probabilities.append(np.mean(solved))

data=[]
for i in range(len(times)):
    time = times[i]
    prob = probabilities[i]
    print(time)
    print(prob)
    data.append((dims[i], log(float(time), float(prob))))

var("x,y,n")
f=3/(4*e)*n*log(n,2)+x*n+y
f=f.function(n)
print(find_fit(data, f))

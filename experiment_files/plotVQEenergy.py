import matplotlib.pyplot as plt
import numpy as np
import sys

energies = []
opt_energies = []
hit_rates = []

detailed_indexes = []
ground_overlaps = []

w_size = 50
def window_convergence(x):
    n=len(x)
    mu=np.mean(x)

    return (np.sum(np.power(x-mu,2)/(n-1)))**(1/2)

with open(sys.argv[1], encoding='utf8') as f:
    i=0
    for line in f:
        splitted=line.split(" ")
       	
        if len(splitted) == 1:
            energy = float(line)
            energies.append(energy)
        elif len(splitted) == 2:
            energy = float(splitted[0])
            ground_overlap = float(splitted[1])
            ground_overlaps.append(ground_overlap)
            energies.append(energy)
        else:
            energy = float(splitted[0])
            opt_energy = float(splitted[1])
            hit_rate = float(splitted[2])
        	
            opt_energies.append(opt_energy)
            hit_rates.append(hit_rate)
            detailed_indexes.append(i)
        	
        i+=1

x=np.linspace(1,len(energies),len(energies))

color = 'tab:red'
fig, ax1 = plt.subplots()
ax1.set_xlabel('# iterations')
ax1.set_ylabel('energy', color=color)        
ax1.plot(x, energies, 'ro', markersize=1)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis


if len(ground_overlaps) > 0:
    color = 'tab:blue'
    ax2.set_ylabel('percent', color=color)  # we already handled the x-label with ax1
    #ax2.plot(detailed_indexes, hit_rates, 'bo', markersize=6)
    ax2.plot(x, ground_overlaps, 'bo', markersize=1)
    ax2.tick_params(axis='y', labelcolor=color)

    rmsds=[]
    rmsds2=[]

    ground_overlaps=np.array(ground_overlaps)
    energies=np.array(energies)
    for i in range(w_size, len(ground_overlaps)):
        window = ground_overlaps[i-w_size:i]
        window2 = energies[i-w_size:i]
        rmsds.append(window_convergence(window))
        rmsds2.append(window_convergence(window2))
    
    r=0.05
    t1=np.max(rmsds)*r
    t2=np.max(rmsds2)*r
    
    for i in range(len(rmsds)):
        if rmsds[i] < t1 and rmsds2[i] < t2:
            plt.bar(i, 0.5)
            print("KUA ", i)
            break


#ax3 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

#color = 'tab:orange'
#ax3.set_ylabel('energy', color=color)  # we already handled the x-label with ax1
#ax3.plot(detailed_indexes, opt_energies, '*', color=color, markersize=6)
#ax3.tick_params(axis='y', labelcolor=color)
#plt.plot(detailed_indexes, hit_rates, 'o', color=color, markersize=1)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title(sys.argv[1])
plt.show()

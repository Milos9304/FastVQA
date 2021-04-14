import matplotlib.pyplot as plt
import numpy as np
import sys

energies = []

with open(sys.argv[1], encoding='utf8') as f:
    for line in f:
        energy = float(line)
        energies.append(energy)
        
plt.plot(np.linspace(1,len(energies),len(energies)), energies, 'ro', markersize=1)
plt.show()

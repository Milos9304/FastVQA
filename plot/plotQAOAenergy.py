import matplotlib.pyplot as plt
import numpy as np

energies = []

with open('statsfile.txt', encoding='utf8') as f:
    for line in f:
        energy = float(line)
        energies.append(energy)
        
plt.plot(np.linspace(1,len(energies),len(energies)), energies, 'ro')
plt.show()
        
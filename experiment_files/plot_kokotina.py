import matplotlib.pyplot as plt
import numpy as np

with open('kokotina.txt') as f:
    lines = list(map(float, f.readlines()))

    print(lines[-1])
    ground_state = int(lines[-1])
    lines = lines[:-1]

x=np.linspace(0,len(lines)-1,len(lines),dtype=int)
bs=plt.bar(x, lines, width=0.1)
bs[ground_state].set_color('r')
plt.tight_layout()
plt.show()

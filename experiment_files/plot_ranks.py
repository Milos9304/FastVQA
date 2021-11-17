import csv
import numpy as np
import sys
import os
import glob
import subprocess
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from glob import glob

dims=[15,16,17,18,19]
hit_rates={}

for dim in dims:
    dim_i=dim-10
    
    iis = []

    #for filename in glob.glob(os.path.join(sys.argv[1],"statsfile*_"+str(dim)+".txt")):

    for i in range(128):
        iis.append(i)

    for r_dir in glob('rank_'+str(dim)):#('r*'):
        hit_rates[r_dir[5:]]=[]
        for i in range(128):

            filename = r_dir+"/statsfile_"+str(i)+"_"+str(dim)+".txt"
            with open(filename) as f:

                line = subprocess.check_output(['tail', '-1', filename])[:-1].decode("utf-8")
                splitted=line.split(' ')
                opt_energy=splitted[0]
                hit_rate=splitted[1]
                x_vect=splitted[2:]

                hit_rates[r_dir[5:]].append(float(hit_rate))
print("done")
scatter_plots=[]
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=.25)
ax.set_ylim([0, 1])

nbSamples = 5000
r_init = 15
ri=str(r_init)
colors = list(map(lambda x: 'green' if x > 1/nbSamples else 'red', hit_rates[ri]))
for i in range(len(iis)):
    scatterPlot, = plt.plot(iis[i], hit_rates[ri][i], 'o', markersize=1, color=colors[i])
    scatter_plots.append(scatterPlot)

text = ax.text(0.05, 1.05, "Rank " +ri+ "success rate with " + str(nbSamples) + " measures and cvar=1.75:  " + str(colors.count('green')/len(colors)), transform=ax.transAxes)


def update(val):

    global nbSamples
    global ri

    nbSamples = val;

    colors=list(map(lambda x: 'green' if x > 1/nbSamples else 'red', hit_rates[ri]))
    text.set_text("Rank " +ri+ "success rate with " + str(nbSamples) + " measures and cvar=1.75:  " + str(colors.count('green')/len(colors)))

    for i in range(len(scatter_plots)):
        scatter_plots[i].set_color(colors[i])
    fig.canvas.draw_idle()
    
    print(ri)
    
def update_r(r):

    global nbSamples
    global ri

    ri=str(r)

    colors=list(map(lambda x: 'green' if x > 1/nbSamples else 'red', hit_rates[ri]))
    text.set_text("Rank " +ri+ "success rate with " + str(nbSamples) + " measures and cvar=1.75:  " + str(colors.count('green')/len(colors)))

    for i in range(len(scatter_plots)):
        scatter_plots[i].set_ydata(hit_rates[ri][i])
        scatter_plots[i].set_color(colors[i])
    fig.canvas.draw_idle()
    
axcolor="lightgoldenrodyellow"
slider = Slider(
        ax=plt.axes([0.25,0.1,0.65,0.03],facecolor=axcolor),
        label="nbSamples",
        valmin=50,
        valmax=50000,
        valinit=nbSamples,
        valstep=1)
        
slider_r = Slider(
        ax=plt.axes([0.25,0.15,0.65,0.03],facecolor=axcolor),
        label="rank",
        valmin=dims[0],
        valmax=dims[-1],
        valinit=r_init,
        valstep=1)

slider.on_changed(update)
slider_r.on_changed(update_r)

plt.show()
#with open('paper_experiment/out_higher_dims_small_dims_matrices.csv', newline='') as f:



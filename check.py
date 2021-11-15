import csv
import numpy as np
import sys
import os
import glob
import subprocess
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

dim=15
dim_i=dim-10

results=[]
exclude=[]

with open('paper_experiment/out_higher_dims_small_dims.csv', newline='') as f:
    reader = csv.reader(f, delimiter=';')
    for row in reader:
        splitted=row[dim_i][1:].split(',')
        results.append(int(round(float(splitted[0]))))
        #print(splitted)
        splitted[1]=splitted[1][2:]
        splitted[-1]=splitted[-1][:-2]
        solution=list(map(float, splitted[1:]))
        for sol in solution:
            if sol < 0:
                solution=list(map(lambda x: x*(-1), solution))
                break
        for sol in solution:
            #print(sol)
            if sol > 1 or sol < 0:
                exclude.append(len(results)-1)
                break

#print(exclude)
#print(len(exclude))

"""with open('paper_experiment/out_higher_dims_small_dims_matrices.csv', newline='') as f:
  reader = csv.reader(f, delimiter=';')
  for row in reader:
    matrices.append(row)

print(len(matrices))

a=np.empty((10,180))
for m in matrices:
    
    a[i,:]=np.array(list(map(int, r[1:-1].split(','))))
    i+=1

res=np.matmul(np.array([0,0,1,0,0,0,0,0,0,0]),a)
print(res)
print()
print(np.sum(np.square(res)))
"""
#i = 0
success=0
fail=0
hit_rates=[]
iis = []

#for filename in glob.glob(os.path.join(sys.argv[1],"statsfile*_"+str(dim)+".txt")):
for i in range(128):
    if i in exclude:
        continue

    iis.append(i)

    filename2 = "experiment_files/just_picking_second_highest/statsfile_"+str(i)+"_"+str(dim)+".txt"
    filename = "experiment_files/statsfile_"+str(i)+"_"+str(dim)+".txt"
    with open(filename) as f:
        line = subprocess.check_output(['tail', '-1', filename])[:-1].decode("utf-8")
        splitted=line.split(' ')
        opt_energy=splitted[0]
        hit_rate=splitted[1]
        x_vect=splitted[2:]

        hit_rates.append(float(hit_rate))

        with open(filename) as f:
            line = subprocess.check_output(['tail', '-1', filename2])[:-1].decode("utf-8")
            opt2=line.split(' ')[0]

            if results[i] == int(round(float(opt_energy))):
                print(i, ": ", results[i],"=",int(round(float(opt_energy))), "    old: ",int(round(float(opt2))) )
                success+=1
            else:
                if results[i] > int(round(float(opt_energy))):
                    print("dopici")
                print(i, ": ", results[i],"|",int(round(float(opt_energy))), "    old: ",int(round(float(opt2))) )
                fail+=1
print("Success rate = ", float(success)/(success+fail))

print(hit_rates)

scatter_plots=[]
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=.25)

nbSamplesInit=5000;
colors = list(map(lambda x: 'green' if x > 1/nbSamplesInit else 'red', hit_rates))
for i in range(len(iis)):
    scatterPlot, = plt.plot(iis[i], hit_rates[i], 'o', markersize=1, color=colors[i])
    scatter_plots.append(scatterPlot)

text = ax.text(0.05, 1.05, "Success rate with " + str(nbSamplesInit) + " measures:  " + str(colors.count('green')/len(colors)), transform=ax.transAxes)

unsucessful = []
for i in range(len(colors)):
    if colors[i] == "red":
        unsucessful.append(iis[i])

print("Unsuccessful: ", unsucessful)

def update(val):
    #scatterPlot.set_color(list(map(lambda x: 'green' if x > 1/val else 'red', hit_rates)))
    colors=list(map(lambda x: 'green' if x > 1/val else 'red', hit_rates))
    text.set_text("Success rate with " + str(val) + " measures:  " + str(colors.count('green')/len(colors)))

    unsucessful = []
    for i in range(len(colors)):
        if colors[i] == "red":
            unsucessful.append(iis[i])

    print("Unsuccessful: ", unsucessful)


    for i in range(len(scatter_plots)):
        scatter_plots[i].set_color(colors[i])
    fig.canvas.draw_idle()

axcolor="lightgoldenrodyellow"
slider = Slider(
        ax=plt.axes([0.25,0.1,0.65,0.03],facecolor=axcolor),
        label="nbSamples",
        valmin=50,
        valmax=50000,
        valinit=nbSamplesInit,
        valstep=1)

slider.on_changed(update)

plt.show()
#with open('paper_experiment/out_higher_dims_small_dims_matrices.csv', newline='') as f:



import csv
import numpy as np
import sys
import os
import glob
import subprocess


dim=15
dim_i=dim-10

results=[]
with open('paper_experiment/out_higher_dims_small_dims.csv', newline='') as f:
    reader = csv.reader(f, delimiter=';')
    for row in reader:
        splitted=row[dim_i][1:].split(',')
        results.append(int(round(float(splitted[0]))))
        print(splitted)
        splitted[1]=splitted[1][2:]
        splitted[-1]=splitted[-1][:-2]
        print(list(map(float, splitted[1:])))
        if list(map(float, splitted[1:])).count(1.) > 1: #+ list(map(float, splitted[1:])).count(-1.) > 1:
            print(row)
        
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
#for filename in glob.glob(os.path.join(sys.argv[1],"statsfile*_"+str(dim)+".txt")):
for i in range(128):
    filename2 = "experiment_files/just_picking_second_highest/statsfile_"+str(i)+"_"+str(dim)+".txt"
    filename = "experiment_files/statsfile_"+str(i)+"_"+str(dim)+".txt"
    with open(filename) as f:
        line = subprocess.check_output(['tail', '-1', filename])[:-1].decode("utf-8")
        splitted=line.split(' ')
        opt_energy=splitted[0]
        hit_rate=splitted[1]
        x_vect=splitted[2:]

        with open(filename) as f:
            line = subprocess.check_output(['tail', '-1', filename2])[:-1].decode("utf-8")
            opt2=line.split(' ')[0]

            if results[i] == int(round(float(opt_energy))):
                print(results[i],"=",int(round(float(opt_energy))), "    old: ",int(round(float(opt2))) )
                success+=1
            else:
                if results[i] > int(round(float(opt_energy))):
                    print("dopici")
                print(results[i],"|",int(round(float(opt_energy))), "    old: ",int(round(float(opt2))) )
                fail+=1
print("Success rate = ", float(success)/(success+fail))

#with open('paper_experiment/out_higher_dims_small_dims_matrices.csv', newline='') as f:


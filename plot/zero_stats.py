import sys
import glob
import os
import subprocess

counter=0
zero_sols=0

for filename in glob.glob(os.path.join(sys.argv[1],"statsfile*")):
    with open(filename) as f:
        line = subprocess.check_output(['tail', '-1', filename])[:-1].decode("utf-8")
        splitted=line.split(' ')
        opt_energy=splitted[1]
        #print(opt_energy)
        if float(opt_energy) == 0:
            zero_sols+=1
        counter+=1

print("Ratio of zero solutions: ", float(zero_sols)/counter)

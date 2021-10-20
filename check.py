import csv
import numpy as np
with open('paper_experiment/out_higher_dims_small_dims_matrices.csv', newline='') as f:
  reader = csv.reader(f, delimiter=';')
  for row in reader:
    rr=row
    break

print(len(rr))

a=np.empty((10,180))
i=0
for r in rr:
    if i == 10:
        break
    a[i,:]=np.array(list(map(int, r[1:-1].split(','))))
    i+=1

res=np.matmul(np.array([0,0,1,0,0,0,0,0,0,0]),a)
print(res)
print()
print(np.sum(np.square(res)))

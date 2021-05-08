#!/usr/bin/env python
# coding: utf-8

# In[1]:

from fpylll import FPLLL, IntegerMatrix, LLL, SVP
import os
from glob import glob

FPLLL.set_random_seed(1997)

save_directory = "littleSombrero"


# In[ ]:
q=65537
d_full = 180
rank_reduced_min = 10
rank_reduced_max = 40

#num_instances = 20

for rank in range(rank_reduced_min, rank_reduced_max+1):
    path = "./"+save_directory+"/rank_"+str(rank)
    
    if not os.path.exists(path):
        os.mkdir(path)
    else:
        files = glob(path+"/*")
        for f in files:
            os.remove(f)
        
for i in range (20,50):#(0, num_instances):

    with open(save_directory+"/"+str(i)+"_info.txt", 'w') as f_info:  
    
        print("Instance " + str(i))

        A = LLL.reduction(IntegerMatrix.random(d_full, "qary", k=d_full//2, q=q))

        for rank in range(rank_reduced_min, rank_reduced_max+1):
            lattice_name = "q_d_"+str(d_full)+"_q_"+str(q)+"_r_"+str(rank)+"_i_"+str(i)
            B = A[:rank]

            print([norm(LLL.reduction(B)[0]), norm(vector(SVP.shortest_vector(B))).n()], file=f_info)
            
            with open(save_directory+"/rank_"+str(rank)+"/"+lattice_name, 'w') as f:  
                print(B, file=f)


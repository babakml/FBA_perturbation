import os
import cbmpy as cbm
import numpy as np 
import random
 


modelDir = os.getcwd()                                                               
cmod=cbm.CBRead.readSBML3FBC('mut-chem.xml', modelDir)
cmod.createGeneAssociationsFromAnnotations()



#fva
f, n = cbm.FluxVariabilityAnalysis(cmod, optPercentage=100)


#f[:,0]: opt, f[:,2] : min, f[:,3]: max, f[:,4]: dif




#performing the perturbation
 
prt_ind =[]

leng = range(709)

 
for i in leng:
    if f[i,4] > 0.000001:
        prt_ind.append(i)
        

fba_sol = [] #containing flux distributions after perturbation

rand_val = [] #random values used to fix each reaction
r=range(10)
for i in prt_ind:
    min_org = f[i,2]
    max_org = f[i,3]
    
    for j in r:
        rand = random.uniform(float(f[i,2]),float(f[i,3]))
        cmod.setReactionBounds(n[i], rand,rand)
        ff, nn = cbm.FluxVariabilityAnalysis(cmod, optPercentage=100)
        fba_sol.append(ff[:,0])
        rand_val.append(rand)
        cmod.setReactionBounds(n[i], min_org,max_org)
        
        
    np.hstack((fba_sol,ff[:,0]))
fba_np = np.array(fba_sol)          
np.savetxt("fba_np.csv", fba_np, delimiter=",")
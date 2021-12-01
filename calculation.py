import os
import cbmpy as cbm
import numpy as np
import random


modelDir = os.getcwd()                                                               
cmod=cbm.CBRead.readSBML3FBC('mut-chem.xml', modelDir)
cmod.createGeneAssociationsFromAnnotations()



#fva
f, n = cbm.FluxVariabilityAnalysis(cmod, optPercentage=100)


#finding reactions that got significantly affected by perturbation

fba_np = np.genfromtxt('fba_np_test.csv', delimiter=',')
rand_np = np.genfromtxt('rand_val_test.csv', delimiter=',')
r=range(10)

opt_som = []
l=1
c=1

for j in range(len((fba_np))):
    for i in range(len((fba_np[0]))):
        new_val_abs = abs(fba_np[j, i])
        dif_abs = abs(new_val_abs - f[i,0])
        tol1 = 5*abs(f[i,0])/100
        tol2 = 0.000001
        if abs(f[i,0]) !=0:
            tol = tol1
        else:
                tol = tol2
                
        if dif_abs > tol:
            opt_som.append([c,l,i,fba_np[j, i],f[i,0],tol])
        
        if i == 708:
            l = l+1
        if j/c == 10:
            c=c+1
            l=1
        

            
            
            
            
np.savetxt("final.csv", opt_som, delimiter=",")
        
        
            
    
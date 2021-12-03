import os
import cbmpy as cbm
import numpy as np
import random


def calculation(model, pert_result):
    """calculates the number of significantly change flux values 
       and saves the results in `final.csv` which can be used for 
       statistical analysis either in Matlab or python.

       :param model: the model filename to be used (is assumed to 
                     be in the current directory)
       :type model: str
       
       :param pert_result: csv file of pertubations obtained by running `perturb`
       :type pert_result: str 
    """
    
    modelDir = os.getcwd()                                                               
    cmod=cbm.CBRead.readSBML3FBC(model, modelDir)
    cmod.createGeneAssociationsFromAnnotations()

    #fva
    f, n = cbm.FluxVariabilityAnalysis(cmod, optPercentage=100)

    #finding reactions that got significantly affected by perturbation
    fba_np = np.genfromtxt(pert_result, delimiter=',')
    #rand_np = np.genfromtxt('rand_val_test.csv', delimiter=',')
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
        
            if i == len((fba_np[0]))-1:
                l = l+1
            if j/c == 10:
                c=c+1
                l=1
        
    np.savetxt("final.csv", opt_som, delimiter=",")
    print (opt_som)

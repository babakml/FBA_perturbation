import os
import cbmpy as cbm
import numpy as np
import random


def calculation(model, pert_result, method):
    
    """calculates the number of significantly change flux values 
       and saves the results in `final.csv` which can be used for 
       statistical analysis either in Matlab or python.

       :param model: the model filename to be used (is assumed to 
                     be in the current directory)
       :type model: str
       
       :param pert_result: csv file of pertubations obtained by running `perturb`
       :type pert_result: str 
    
       :param method: the method to be used either 'fba' or 'fva'
       :type method: str
    
        Example:
    
        calculation(model = 'mut-chem.xml', pert_result = 'fba_sol_fin.csv', method = 'fba')     
    
    """
    
    modelDir = os.getcwd()                                                               
    cmod=cbm.CBRead.readSBML3FBC(model, modelDir)
    cmod.createGeneAssociationsFromAnnotations()

    #fva
    fva_res, names = cbm.FluxVariabilityAnalysis(cmod, optPercentage=100)
    #fba
    fba_sol = cbm.analyzeModel(cmod, return_lp_obj=True) 
    sol_dist = fba_sol.solution.get_values()

    #finding reactions that got significantly affected by perturbation
    fba_np = np.genfromtxt(pert_result, delimiter=',')
    
    #rand_np = np.genfromtxt('rand_val_test.csv', delimiter=',')
    
    rng = range(10) #defining the number of perturbations

    opt_som = []
    
    c_count = 1 #counter for the perturbed reactions
    l_count = 1 #counter for perturbation number
    
    
    #finding the significant flux changes and sorting the results => perturbed reaction, the number of perturbation, the affected reaction (significant         flux change)
    for j in range(len((fba_np))):
        for i in range(len((fba_np[0]))):

                        
            if method == 'fva':
                new_val = fba_np[j, i]
                dif_abs = abs(new_val - fva_res[i,0])
                
            if method == 'fba':
                new_val = fba_np[j, i]
                dif_abs = abs(new_val - sol_dist[i])
                
            
            #defining the tolerance for non-zero fluxe values (tol1) and zero flux values (tol2)
            tol1 = 5*abs(fva_res[i,0])/100
            tol2 = 0.000001
            if abs(fva_res[i,0]) !=0:
                tol = tol1
            else:
                tol = tol2
            
            #checking whether the flux change is larger thant the tolerance    
            if dif_abs > tol:
                opt_som.append([c_count, l_count, i+1, fba_np[j, i], fva_res[i,0], tol]) #appending the significantly affected reactions, the number of                                                                                              perturbation  and the respective flux values
        
            if i == len((fba_np[0]))-1:
                
                l_count = l_count +1
                
            if j/c_count == 10:
                
                c_count = c_count +1
                l_count = 1
        
    np.savetxt("final.csv", opt_som, delimiter=",")
    print (opt_som)

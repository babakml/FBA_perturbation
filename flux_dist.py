import os
import cbmpy as cbm
import numpy as np 
import random
 
def perturb(model, method):
    
    """performs the perturbation, collecting flux distribution using FBA or FVA, 
    
    and saves the results in `fba_sol_fin.csv` and `rand_val_fin.csv`. 
    The file `fba_sol_fin.csv` can be used either by Matlab or python to 
    perform the down stream analysis.
    
    :param model: the model to be loaded (assumed to be in current directory)
    :type model: str 
    
    :param method: the method to be used either 'fba' or 'fva'
    :type method: str
    
    Example:
    
    perturb(model = 'mut-chem.xml', method = 'fba')
    
    performs the perturbation on the E. faecalis mutant model (integrated with medium composition and metabolic data),
    and collects the flux distributions using FBA
    
    """
    
    modelDir = os.getcwd()                                                               
    cmod=cbm.CBRead.readSBML3FBC(model, modelDir)
    cmod.createGeneAssociationsFromAnnotations()

    #fva
    fva_res, names = cbm.FluxVariabilityAnalysis(cmod, optPercentage=100)


    #f[:,0]: opt, f[:,2] : min, f[:,3]: max, f[:,4]: dif




    #preparing for perturbation
    prt_ind =[]

    rxns=cmod.getReactionIds()

    leng = range(len(rxns))

    #filtering variable reactions with the interval size larger than 0.000001 to use for perturbation
    for ind in leng:
        if fva_res[ind,4] > 0.000001:
            prt_ind.append(ind) #appending the index of reactions that will be pertubed
        
    sol_dist = []
    flux_dist = [] #to save the flux distributions after perturbation
    rand_val = [] #to save the random values used to fix each reaction
    dist_all = np.empty((10,10), int)

    #performing perturbation

    rng = range(10) #defining the number of perturbations for each reaction
    
    for reac in prt_ind:
    
        min_org = fva_res[reac,2] #saving the original lower bound
        max_org = fva_res[reac,3] #saving the original upper bound
    
        for j in rng:
            
            rand = random.uniform(float(fva_res[reac,2]),float(fva_res[reac,3])) #creating a random value
            cmod.setReactionBounds(names[reac], rand,rand) #fixing the reaction at the random value
        
            #collecting flux distributions using FBA
            if method == 'fba':
                
                fba_sol = cbm.analyzeModel(cmod, return_lp_obj=True) #FBA
                sol_dist = fba_sol.solution.get_values()             #collecting flux distribution profile
                cmod.setReactionBounds(names[reac], min_org,max_org) #bringing the original upper and lower bound back
                flux_dist.append(sol_dist)
                rand_val.append(rand)
                

        
            #collecting flux distributions using FVA
            if method == 'fva':

                
                fva_res_new, name_n = cbm.FluxVariabilityAnalysis(cmod, optPercentage=100) #FVA
                sol_dist=fva_res_new[:,0]                                                  #collecting flux distribution profile
                cmod.setReactionBounds(names[reac], min_org,max_org)                       #bringing the original upper and lower bound back
                flux_dist.append(sol_dist)
                rand_val.append(rand)
                

    #saving the result
    np.savetxt("fba_sol_fin.csv", flux_dist, delimiter=",")
    np.savetxt("rand_val_fin.csv", rand_val, delimiter=",")
    print (flux_dist, rand_val)

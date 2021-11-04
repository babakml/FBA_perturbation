import os
import cbmpy as cbm
import numpy as np
import random
os.chdir('/Users/babak/Documents/model_ensamble/flux_distribution_steady_state/mut_chem_cbm/')   


modelDir = os.getcwd()                                                               
cmod=cbm.CBRead.readSBML3FBC('mut-chem.xml', modelDir)
cmod.createGeneAssociationsFromAnnotations()



#fva
f, n = cbm.FluxVariabilityAnalysis(cmod, optPercentage=100)


###finding various statistical measures

prt_ind =[]
leng = range(709)
for i in leng:
    if f[i,4] > 0.000001:
        prt_ind.append(i)
        
        
fva = np.count_nonzero(f[:,4])
stable = len(n) - fva

#finding how many reactions were affected by each reaction


optsom_np = np.genfromtxt('final.csv', delimiter=',')
#optsom_np = numpy.array(opt_som) 

mx = len(optsom_np)

ind = optsom_np[:,0]
ind_u = np.unique(ind)

num_aff = []
for i in ind_u:
    aff = np.count_nonzero(optsom_np[:,0] == i)
    num_aff.append([i, aff])

num_aff = np.array(num_aff)    
aff_avg = np.average(num_aff[:,1])/10
aff_avg_std = np.std(num_aff[:,1])
aff_max = max(num_aff[:,1])
aff_min = min(num_aff[:,1])



reac_list = cmod.getReactionIds()

reac_sen = []
for i in range(len(cmod.reactions)):
    i2 = i+1
    
    sen = np.count_nonzero(optsom_np[:,2] == i2)
    aff_ind = np.where(optsom_np[:,2]==i2)
    b=aff_ind[0].size==0
    if b is False:
        aff_ind = aff_ind[0]
        aff_u = np.unique(optsom_np[aff_ind[0],0])
        aff_num = len(aff_u)
        reac_sen.append([i2, reac_list[i], aff_u[0], sen])
    
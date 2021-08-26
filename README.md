# FBA_perturbation
contact: babak.loghmani@bioquant.uni-heidelberg.de

a toolbox to perturb FBA flux distribution in constraint-based metabolic networks

This toolbox is designed to discover the uncertainties in FBA flux distribution by perturbing the system using random values.
The toolbox is developed for MATLAB (mathworks.com) and uses cobratoolbox (https://opencobra.github.io/cobratoolbox/stable/) to collect flux distribution profiles.

Make sure you have cobratoolbox installed.
To obtain flux distribution profiles using 'fastFVA' function, IBM_CPLEX (IBM.com) must be installed.
To obtain flux distribution profiles using 'optimizeCbModel', either IBM_CPLEX or glpk can be used.

To perform perturbation: 
 
 [fbasol_fin, randval_fin] = perturb(model_n, solver_n) #enter model and solver name as strings
 
To perform statistical analysis:
 table = stat(model_n, per_result, print)

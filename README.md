# FBA_perturbation
contact: babak.loghmani@bioquant.uni-heidelberg.de

a toolbox to perturb FBA flux distribution in constraint-based metabolic networks

This toolbox is designed to discover the uncertainties in FBA flux distribution by perturbing the system using random values.
The toolbox is developed for [MATLAB](mathworks.com) and uses the [cobratoolbox](https://opencobra.github.io/cobratoolbox/stable/) to collect flux distribution profiles. A python version is also developed based o Pysces CBMpy (https://github.com/SystemsBioinformatics/cbmpy).

## Requirements

 * On MATLAB:
 
   * Make sure you have [cobratoolbox](https://opencobra.github.io/cobratoolbox/stable/installation.html) installed.
   * To obtain flux distribution profiles using the `fastFVA` function, [IBM\_CPLEX](https://www.ibm.com/support/pages/cplex-optimization-studio-v128) must be installed. The cobra toolbox works only for cplex versions up to 12.8. 
   * To obtain flux distribution profiles using `optimizeCbModel`, either IBM\_CPLEX or glpk can be used.



 * On python 2:
   * Make sure you have Pysces [CBMpy](https://pythonhosted.org/cbmpy/install_doc.html) installed.
   * Make sure you use [IBM_Cplex](https://www.ibm.com/support/pages/cplex-optimization-studio-v128) as the solver.

## Executing the Matlab functions

* To perform the perturbations use: 
 
 	`[fbasol_fin, randval_fin] = perturb(model_n, solver_n, method)`

 	enter model and solver and method of choice (fba or fva) name as strings.

* To find significant flux changes and sort the output:

	`[per_result] = calculation(model_n, fbasol_fin, randval_fin, method)`

	enter the model and the result of perturbation as the inputs, as well as the method of choice (fba or fva)
 
* To perform statistical analysis:

	`table = stat(model_n, per_result, print)`

* So for example: 

	`[fbasol_fin, randval_fin] = perturb('mut-chem.xml' , 'ibm_cplex', 'fba')`


* In case the method of use is fva, this will take quite a while to run and produce files:

  * `fbasol_fin` file contains significantly different flux values in a sorted format

  * `randval_fin` file contains the random values at which reactions were fixed.
  
* To find significant flux changes:

	* [per_result] = calculation('mut-chem.xml','fbasol_fin.dat','randval_fin.dat','fba')

  * `per_result` file contains all the significant flux changes in a sorted format. 

* Then to calculate the statistics for this model you'd run:

	`table = stat('mut-chem.xml', 'per_result', 'print')`

* Which would produce a result like:

		table =

		2×23 table

			out1           out2         out3         out4          out5                   out6                   out7                      out8                     out9                    out10                    out11                      out12                       out13              out14                out15                out16                   out17                    out18                    out19              out20             out21                    out22                    out23        
		____________    __________    ________    ___________    ________    ______________________________    _________    __________________________________    _________    ________________________________    __________    ____________________________________    ____________    _________________    _________________    __________________    _____________________    ______________________    ____________________    _________    ____________________    _____________________    _____________________

		'model name'    'variable'    'stable'    'sensitive'    'robust'    'affecting-avg(reaction-wise)'    'std-aff'    'affecting-avg(perturbation-wise)'    'std-aff'    'avg-sensitivity(reaction-wise)'    'std-sen'     'avg-sensitivity(perturbation-wise)'    'std-sen'       'max-sensitivity'    'min-sensitivity'    'robust with flux'    'robust without flux'    'robust-w-variability'    'avg-affected by ex'    'std-ex'     'max-affected by ex'    'max-affected by all'    'min-affected by all'
		'mut-ms.xml'    [     311]    [   398]    [      294]    [   415]    [                    150.1342]    [34.4665]    [                        123.3299]    [22.5815]    [                      152.1769]    [106.2194]    [                        1.2501e+03]    [1.0147e+03]    [            298]    [              2]    [             160]    [                255]    [                  17]    [          141.3600]    [30.1322]    [               213]    [           153.8000]    [            33.5000]

## Python files description

#Python 2.7

The file `flux_dist.py` performs the perturbation, collecting flux distribution using FVA, and saves the results in `fba_sol_fin.csv` and `rand_val_fin.csv`. The file `fba_sol_fin.csv` can be used by calculation function in python to perform the down stream analysis.

example:

	from flux_dist import perturb
	
	perturb(model='mut-chem.xml', method='fba')


The file `calculation.py` calculates the number of significantly change flux values and saves the results in `final.csv` which can be used for statistical analysis either in Matlab or python.

example:

	from calculation import calculation 
	
	calculation(model='mut-chem.xml', pert_result='fba_sol_fin.csv', method='fba')


The file `stat_asl.py` uses `final.csv` as input and performs different statistical analyses. This file is still to be completed and does not yet provide all the statistical measures that are provided by `stat.m` in Matlab. For now, to obtain the full list of statistical measures, the `final.csv` file can be used by stat.m in Matlab.

#Python 3

The files `flux_dist_p3.py`, `calculation_p3.py` and `stat_asl_p3.py` are adjusted for python 3. Make sure that the cplex version is compatible with your python 3 version. Pysces CBMpy does not support all the modules for glpk on python 3, so make sure that you have a compatible cplex version.

## List of models

* `wt_unconstrained`: E.faecalis model with no biological constraint
* `wt_medium`: E. faecalis model constrained with medium compistion data
* `wt-chem`: E. faecalis model constrained with medium composition data + metabolic data
* `wt_ms`: E. faecalis model constrained with medium composition data + metabolic data + proteome data
* `mut_unconstrained`:△glnA E.faecalis model with no biological constraint
* `mut-medium`: △glnA E. faecalis model constrained with medium compistion data
* `mut-chem`: △glnA E. faecalis model constrained with medium composition data + metabolic data
* `mut-ms`: △glnA E. faecalis model constrained with medium composition data + metabolic data + proteome data

## License
This project is licensed under the BSD license: 

	Copyright (c) 2021, Seyed Babak Loghmani
	All rights reserved. 
	
	Redistribution and use in source and binary forms, with or without 
	modification, are permitted provided that the following conditions are 
	met: 
	
	Redistributions of source code must retain the above copyright notice, 
	this list of conditions and the following disclaimer. Redistributions in 
	binary form must reproduce the above copyright notice, this list of 
	conditions and the following disclaimer in the documentation and/or 
	other materials provided with the distribution. THIS SOFTWARE IS 
	PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
	EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
	PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
	EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
	PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
	PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
	

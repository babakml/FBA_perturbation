# FBA_perturbation
contact: babak.loghmani@bioquant.uni-heidelberg.de

a toolbox to perturb FBA flux distribution in constraint-based metabolic networks

This toolbox is designed to discover the uncertainties in FBA flux distribution by perturbing the system using random values.
The toolbox is developed for [MATLAB](mathworks.com) and uses the [cobratoolbox](https://opencobra.github.io/cobratoolbox/stable/) to collect flux distribution profiles.

## Requirements

 * Make sure you have cobratoolbox installed.
 * To obtain flux distribution profiles using the `fastFVA` function, [IBM\_CPLEX](IBM.com) must be installed. The cobra toolbox works only for cplex versions up to 12.8. 
 * To obtain flux distribution profiles using `optimizeCbModel`, either IBM\_CPLEX or glpk can be used.

## Usage
To perform the perturbations use: 
 
 	[fbasol_fin, randval_fin] = perturb(model_n, solver_n) #enter model and solver name as strings
 
To perform statistical analysis:

	table = stat(model_n, per_result, print)

So for example: 

	[fbasol_fin, randval_fin] = perturb('mut-ms.xml' , 'ibm_cplex')

this will take quite a long time to run and produce files:

* `per_result*` files the file containing significantly different flux values in a sorted format

Then to calculate the statistics for this model you'd run:

	table = stat('mut-ms.xml', 'per_result', 'print')

Which would produce a result like:

	table =

  	2×23 table

        out1           out2         out3         out4          out5                   out6                   out7                      out8                     out9                    out10                    out11                      out12                       out13              out14                out15                out16                   out17                    out18                    out19              out20             out21                    out22                    out23        
    ____________    __________    ________    ___________    ________    ______________________________    _________    __________________________________    _________    ________________________________    __________    ____________________________________    ____________    _________________    _________________    __________________    _____________________    ______________________    ____________________    _________    ____________________    _____________________    _____________________

    'model name'    'variable'    'stable'    'sensitive'    'robust'    'affecting-avg(reaction-wise)'    'std-aff'    'affecting-avg(perturbation-wise)'    'std-aff'    'avg-sensitivity(reaction-wise)'    'std-sen'     'avg-sensitivity(perturbation-wise)'    'std-sen'       'max-sensitivity'    'min-sensitivity'    'robust with flux'    'robust without flux'    'robust-w-variability'    'avg-affected by ex'    'std-ex'     'max-affected by ex'    'max-affected by all'    'min-affected by all'
    'mut-ms.xml'    [     311]    [   398]    [      294]    [   415]    [                    150.1342]    [34.4665]    [                        123.3299]    [22.5815]    [                      152.1769]    [106.2194]    [                        1.2501e+03]    [1.0147e+03]    [            298]    [              2]    [             160]    [                255]    [                  17]    [          141.3600]    [30.1322]    [               213]    [           153.8000]    [            33.5000]


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
	

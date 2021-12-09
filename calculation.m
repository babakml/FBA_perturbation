function [per_result] = calculation(model_n, fbasol_fin, randval_fin, method)

% calculation finds the optimal values in flux distribution profiles following perturbation that
% are significantly different (>0.05 0r >0.000001) from the original FBA/FVA flux distribution
% profiles.

%USAGE: 
%    [per_result] = calculation(model, fbasol_fin, randval_fin)

%INPUT

%model: metabolic model in SBML format

%fbasol_fin: flux distribution profiles following perturbation

%randval_fin: random values at which variable reactions were fixed

%method: 'fba' or 'fva', to determine the method of use to obtain flux

%OUTPUT:

%per_result: significantly different flux values in a sorted format

%Example:

%calculation [per_result] = calculation('mut-chem.xml','fbasol_fin.dat','randval_fin.dat','fba')

% Authors:

% Seyed Babak Loghmani

% Last updated: August 2021

model = readCbModel(model_n);
load(fbasol_fin)
load(randval_fin)

%FVA
[minFluxF1, maxFluxF1, optsol, ret, fbasol, fvamin, fvamax, statussolmin, statussolmax] = fastFVA(model);

if method == 'fva'
    opt_sol = fbasol;
elseif method == 'fba'
    sol_dist = optimizeCbModel(model);
    opt_sol = sol_dist.v;
end

%finding reactions with fva intervals larger that 0.000001
fva_n = maxFluxF1 - minFluxF1;
r=fva_n > 0.000001;
fva_n_f=find(r);

n = numel(model.rxns);
opt_som = [];
sols = [];
per_result = [];
%opt_sol = fbasol(1:n);
opt_sol_abs = abs(opt_sol);

% defining the 5% tolerance
tol = 5*opt_sol_abs/100;
s = size(randval_fin);
f = 0;
rxn_n = numel(model.rxns);


%finding and sorting significant flux changes, assigning their indices and random values in which they
%got fixed

num_f_p = numel(fva_n_f)/50;
num_f = ceil(num_f_p);
rg = 1:num_f;
c=1;
file_names = [];

for i = 1:numel(fva_n_f)
    k = i-1;
    
    l = (k*n) + 1;%lowe3r end of the data set
    h = (i*n);%higher end of the data set

    
if i == rg(c)*50%counter to save the data regarding ecery 50 reactions at once
    f = 1
    c = c+1;
elseif c == max(rg)
    
    if i == numel(fva_n_f)
       f = 1
    end 
end


if i ==1
    
    sols = fbasol_fin(1:rxn_n,1:10);%extracting the data for the first reaction
else
    sols = fbasol_fin(l:h,1:10);%extracting the data for the other reactions
end

for j = 1:10

    dif = opt_sol - sols(:,j);
    dif_abs = abs(dif);
    for k = 1:rxn_n
        
        %finding changes in flux values larger than 0.000001 for
        %reactions whose original value was 0
        
        if tol(k)==0
           if dif_abs(k) > 0.000001
              opt_som = [i,j, k, sols(k,j), opt_sol(k),tol(k), randval_fin(i,j)];
              per_result = [per_result; opt_som];
            
           end
           
        else
            
           %finding changes in flux values larger than 5% tolernace
           
           if dif_abs(k) > tol(k)
               
               %saving significantly different flux values in per_result
               
               opt_som = [i,j, k, sols(k,j), opt_sol(k),tol(k), randval_fin(i,j)];
               per_result = [per_result; opt_som];
           end

        end
        
    end
end

%saving sorted flux values in multiple files for the sake of memory usage

if f ==1
filename = sprintf('%s%d','per_result_new',i, '.dat');
save(filename,'per_result','-ascii','-double');

f_str = string(filename)
file_names = [file_names, f_str];

f=0;
clear per_result

per_result =[];

    end
end


%merging all the files into one comprehensive file
for i = 1:numel(file_names);
    
    file = load(file_names(i));
    per_result = [per_result;file];
end
save('per_result','per_result','-ascii','-double');
            
        
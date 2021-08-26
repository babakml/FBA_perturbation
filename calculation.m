function [per_result] = calculation(model, fbasol_fin, randval_fin)

% calculation finds the optimal values in flux distribution profiles following perturbation that
% are significantly different (>0.05 0r >0.000001) from the original FBA/FVA flux distribution
% profiles.

%USAGE: [per_result] = calculation(model, fbasol_fin, randval_fin)

%INPUT

%model: metabolic model in SBML format

%fbasol_fin: flux distribution profiles following perturbation

%randval_fin: random values at which variable reactions were fixed

%OUTPUT:

%per_result: significantly different flux values in a sorted format



% Authors:

% Seyed Babak Loghmani

% Last updated: August 2021



[minFluxF1, maxFluxF1, optsol, ret, fbasol, fvamin, fvamax, statussolmin, statussolmax] = fastFVA(model);

fva_n = maxFluxF1 - minFluxF1;
r=fva_n > 0.000001;
fva_n_f=find(r);

n = numel(model.rxns);
opt_som = [];
sols = [];
per_result = [];
opt_sol = fbasol(1:n);
opt_sol_abs = abs(opt_sol);
tol = 5*opt_sol_abs/100;
s = size(randval_fin);
f = 0;
rxn_n = numel(model.rxns);


%sorting reactions, assigning their indices and random values in which they
%got fixed

num_f_p = numel(fva_n_f)/50;
num_f = ceil(num_f_p);
rg = 1:num_f;
c=1;

file_names = [];

for i = 1:numel(fva_n_f)
    k = i-1;
    l = (k*n) + 1;

    h = (i*n);

if i == rg(c)*50
    f = 1
    c = c+1;
elseif c == max(rg)
    if i == numel(fva_n_f)
       f = 1
    end
    
end


if i ==1
    
    sols = fbasol_fin(1:rxn_n,1:10);
else
    sols = fbasol_fin(l:h,1:10);
end

for j = 1:10
    sols_abs = abs(sols);
    dif = opt_sol_abs - sols_abs(:,j);
    dif_abs = abs(dif);
    for k = 1:rxn_n
        
        if dif_abs(k) > tol(k)
            if tol(k)==0
                if sols_abs(k,j) > 0.000001
                     opt_som = [i,j, k, sols(k,j), opt_sol(k),tol(k), randval_fin(i,j)];
                     per_result = [per_result; opt_som];
            
                end
            else
                opt_som = [i,j, k, sols(k,j), opt_sol(k),tol(k), randval_fin(i,j)];
                     per_result = [per_result; opt_som];
            end
            
            
            elseif dif_abs(k) ~= 0

        end
        
    end
end
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

for i = 1:numel(file_names);
    
    file = load(file_names(i));
    per_result = [per_result;file];
end
save('per_result','per_result','-ascii','-double');
            
        
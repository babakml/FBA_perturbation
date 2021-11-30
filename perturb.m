function [fbasol_fin, randval_fin] = perturb(model_n, solver_n)

%pertub function discovers the effect of random perturbation in
%FBA/FVA results. It selects 10 random values in the feasible interval of
%each variable reaction (with an interval larger than 0.000001), fixes
%reactions in those values and collect the flux distribution profile of the
%whole metabolic network

%USAGE

%    [fbasol_fin, randval_fin] = perturb(model_n, solver_n)


%INPUT

% model_n: metabolic model in SBML format
% solver_n: solver name

% OUTPUT

% fbasol_fin: flux distribution profiles following 10 perturbation in each
% variable reaction

% randval_fin: random values at which variable reactions were fixed

% Authors:

% Seyed Babak Loghmani

% Last updated: August 2021




initCobraToolbox;

%defining the solver
changeCobraSolver(solver_n);
model = readCbModel(model_n);

[minFluxF1, maxFluxF1, optsol, ret, fbasol, fvamin, fvamax, statussolmin, statussolmax] = fastFVA(model);



end
da = maxFluxF1;
db = minFluxF1;
fbasol_fin = [];
randval_fin =[];
dif = (da - db);

rxn_n = numel(model.rxns);

%runing the model while fixing each reaction in 10 random values
for i =1:rxn_n
    fbasol_m = [];
            randval_m = [];
            dif = (da(i) - db(i));
    if dif > 0.000001       
    if dif ~= 0
        for j = 1:10
        
        %finding a random value
        randval = random('unif',db(i),da(i));
        %saving the original lower and upper bounds
        pu = model.ub(i);
        pl = model.lb(i);
        %fixing the reaction at the randomly selected value
        model.ub(i) = randval;
        model.lb(i) = randval;
        i
        da(i)
        db(i)
        da(i) - db(i)
        %collecting flux distributions using FBA
        
         sol_dist = optimizeCbModel(model);
         fbasol_m = [fbasol_m, sol_dist.v];
         
        %flux distributions using FVA
        
        %[minFluxF1, maxFluxF1, optsol, ret, fbasol, fvamin, fvamax, statussolmin, statussolmax] = fastFVA(model);
        %fbasol_m = [fbasol_m, fbasol];
        
        %saving the random value in randval_m
        randval_m = [randval_m, randval];
        
        
        %bringing the original uppre and lower bound back
        model.ub(i) = pu;
        model.lb(i) = pl;

        end
        
        %saving the flux distribution profiles & random values for 10 perturbation in one
        %reaction in fbasol_fin & randval_fin, respectively
        fbasol_fin = [fbasol_fin; fbasol_m];
        randval_fin = [randval_fin; randval_m];

        
        
        
    end
    end
end

save fbasol_fin.dat fbasol_fin -ascii -double
save randval_fin.dat randval_fin -ascii -double

calculation(model, fbasol_fin, randval_fin)
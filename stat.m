function table = stat(model_n, per_result, print)

% stat function calculates various statiscal measures on the flux values
% following perturbation in metabolic network

%USAGE: table = stat(model_name, per_result, print)


%INPUT

%model_n: metabolic model in SBML format

%per_result: significantly different flux values in a sorted format

%print: printing the statistics in an excel file (default: false)

%OUTPUT:

%table: statistical categories obtained from analysis of flux distribution
%profiles

%Example:

%table = stat('mut-chem.xml', 'per_result') using the perturbation results
%from MATLAB version

%table = stat('mut-chem.xml', 'final.csv') using the perturbation result
%from the python version



% Authors:

% Seyed Babak Loghmani

% Last updated: August 2021


if (nargin < 3)
    print = false;
end

model = readCbModel(model_n);
final = load(per_result);

[minFluxF1, maxFluxF1, optsol, ret, fbasol, fvamin, fvamax, statussolmin, statussolmax] = fastFVA(model);
fva_n = maxFluxF1 - minFluxF1;
fva_ind = find(fva_n);
stable_ind = find(~fva_n);
fva=numel(find(fva_n));
stable = numel(stable_ind);
r=fva_n > 0.000001;
fva_n_f=find(r);

rxn_n = numel(model.rxns);
rxn_n2 = rxn_n +1;


reac =[];
reac_no=[];
rand_no1=[];
rand = [];
rand_no = [];
a = final(:,1);
sz = size(final);
mx = final(sz(1));
robust_var= [];
%finding how many reactions were affected by each reaction
for i = 1:mx
r_a = find(ismember(a, i));
r_a_min = min(r_a);
r_a_max = max(r_a);
 r_a_n=numel(r_a);
    r_af = final(r_a_min:r_a_max,2);
    for j= 1:10
        r_af_i = find(ismember(r_af, j));
        d1 = numel(r_af_i);
        rand =[i,j,d1];
        rand_no =[rand_no;rand];
    end
    r_a_al = final(r_a_min:r_a_max,3);
    r_a_al_u=unique(r_a_al);
    r_a_num = numel(r_a_al_u);
    r_a_n2 = r_a_n/10;
    reac =[i,r_a_n,r_a_n2,r_a_num];
    reac_no =[reac_no;reac]; %size equal to fva result
  
end
re = reac_no(:,3);
re2 = reac_no(:,4);
mean_re = mean(re); %affecting avg(perturbation wise)
mean_re2 = mean(re2); %affecting avg(reaction wise)
std_dev_meanre= std( re ); %standard deviation
std_dev_meanre2= std( re2 ); %standard deviation
min_re = min(re); %affecting min
max_re = max(re); %affecting max

%finding sensitivity of each reaction
%column= affected, row = affecting
reac_sen =[];
ou = [];
robust=[];
pul = 0;
mat = cell(rxn_n2:rxn_n2)
mat(2:rxn_n2) = model.rxns;
mat(1,2:rxn_n2)= model.rxns;
for i = 1:rxn_n
    rc = find(ismember(final(:,3), i));
    t = numel(rc);
    ou = [];
   
    if isempty(rc) == 1
        robust=[robust; i];
    else
    for j =1:numel(rc)
    
            ind = rc(j);
            
            ou =[ou; final(ind, 1)];
            
    end
    ou_u = unique(ou);
    for s = 1:numel(ou_u)
        s_id = find(ismember(ou, ou_u(s)));
        s_id_u = numel(s_id);
        nb = ou_u(s);
       i2 = fva_n_f(nb)+1;
       i3 = i+1;
     
       mat{i2, i3} = s_id_u;
    end
    gu = numel(ou);
    n_ou = numel(ou_u);
    reac_sen = [reac_sen; i, model.rxns(i), n_ou, t]; %sensitive to perturbation
    
    end
end

sz_sen = size(reac_sen);
szs = sz_sen(1);


sensitivity=[];
for i = 1:szs
u=reac_sen{i,1};
sensitivity=[sensitivity;u];
end

num_sen = numel(sensitivity);
num_robust = numel(robust);

    for i = 1:rxn_n2
        mat{i,i}=0;
    end

    %average and maximum and minimum number of sensitivities
    
    sen1 = cell2mat(reac_sen(:,3));
    sen2 = cell2mat(reac_sen(:,4));
    msen = mean(sen1); %average sensitivity
    msen2 = mean(sen2); %average sensitivity(perturbation wise)
    std_dev_sen= std( sen1 ); %standard deviation
    std_dev_sen2= std( sen2 ); %standard deviation 2
    maxsen = max(sen1); %max sensitivity
    minsen = min(sen1); %min sensitivity
    
    
    %robusts with and without flux

    du= find(fbasol);
    pak=find(ismember(du,robust));
    pak1=du(pak);
    robust_wf=model.rxns(pak1); %robust with flux
    
    %check with FVA result to see if they had variabilty or not
    
    robust_var = [];
    for i = 1:numel(robust)
        ek = find(ismember(fva_ind, robust(i)));
        if numel(ek) == 1
            robust_var = [robust_var; robust(i)];
        end
    end
         
    num_robust_var = numel(robust_var);
    num_robust_wf = numel(robust_wf);
    num_robust_wof = numel(model.rxns) - num_robust_wf - num_sen;

        
    %variability caused by exchange reactions
    num_imp = [];
    for i = 2:rxn_n2
        ji = 0;
        e = findstr(mat{i,1}, 'Ex');
        if numel(e) ==1
            g = mat(i,:);
            for j = 2:rxn_n2
                if j~=i
                    gg = g{j};
                    if numel(gg) == 1
                        ji = ji + 1;
                    end
                end
            end
        end
        num_imp = [num_imp; i-1,ji];
    end
    
    %maximum and average number that reactions got affected by exchange
    %reactions
nimp = num_imp(:,2);
f_nimp = find(nimp);
nimp_e=nimp(f_nimp);
avg_ex = mean(nimp_e); %avg affected by Ex
std_dev_ex= std( nimp_e ); %standard deviation
max_ex = max(nimp_e); %max affected by Ex

out={'model name','variable','stable','sensitive','robust','affecting-avg(reaction-wise)','std-aff','affecting-avg(perturbation-wise)','std-aff','avg-sensitivity(reaction-wise)','std-sen','avg-sensitivity(perturbation-wise)','std-sen','max-sensitivity','min-sensitivity','robust with flux','robust without flux','robust-w-variability','avg-affected by ex','std-ex','max-affected by ex','max-affected by all','min-affected by all';model.description,fva,stable,num_sen,num_robust,mean_re2,std_dev_meanre2,mean_re,std_dev_meanre,msen,std_dev_sen,msen2,std_dev_sen2,maxsen,minsen,num_robust_wf,num_robust_wof,num_robust_var,avg_ex,std_dev_ex,max_ex,max_re,min_re};
table=cell2table(out);
if (print)
    writetable(table,'statistics.xlsx','WriteVariableNames',0)
end





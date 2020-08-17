function model_delNet = addNetDeltaFlux(model_del, ori_rxnNames)
% Function for taking a model that has been created for delta and adding
% net flux variables for each reaction. The Net Delta Flux is sum of
% forward and reverse reactions. Also we add lower and upper bounds for Net
% Delta Fluxes. 
%
%
% Inputs:
% 1. model_del - Delta model with bounds already constrained. Please refer
% to createDeltaModel.m first before using this
% 2. ori_rxnNames - Original reaction names as defined by a GEM.
%
%
% Outputs:
% 1. model_delNet - Model with Net flux for reversible reactions


% Step 1: Identify reversible reactions and add Netflux variable to combine
% the overall material flow thru forward and reverse reactions.

model = model_del;

[num_eqs, num_vars] = size(model.S);

model.A = sparse(model.S);
model.varNames = model.rxns;
model.vtype = repelem('C', num_vars)';
model.varlb = model.lb; %lb in MILP
model.varub = model.ub; %ub in MILP
model.obj = model.c;
model.eqNames = model.mets;
model.eqtype = repelem('=', num_eqs)'; %csense in MILP
model.rhs = model.b;
model.modelsense = 'max';
model.x0 = [];

num_rev = numel(find(model.rev))/2;
num_irrev = size(model.S,2)-num_rev;
rev_idx = model.match(size(model.S,2)-num_rev+1:end);

for i = 1:num_irrev
    model.A(num_eqs+1, num_vars+1) = -1;
    model.A(num_eqs+1,i) = 1;
    model.eqtype(num_eqs+1) = '=';
    model.rhs(num_eqs+1) = 0;
    model.eqNames{(num_eqs+1),1} = strcat('NF_', ori_rxnNames{i,1});
    model.varNames{(num_vars+1),1} = strcat('NF_', ori_rxnNames{i,1});
    model.varlb(num_vars+1) = model.lb(i);
    model.varub(num_vars+1) = model.ub(i);
    model.vtype(num_vars+1) = 'C';
    model.obj(num_vars+1) = 0;
    if ~isempty(find(rev_idx==i))
        model.A(num_eqs+1,num_irrev+find(rev_idx==i)) = -1;
        %model.varlb(num_vars+1) = model.lb(i)-(abs((epsilon/100)*Flux(i)));
        %model.varub(num_vars+1) = model.ub(i)+(abs((epsilon/100)*Flux(i)));
        model.varlb(num_vars+1) = -(model.ub(i)-model.lb(i));
        model.varub(num_vars+1) = model.ub(i)-model.lb(i);
    end
    [num_eqs, num_vars] = size(model.A);
end
model_delNet = model;
end
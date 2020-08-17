function [model_del, nochange_idx] = createDeltaModel_rev(model, essential_idx, zeroflow_idx,maxflux_val)
% For the inputted model and index of essential reactions create a delta 
% model. In this version, conversion to a reversible model is skipped to
% aid in gene based model. In the next step, the bounds of the variables in 
% the delta model will be set based on the minflux_val. Finally the bounds
% of the essential reactions will be set to 0 as these reactions will not
% be allowed to change. These are known zero delta reaction fluxes. Also we
% will set zero flow flux reaction bounds to zero. 


% Step 1: Convert to irreversible model (based on Cobra's implementation
   %model_tmp = convertToIrreversible(model);
   model_tmp = model;
% Step 2: Look at the minflux_val and convert all lower bounds and upper
% bounds 

bound_val = 10^ceil(log10(maxflux_val));

model_tmp.lb = -(bound_val*ones(size(model_tmp.ub,1),1)+bound_val*ones(size(model_tmp.lb,1),1));
model_tmp.ub = (bound_val*ones(size(model_tmp.ub,1),1)+bound_val*ones(size(model_tmp.lb,1),1));

% Step 3: Convert essential reaction bounds to zero
%rev_essential_idx = nonzeros(model_tmp.match(essential_idx));
%essential_idx = union(essential_idx, rev_essential_idx);

model_tmp.lb(essential_idx) = 0;
model_tmp.ub(essential_idx) = 0;

%rev_zeroflow_idx = nonzeros(model_tmp.match(zeroflow_idx));
%zeroflow_idx = union(zeroflow_idx, rev_zeroflow_idx);

model_tmp.lb(zeroflow_idx) = 0;
model_tmp.ub(zeroflow_idx) = 0;

nochange_idx = union(essential_idx, zeroflow_idx);
model_del = model_tmp;

% Step 4: make any objective function coefficents zero
model_del.c(find(model_del.c)) = 0;

[num_eqs, num_vars] = size(model_del.S);

model_del.A = sparse(model_del.S);
model_del.varNames = model_del.rxns;
model_del.vtype = repelem('C', num_vars)';
model_del.varlb = model_del.lb; %lb in MILP
model_del.varub = model_del.ub; %ub in MILP
model_del.obj = model_del.c;
model_del.eqNames = model_del.mets;
model_del.eqtype = repelem('=', num_eqs)'; %csense in MILP
model_del.rhs = model_del.b;
model_del.modelsense = 'max';
model_del.x0 = [];
end




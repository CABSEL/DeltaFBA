function [model_del, nochange_idx] = createDeltaModel(model, essential_idx, zeroflow_idx,maxflux_val)
% For the inputted model and index of essential reactions create a delta
% model. The delta model will first convert the original model into an
% irreversible model. In the next step, the bounds of the variables in
% the delta model will be set based on the minflux_val. Finally the bounds
% of the essential reactions will be set to 0 as these reactions will not
% be allowed to change. These are known zero delta reaction fluxes. Also we
% will set zero flow flux reaction bounds to zero.
%
%
% Inputs:
% 1. Model - Genome scale metabolic model with COBRA defined fields
% 2. essential_idx - Index of reactions in the model which do not change
% and are constant across both conditions. If empty - use []
% 3. zeroflow_idx - Index of reactions in the model which do not have any
% flow and are constrained to 0. If empty - use []
% 4. maxflux_val - Ceiling of value that the delta model can be changed to.
%
%
% Outputs:
% 1. model_del - Model with bounds changes and transformed to irreversible 
% 2. nochange_idx - Index of all reactions that do not change.

% Step 1: Convert to irreversible model (based on Cobra's implementation
model_tmp = convertToIrreversible(model);

% Step 2: Look at the minflux_val and convert all lower bounds and upper
% bounds

bound_val = 10^ceil(log10(maxflux_val));

model_tmp.lb = -bound_val*ones(size(model_tmp.lb,1),1);
model_tmp.ub = bound_val*ones(size(model_tmp.ub,1),1);

% Step 3: Convert essential reaction bounds to zero
rev_essential_idx = nonzeros(model_tmp.match(essential_idx));
essential_idx = union(essential_idx, rev_essential_idx);

model_tmp.lb(essential_idx) = 0;
model_tmp.ub(essential_idx) = 0;

rev_zeroflow_idx = nonzeros(model_tmp.match(zeroflow_idx));
zeroflow_idx = union(zeroflow_idx, rev_zeroflow_idx);

model_tmp.lb(zeroflow_idx) = 0;
model_tmp.ub(zeroflow_idx) = 0;

nochange_idx = union(essential_idx, zeroflow_idx);
model_del = model_tmp;

% Step 4: make any objective function coefficents zero
model_del.c(find(model_del.c)) = 0;
end




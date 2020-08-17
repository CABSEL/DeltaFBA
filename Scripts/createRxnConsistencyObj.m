function [MILPstructure, model_binary] = createRxnConsistencyObj(model_Bin,de_rxns,de_indUP, de_indDOWN, nochange_idx,option, weights, weights_option)
% Function for taking a model with binary use variables for
% delta(mutant-control) and adding objective equations for maximizing
% consistency with all upregulated, all downregulated and at the same time
% minimizing inconsistency.
%
% Inputs:
% 1. model_Bin - Model with binary descriptors for every reaction.
% 2. de_rxns - Differentially regulated reactions
% 3. de_indUP - Index of Upregulated reactions
% 4. de_indDOWN - Index of downregulated reactions
% 5. nochange_idx - Leave [] if empty. Excluding reactions that are not bounded
%    0
% 6. option:
%   Option = 1 - exclude no change reactions from Z3
%   Option = 0 - Include no change reactions from Z3
% 7. Weights - weights are for the objective function to prioritize reactions for
% regulation. 
% 8. Weights option:
%    1 = Include weights in model; 
%    0 = do not include weights in model.
%
% Outputs:
% 1. MILPstructure - Gurobi MILP structure that can be used to run the
% optimization problem
% 2. model_bvinary - Complete model that can used as an input to any
% optimizer that can solve the MILP/QP.

% Note that de indices correspond to de_rxns
% Note that the no_change_idx corresponds to the irreversible model reactions

model = model_Bin;
[num_eqs, num_vars] = size(model.A);

% Creating an equation for maximizing Z1 of all upregulated reactions.
% Since the reaction can be upregulated both in forward and the reverse
% directions, we add coefficients to both directions. The new variable will
% be sum(z1_u). sum(z2_u) will define inconsistency

Uprxns = (de_rxns(de_indUP));
Uweights = weights(de_indUP);
[~,upIdx] = ismember(strcat('Z1_',Uprxns),model.varNames);
up_weights = Uweights;

model.A(num_eqs+1, num_vars+1) = 1;
model.eqtype(num_eqs+1) = '=';
model.rhs(num_eqs+1) = 0;
model.eqNames{(num_eqs+1),1} = 'MaxSumZ1';
model.varNames{(num_vars+1),1} = 'SumZ1_U';
model.varlb(num_vars+1) = 0;
model.varub(num_vars+1) = numel(upIdx);
model.vtype(num_vars+1) = 'C';
model.obj(num_vars+1) = 0;

if weights_option==1
    model.A(num_eqs+1,upIdx) = -up_weights;
else
    model.A(num_eqs+1,upIdx) = -1;
end
model.Zindex.Z1 = upIdx;


model.A(num_eqs+2, num_vars+2) = 1;
model.eqtype(num_eqs+2) = '=';
model.rhs(num_eqs+2) = 0;
model.eqNames{(num_eqs+2),1} = 'MinSumZ2';
model.varNames{(num_vars+2),1} = 'SumZ2_U';
model.varlb(num_vars+2) = 0;
model.varub(num_vars+2) = numel(upIdx);
model.vtype(num_vars+2) = 'C';
model.obj(num_vars+2) = 0;

if weights_option==1
    model.A(num_eqs+2,upIdx+1) = -up_weights;
else
    model.A(num_eqs+2,upIdx+1) = -1;
end


% Creating an equation for maximizing Z2 of all downregulated reactions.
% Since the reaction can be downregulated both in forward and the reverse
% directions, we add coefficients to both directions. The new variable will
% be sum(zz_d). sum(z1_d) will define inconsistency


Downrxns = (de_rxns(de_indDOWN));
Dweights = 1./weights(de_indDOWN);
[~,downIdx] = ismember(strcat('Z2_',Downrxns),model.varNames);
down_weights = Dweights;

model.A(num_eqs+3, num_vars+3) = 1;
model.eqtype(num_eqs+3) = '=';
model.rhs(num_eqs+3) = 0;
model.eqNames{(num_eqs+3),1} = 'MaxSumZ2';
model.varNames{(num_vars+3),1} = 'SumZ2_D';
model.varlb(num_vars+3) = 0;
model.varub(num_vars+3) = numel(downIdx)*10000;
model.vtype(num_vars+3) = 'C';
model.obj(num_vars+3) = 0;

if weights_option==1
    model.A(num_eqs+3,downIdx) = -down_weights;
else
    model.A(num_eqs+3,downIdx) = -1;
end
model.Zindex.Z2 = downIdx;


model.A(num_eqs+4, num_vars+4) = 1;
model.eqtype(num_eqs+4) = '=';
model.rhs(num_eqs+4) = 0;
model.eqNames{(num_eqs+4),1} = 'MinSumZ1';
model.varNames{(num_vars+4),1} = 'SumZ1_D';
model.varlb(num_vars+4) = 0;
model.varub(num_vars+4) = numel(downIdx)*10000;
model.vtype(num_vars+4) = 'C';
model.obj(num_vars+4) = 0;

if weights_option==1
    model.A(num_eqs+4,downIdx-1) = -down_weights;
else
    model.A(num_eqs+4,downIdx-1) = -1;
end

%Creating an equation for Z1u and Z2dd - Max
model.A(num_eqs+5, num_vars+5) = 1;
model.eqtype(num_eqs+5) = '=';
model.rhs(num_eqs+5) = 0;
model.eqNames{(num_eqs+5),1} = 'MaxSumZ1_Z2';
model.varNames{(num_vars+5),1} = 'SumZ_DE';
model.varlb(num_vars+5) = 0;
model.varub(num_vars+5) = (numel(downIdx)+numel(upIdx))*10000;
model.vtype(num_vars+5) = 'C';
model.obj(num_vars+5) = 0;

model.A(num_eqs+5,num_vars+1) = -1;
model.A(num_eqs+5,num_vars+3) = -1;


%Creating an equation for Z1 and Z2u - Min
model.A(num_eqs+6, num_vars+6) = 1;
model.eqtype(num_eqs+6) = '=';
model.rhs(num_eqs+6) = 0;
model.eqNames{(num_eqs+6),1} = 'MinSumZ1_Z2';
model.varNames{(num_vars+6),1} = 'SumZ_DE_MIN';
model.varlb(num_vars+6) = 0;
model.varub(num_vars+6) = (numel(downIdx)+numel(upIdx))*10000;
model.vtype(num_vars+6) = 'C';
model.obj(num_vars+6) = 0;

model.A(num_eqs+6,num_vars+2) = -1;
model.A(num_eqs+6,num_vars+4) = -1;

%Creating an equation for Objective
model.A(num_eqs+7, num_vars+7) = 1;
model.eqtype(num_eqs+7) = '=';
model.rhs(num_eqs+7) = 0;
model.eqNames{(num_eqs+7),1} = 'Obj';
model.varNames{(num_vars+7),1} = 'Obj';
model.varlb(num_vars+7) = 0;
model.varub(num_vars+7) = (numel(downIdx)+numel(upIdx))*10000;
model.vtype(num_vars+7) = 'C';
model.obj(num_vars+7) = 1;

model.A(num_eqs+7,num_vars+5) = -1;
model.A(num_eqs+7,num_vars+6) = 1;

num_rev = numel(find(model.rev))/2;
num_irrev = size(model.S,2)-num_rev;
ori_rxns = num_rev+num_irrev;
% Creating an equation for Z3 of all unchanged reactions.

z3idx = (ori_rxns+num_irrev+3:3:numel(model.varNames)-7);
nochange_z3 = ori_rxns+num_irrev+((nochange_idx-1)*3)+3;

if option==1
    z3idx = setdiff(z3idx,nochange_z3);
end

z3idx = setdiff(z3idx,union(upIdx+2,downIdx+1));
model.A(num_eqs+8, num_vars+8) = 1;
model.eqtype(num_eqs+8) = '=';
model.rhs(num_eqs+8) = 0;
model.eqNames{(num_eqs+8),1} = 'MaxSumZ3';
model.varNames{(num_vars+8),1} = 'SumZ3';
model.varlb(num_vars+8) = 0;
model.varub(num_vars+8) = numel(z3idx)*10000;
model.vtype(num_vars+8) = 'C';
model.obj(num_vars+8) = 0;

model.A(num_eqs+8,z3idx) = -1;
model.Zindex.Z3 = z3idx;

% All Z equation
model.A(num_eqs+9, num_vars+9) = 1;
model.eqtype(num_eqs+9) = '=';
model.rhs(num_eqs+9) = 0;
model.eqNames{(num_eqs+9),1} = 'MaxSumZ';
model.varNames{(num_vars+9),1} = 'SumZ';
model.varlb(num_vars+9) = 0;
model.varub(num_vars+9) = (numel(z3idx)+numel(downIdx)+numel(upIdx))*10000;
model.vtype(num_vars+9) = 'C';
model.obj(num_vars+9) = 0;

model.A(num_eqs+9,num_vars+8) = -1;
model.A(num_eqs+9,num_vars+5) = -1;

model_binary = model;

MILPstructure.A = model.A;
MILPstructure.varNames = model.varNames;
MILPstructure.vtype = model.vtype;
MILPstructure.lb = model.varlb;
MILPstructure.ub = model.varub;
MILPstructure.obj = model.obj;
MILPstructure.eqNames = model.eqNames;
MILPstructure.sense = model.eqtype;
MILPstructure.rhs = model.rhs;
MILPstructure.modelsense = model.modelsense;
MILPstructure.x0 = model.x0;
end

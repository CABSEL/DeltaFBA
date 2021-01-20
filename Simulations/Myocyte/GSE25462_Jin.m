%% Clear Workspace
clear all;
clc

%% Set paths, initiate cobra and solver parameters
addpath(genpath('D:\GitHub\D_FBA'));
addpath(genpath('C:\Users\ravis\Desktop\GPR Training\cobratoolbox'));
initCobraToolbox()
[cobraParams,solverParams] = parseSolverParameters('MILP');
cobraParams.timeLimit = 3000;
cobraParams.printLevel = 1;

%% Load metabolic model
load('D:\Github\D_FBA\GEMs\iMyocyte2419.mat');
model = iMyocyte2419;

%% Set maxflux val
maxflux_val = 200;

%% Convert to irreversible model
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

%% Check for essential reactions that are constrained
essential_idx = union(find(model.ub~=0 & model.ub<1000), find(model.lb~=0 & model.lb>-1000));
% Note = Check essential reactions and prune set

%% Check for reactions that have no flux flows
zeroflow_idx = find(model.lb==0 & model.ub==0); %zero fluxes by model definition

%% Create delta model by setting bounds of all reactions
[model_del, nochange_idx] = createDeltaModel(model, [], zeroflow_idx,maxflux_val);

%% Set parameters for Delta FBA
epsilon = 10;
M_prime =1e5;

%% Add NetFlux Variable
[model_delNet] = addNetDeltaFlux(model_del, model.rxns);

%% Incorporation of gene expression data - Read gene expression data
[num,txt,~] = xlsread('D:\GitHub\D_FBA\Data\Myocytes\GSE25462_DE.xlsx');
genes = txt(2:end,9);
genes_data = num(:,5);
genes_data = 2.^genes_data;

expressionData.gene = model.genes;
expressionData.value = zeros(numel(model.genes),1);

for i = 1:numel(model.genes)
    [~,ind]=ismember(model.genes{i},genes);
    if ind>0
        expressionData.value(i)=genes_data(ind);
    end
end

RxnExpr = mapExpressionToReactions(model,expressionData);
RxnExpr(find(RxnExpr==(-1))) = 0;
de_exprrxns = mapExpressionToReactions(modelIrrev,expressionData);
de_rxns = modelIrrev.rxns;

sel = find(~(de_exprrxns==(-1) | de_exprrxns==0));
de_exprrxns = de_exprrxns(sel);
de_rxns = de_rxns(sel);

%% Define up and down regulated reactions
cut_off = 0.05;
[~,de_indUP]=maxk(de_exprrxns,round(cut_off*size(find(de_exprrxns>1),1)));
[~,de_indDOWN]=mink(de_exprrxns,round(cut_off*size(find(de_exprrxns<1),1)));

regRxns=[de_rxns(de_indUP); de_rxns(de_indDOWN)];
regRxnsRatio = [de_exprrxns(de_indUP) ; de_exprrxns(de_indDOWN)];

%% Create binary variables for all reactions
% Choose 1 of the following two statements

% Statement A : Relaxed conditions
  model_Bin = createBinaryUseVariable(model_delNet, epsilon, M_prime, regRxns, ones(numel(regRxns),1));

% Statement B : Stringent conditions (when using this condition, please
% make sure that generatio isnt out of bounds
% regRxnsRatio(find(abs(regRxnsRatio)<0.1)) = 0.1;
% regRxnsRatio(find(abs(regRxnsRatio)>100)) = 100;
% model_Bin = createBinaryUseVariable_leg(model_delNet, epsilon, M_prime, regRxns, regRxnsRatio);


%% Creating objectives for maximizing the consistency of DE reactions
% Choose 1 of the following two statements

% Statement A : Without weighting the objective function with gene
% expression changes
  [MILPstructure, model_binary] = createRxnConsistencyObj(model_Bin,de_rxns,de_indUP, de_indDOWN, nochange_idx, 1, de_exprrxns, 0);

% Statement B : With weighting the objective function with gene
% expression changes
% de_exprrxns(find(abs(de_exprrxns)<0.1)) = 0.1;
% de_exprrxns(find(abs(de_exprrxns)>100)) = 100;
% [MILPstructure, model_binary] = createRxnConsistencyObj(model_Bin,de_rxns,de_indUP, de_indDOWN, nochange_idx, 1, de_exprrxns, 1);

[~,NFind]=ismember(strcat('NF_',model.rxns),MILPstructure.varNames);



%% Solve for 1. Maximizing consistency in DE reactions
%Setup parameters for optimizer in gurobi
params.OutputFlag = 1;
params.DisplayInterval = 5;
params.TimeLimit = cobraParams.timeLimit;

sol_Z_de = gurobi(MILPstructure, params);

%% Solve for 2. Minimizing the 2-norm of dFBA results
repMILP = MILPstructure;
[~,num_vars] = size(repMILP.A);
repMILP.lb(find(repMILP.obj)) = sol_Z_de.x(find(repMILP.obj));
%repMILP.ub(find(repMILP.obj)) = sol_Z_de.x(find(repMILP.obj));

Qmat = sparse(num_vars,num_vars);
for i = 1:numel(NFind)
    Qmat(NFind(i), NFind(i)) = 1;
end

repMILP = rmfield(repMILP,'obj');
repMILP.Q = Qmat;
repMILP.modelsense = 'min';
sol_rep = gurobi(repMILP, params);

[MILP] = create1normSolutionMILP(MILPstructure, NFind, sol_Z_de.objval);
sol_rep2 = gurobi(MILP,params);
P_D2 = sol_rep2.x(NFind);
P_D = sol_rep.x(NFind);
R_D = RxnExpr;

xlswrite('D:\GitHub\D_FBA\Simulations\Myocyte\Results_5_2_GSE25462.xlsx',P_D,1);
xlswrite('D:\GitHub\D_FBA\Simulations\Myocyte\Results_5_2_GSE25462.xlsx',R_D,2);
xlswrite('D:\GitHub\D_FBA\Simulations\Myocyte\Results_5_1_GSE25462.xlsx',P_D2,1);
xlswrite('D:\GitHub\D_FBA\Simulations\Myocyte\Results_5_1_GSE25462.xlsx',R_D,2);
save('D:\GitHub\D_FBA\Simulations\Myocyte\GSE25462_5.mat');

%% Clear Workspace
clear all;
clc

%% Set paths, initiate cobra and solver parameters
addpath(genpath('D:\GitHub\D_FBA'));
addpath(genpath('C:\Users\ravis\Desktop\GPR Training\cobratoolbox'));
initCobraToolbox()
[cobraParams,solverParams] = parseSolverParameters('MILP');
cobraParams.timeLimit = 3000;
cobraParams.printLevel = 1;

%% Load metabolic model
load('D:\Github\D_FBA\GEMs\iMyocyte2419.mat');
model = iMyocyte2419;

%% Set maxflux val
maxflux_val = 200;

%% Convert to irreversible model
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

%% Check for essential reactions that are constrained
essential_idx = union(find(model.ub~=0 & model.ub<1000), find(model.lb~=0 & model.lb>-1000));
% Note = Check essential reactions and prune set

%% Check for reactions that have no flux flows
zeroflow_idx = find(model.lb==0 & model.ub==0); %zero fluxes by model definition

%% Create delta model by setting bounds of all reactions
[model_del, nochange_idx] = createDeltaModel(model, [], zeroflow_idx,maxflux_val);

%% Set parameters for Delta FBA
epsilon = 10;
M_prime =1e5;

%% Add NetFlux Variable
[model_delNet] = addNetDeltaFlux(model_del, model.rxns);

%% Incorporation of gene expression data - Read gene expression data
[num,txt,~] = xlsread('D:\GitHub\D_FBA\Data\Myocytes\GSE25462_DE.xlsx');
genes = txt(2:end,9);
genes_data = num(:,5);
genes_data = 2.^genes_data;

expressionData.gene = model.genes;
expressionData.value = zeros(numel(model.genes),1);

for i = 1:numel(model.genes)
    [~,ind]=ismember(model.genes{i},genes);
    if ind>0
        expressionData.value(i)=genes_data(ind);
    end
end

RxnExpr = mapExpressionToReactions(model,expressionData);
RxnExpr(find(RxnExpr==(-1))) = 0;
de_exprrxns = mapExpressionToReactions(modelIrrev,expressionData);
de_rxns = modelIrrev.rxns;

sel = find(~(de_exprrxns==(-1) | de_exprrxns==0));
de_exprrxns = de_exprrxns(sel);
de_rxns = de_rxns(sel);

%% Define up and down regulated reactions
cut_off = 0.25;
[~,de_indUP]=maxk(de_exprrxns,round(cut_off*size(find(de_exprrxns>1),1)));
[~,de_indDOWN]=mink(de_exprrxns,round(cut_off*size(find(de_exprrxns<1),1)));

regRxns=[de_rxns(de_indUP); de_rxns(de_indDOWN)];
regRxnsRatio = [de_exprrxns(de_indUP) ; de_exprrxns(de_indDOWN)];

%% Create binary variables for all reactions
% Choose 1 of the following two statements

% Statement A : Relaxed conditions
  model_Bin = createBinaryUseVariable(model_delNet, epsilon, M_prime, regRxns, ones(numel(regRxns),1));

% Statement B : Stringent conditions (when using this condition, please
% make sure that generatio isnt out of bounds
% regRxnsRatio(find(abs(regRxnsRatio)<0.1)) = 0.1;
% regRxnsRatio(find(abs(regRxnsRatio)>100)) = 100;
% model_Bin = createBinaryUseVariable_leg(model_delNet, epsilon, M_prime, regRxns, regRxnsRatio);


%% Creating objectives for maximizing the consistency of DE reactions
% Choose 1 of the following two statements

% Statement A : Without weighting the objective function with gene
% expression changes
  [MILPstructure, model_binary] = createRxnConsistencyObj(model_Bin,de_rxns,de_indUP, de_indDOWN, nochange_idx, 1, de_exprrxns, 0);

% Statement B : With weighting the objective function with gene
% expression changes
% de_exprrxns(find(abs(de_exprrxns)<0.1)) = 0.1;
% de_exprrxns(find(abs(de_exprrxns)>100)) = 100;
% [MILPstructure, model_binary] = createRxnConsistencyObj(model_Bin,de_rxns,de_indUP, de_indDOWN, nochange_idx, 1, de_exprrxns, 1);

[~,NFind]=ismember(strcat('NF_',model.rxns),MILPstructure.varNames);



%% Solve for 1. Maximizing consistency in DE reactions
%Setup parameters for optimizer in gurobi
params.OutputFlag = 1;
params.DisplayInterval = 5;
params.TimeLimit = cobraParams.timeLimit;

sol_Z_de = gurobi(MILPstructure, params);

%% Solve for 2. Minimizing the 2-norm of dFBA results
repMILP = MILPstructure;
[~,num_vars] = size(repMILP.A);
repMILP.lb(find(repMILP.obj)) = sol_Z_de.x(find(repMILP.obj));
%repMILP.ub(find(repMILP.obj)) = sol_Z_de.x(find(repMILP.obj));

Qmat = sparse(num_vars,num_vars);
for i = 1:numel(NFind)
    Qmat(NFind(i), NFind(i)) = 1;
end

repMILP = rmfield(repMILP,'obj');
repMILP.Q = Qmat;
repMILP.modelsense = 'min';
sol_rep = gurobi(repMILP, params);

[MILP] = create1normSolutionMILP(MILPstructure, NFind, sol_Z_de.objval);
sol_rep2 = gurobi(MILP,params);
P_D2 = sol_rep2.x(NFind);
P_D = sol_rep.x(NFind);
R_D = RxnExpr;

xlswrite('D:\GitHub\D_FBA\Simulations\Myocyte\Results_25_2_GSE25462.xlsx',P_D,1);
xlswrite('D:\GitHub\D_FBA\Simulations\Myocyte\Results_25_2_GSE25462.xlsx',R_D,2);
xlswrite('D:\GitHub\D_FBA\Simulations\Myocyte\Results_25_1_GSE25462.xlsx',P_D2,1);
xlswrite('D:\GitHub\D_FBA\Simulations\Myocyte\Results_25_1_GSE25462.xlsx',R_D,2);
save('D:\GitHub\D_FBA\Simulations\Myocyte\GSE25462_25.mat');
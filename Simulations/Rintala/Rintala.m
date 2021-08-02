%% Clear Workspace
clear all;
clc
root = 'D:\DeltaFBA';
cobrapath = 'C:\Users\ravis\Desktop\GPR Training\cobratoolbox';

%% Set paths, initiate cobra and solver parameters
addpath(genpath(root));
addpath(genpath(cobrapath));
initCobraToolbox()
[cobraParams,solverParams] = parseSolverParameters('MILP');
cobraParams.timeLimit = 1000;
cobraParams.printLevel = 1;

%% Load metabolic model
load('GEMs\iMM904.mat');

%% Read model and optimze Wild type model
model = iMM904;
model = generateRules(model);
FBA = optimizeCbModel(model);

%% Find minimized model flux results to Set maximum flux of a single reaction

[MinimizedFlux modelIrrev]= minimizeModelFlux(model);
if isfield(FBA, 'obj')==1
    obj = FBA.obj;
else isfield(FBA,'f')== 1
    obj = FBA.f;
end

modelIrrev.lb(find(model.c)) = obj;
modelIrrev.ub(find(model.c)) = obj;
minflux = optimizeCbModel(modelIrrev,'min');

if isfield(minflux, 'full')==1
    vec = minflux.full;
else isfield(minflux,'x')== 1
    vec = minflux.x;
end

maxflux_val = max(abs(vec(1:(size(model.rxns,1)+numel(find(model.rev))))));

%% Convert to irreversible model and check results (Sanity Check)
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
FBA_ir = optimizeCbModel(modelIrrev);

%% Check for essential reactions that are constrained
essential_idx = union(find(model.ub~=0 & model.ub<999999), find(model.lb~=0 & model.lb>-999999));
% Note = Check essential reactions and prune set

%% Check for reactions that have no flux flows
zeroflow_idx = find(model.lb==0 & model.ub==0); %zero fluxes by model definition

%% Create delta model by setting bounds of all reactions
[model_del, nochange_idx] = createDeltaModel(model, [], zeroflow_idx,maxflux_val);

%% Set parameters
epsilon = 1;
M_prime =1e6;

%% Add NetFlux Variable
[model_delNet] = addNetDeltaFlux(model_del, model.rxns);

%% Load flux data
[num,txt,~] = xlsread('Data\Rintala\Fluxomics.xlsx');
conditions = num(1,1:end)*100;
conditions = num2cell(conditions);
flux_data = num(2:end,:)./1000;
flux_list = txt(2:end,1);
checkFlux = (numel(find(ismember(flux_list, model.rxns)))==numel(flux_list));

%% Growth rate
biomass = [5.24 4.735 2.96 2.145 1.05];
%biomass = (biomass./5.24).*obj;
%scale = biomass./biomass(1);
o2 = [-2.7 -2.5 -1.7 -1.2 -0.0];

%o2scale = o2.*scale;
%o2scale = o2.*scale;
%o2b = o2.*biomass;
%o2b = (o2b./o2b(1)).*FBA.x(essential_idx(3));

%% Preinitialize Results vectors
M_D = zeros(numel(flux_list),numel(conditions)-1);
P_D = zeros(numel(flux_list),numel(conditions)-1);
R_D = zeros(numel(flux_list),numel(conditions)-1);
fullP = zeros(numel(model.rxns),numel(conditions)-1);
Sol_stat = zeros(numel(conditions)-1,1); %1 is optimization for l2 norm is successful
Sol_stat2 = zeros(numel(conditions)-1,1); % MIPGAP for l2 norm

%% Identify gene expression for all metabolic genes in the condition
[num,txt,~] = xlsread('Data\Rintala\Transcriptomics.xlsx');
genes=txt(1:end,1);
genes_data = num(2:end,[5 4 3 2 1]);


for j = 1:numel(conditions)-1
    
    wt=1;
    mut=j+1;
    geneWT=nan(numel(model.genes),1);
    geneMUT=nan(numel(model.genes),1);
    for i=1:numel(model.genes)
        [~,ind]=ismember(model.genes{i},genes);
        if ind>0
            geneWT(i)=genes_data(ind,wt);
            geneMUT(i)=genes_data(ind,mut);
        end
    end
    
    %% Calculate ratio and filter
    Generatio=geneMUT./geneWT;
    expr_genes = model.genes;
    sel = find(~isnan(Generatio));
    Generatio = Generatio(sel);
    expr_genes = expr_genes(sel);
    expressionData.gene = expr_genes;
    expressionData.value = Generatio;
    
    RxnExpr = mapExpressionToReactions(model,expressionData);
    [~,ids] = ismember(flux_list, model.rxns);
    RxnExpr_mod = RxnExpr(ids);
    RxnExpr_mod(find(RxnExpr_mod==(-1))) = 0;
    
    de_exprrxns = mapExpressionToReactions(modelIrrev,expressionData);
    de_rxns = modelIrrev.rxns;
    sel = find(~(de_exprrxns==(-1) | de_exprrxns==0));
    de_exprrxns = de_exprrxns(sel);
    de_rxns = de_rxns(sel);
    
    %% Define up and down regulated reactions
    cut_off = 0.10;
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
    % regRxnsRatio(find(abs(regRxnsRatio)<0.05)) = 0.05;
    % regRxnsRatio(find(abs(regRxnsRatio)>100)) = 100;
    % model_Bin = createBinaryUseVariable(model_delNet, epsilon, M_prime, regRxns, regRxnsRatio);
    
    
    %% Creating objectives for maximizing the consistency of DE reactions
    % Choose 1 of the following two statements
    
    % Statement A : Without weighting the objective function with gene
    % expression changes
    [MILPstructure, model_binary] = createRxnConsistencyObj(model_Bin,de_rxns,de_indUP, de_indDOWN, nochange_idx, 1, de_exprrxns, 0);
    
    % Statement B : With weighting the objective function with gene
    % expression changes
    % de_exprrxns(find(abs(de_exprrxns)<0.05)) = 0.05;
    % de_exprrxns(find(abs(de_exprrxns)>100)) = 100;
    % [MILPstructure, model_binary] = createRxnConsistencyObj(model_Bin,de_rxns,de_indUP, de_indDOWN, nochange_idx, 1, de_exprrxns, 1);
    
    [~,NFind]=ismember(strcat('NF_',model.rxns),MILPstructure.varNames);
    
    
    %% Add additional constraints - e.g. exchange fluxes
    % Grwoth rate - Biomass 7-8
    %MILPstructure.lb(NFind(find(model.c))) = biomass(mut) - biomass(1);
    %MILPstructure.ub(NFind(find(model.c))) = biomass(mut) - biomass(1);
    
    MILPstructure.lb(NFind(essential_idx(3))) = o2(mut) - o2(1);
    MILPstructure.ub(NFind(essential_idx(3))) = o2(mut) - o2(1);
    
    % All exchange measurements
    e_ids = find(~cellfun(@isempty,regexp(flux_list,'EX_')));
    [~,ids] = ismember(strcat('NF_',flux_list(e_ids)), MILPstructure.varNames);
    MILPstructure.lb(ids) = (flux_data(e_ids,mut)) - (flux_data(e_ids,1));
    MILPstructure.ub(ids) = (flux_data(e_ids,mut)) - (flux_data(e_ids,1));
    
    %% Solve for 1. Maximizing consistency in DE reactions
    %Setup parameters for optimizer in gurobi
    params.OutputFlag = 1;
    params.DisplayInterval = 5;
    params.TimeLimit = cobraParams.timeLimit;
    %params.MIPGap = 1e-3;
    params.IntFeasTol = 1e-09;
    params.FeasibilityTol = cobraParams.feasTol;
    params.OptimalityTol = cobraParams.optTol;
    
    sol_Z_de = gurobi(MILPstructure, params);
    
    %% Solve for 2. Minimizing the 1-norm of dFBA results
    
    repMILP = MILPstructure;
    [~,num_vars] = size(repMILP.A);
    repMILP.lb(find(repMILP.obj)) = sol_Z_de.x(find(repMILP.obj));
    repMILP.ub(find(repMILP.obj)) = sol_Z_de.x(find(repMILP.obj));
    
    Qmat = sparse(num_vars,num_vars);
    for i = 1:numel(NFind)
        Qmat(NFind(i), NFind(i)) = 1;
    end
    
    repMILP = rmfield(repMILP,'obj');
    repMILP.Q = Qmat;
    repMILP.modelsense = 'min';
    sol_rep = gurobi(repMILP, params);
    if isfield(sol_rep,'x')
        sol_rep = sol_rep;
        Sol_stat(j) = 1;
        Sol_stat2(j) = sol_rep.mipgap;
    else
        sol_rep = sol_Z_de;
    end
    
    
    %% Compare Solution with measured flux data
    Measured = (flux_data(:,mut)) - (flux_data(:,1));
    [~,ids] = ismember(strcat('NF_',flux_list), repMILP.varNames);
    Predicted = sol_rep.x(ids);
    sel = setdiff(1:numel(Measured),union(find(Predicted==0), find(Measured==0)));
    Res = CalcPerf(Measured, Predicted);
    
    M_D(:,j) = Measured;
    P_D(:,j) = Predicted;
    R_D(:,j) = RxnExpr_mod;
    fullP(:,j) = sol_rep.x(NFind);
    save(strcat(root,'\Simulations\Rintala\MAT_Files\', num2str(conditions{1,mut}),'.mat'))
    
end

xlswrite(strcat(root,'\Simulations\Rintala\Results.xlsx'),M_D,1);
xlswrite(strcat(root,'\Simulations\Rintala\Results.xlsx'),P_D,2);
xlswrite(strcat(root,'\Simulations\Rintala\Results.xlsx'),R_D,3);
xlswrite(strcat(root,'\Simulations\Rintala\Results.xlsx'),Sol_stat,4);
xlswrite(strcat(root,'\Simulations\Rintala\Results.xlsx'),Sol_stat2,5);
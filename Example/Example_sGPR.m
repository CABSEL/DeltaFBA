% Example to run D_FBA-sGPR on a metabolic model. Here in this example, we show
% a simple example of a E_coli_core model having a dilution change and
% transcriptomic data

%% Clear Workspace and setup Gurobi and COBRA tools
clear all;
clc
cd /usr/local/gurobi811/linux64/matlab
gurobi_setup()
cd /morpheus1/users/ravis
addpath(genpath('/morpheus1/users/ravis/cobratoolbox'));
addpath(genpath('/morpheus1/users/ravis/D_FBA'));
load('/morpheus1/users/ravis/D_FBA/GEMs/e_coli_core.mat');
initCobraToolbox()
[cobraParams,solverParams] = parseSolverParameters('MILP');
cobraParams.timeLimit = 600;
cobraParams.printLevel = 1;

%% Read model and optimze Wild type model
model = e_coli_core;
model = generateRules(model);
model = removeUnusedGenes(model);
FBA = optimizeCbModel(model);

%% Find minimized model flux results
[MinimizedFlux modelIrrev]= minimizeModelFlux(model);
modelIrrev.lb(find(model.c)) = FBA.obj;
modelIrrev.ub(find(model.c)) = FBA.obj;
minflux = optimizeCbModel(modelIrrev,'min');
maxflux_val = max(abs(minflux.full(1:(size(model.rxns,1)+numel(find(model.rev))))));

%% Convert to irreversible model and check results
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
FBA_ir = optimizeCbModel(modelIrrev);

%% Check for essential reactions that are constrained
essential_idx = union(find(model.ub~=0 & model.ub<1000), find(model.lb~=0 & model.lb>-1000));
% Note = Check essential reactions and prune set

%% Check for reactions that have no flux flows
zeroflow_idx = find(model.lb==0 & model.ub==0); %zero fluxes by model definition

%% Create delta model by setting bounds of all reactions
[model_del, nochange_idx] = createDeltaModel_rev(model, [], zeroflow_idx,maxflux_val);
IrGPRs = GPR2cell(model);
m_GPRT = createGPRTmodel(model_del, IrGPRs);

%% Set parameters
epsilon = 0.1;
M_prime =1e5;

%% Incorporation of gene expression data - Read gene expression data
[num,txt,~] = xlsread('/morpheus1/users/ravis/D_FBA/Data/Ishii/Transcriptomics_RTPCR.xlsx');
genes=txt(2:end,1);
conditions = txt(1,2:end);
genes_data=num;

%% Load flux data
[num,txt,~] = xlsread('/morpheus1/users/ravis/D_FBA/Data/Ishii/Fluxomics.xlsx');
checkCond = isequal(conditions, txt(1,2:end));
flux_data = num(2:end,:);
flux_list = txt(3:end,1);
checkFlux = (numel(find(ismember(flux_list, model.rxns)))==numel(flux_list));

%% Load geneName to RegulonDB for gene deletions
[~,g2r,~] = xlsread('/morpheus1/users/ravis/D_FBA/Data/Ishii/Ecoli_genename_2_regulon.xlsx');


%% Run D_FBA for a pair of conditions
wt=1;
mut=2;
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

nsets = zeros(numel(expr_genes),1);
for ns = 1:numel(expr_genes)
    allrx = struct2cell(findRxnsFromGenes(model, expr_genes(ns)));
    allrx = vertcat(allrx{:});
    nsets(ns) = numel(allrx(:,1));
end

Generatio(Generatio>1) = Generatio(Generatio>1).*nsets(Generatio>1);
Generatio(Generatio<1) = Generatio(Generatio<1)./nsets(Generatio<1);

%% Define up and down regulated reactions
cut_off = 1;
[~,de_indUP]=maxk(Generatio,round(cut_off*size(find(Generatio>(1)),1)));
[~,de_indDOWN]=mink(Generatio,round(cut_off*size(find(Generatio<(1)),1)));

regGenes=[expr_genes(de_indUP); expr_genes(de_indDOWN)];
regGeneRatio = [Generatio(de_indUP) ; Generatio(de_indDOWN)];

regRxnsUp = struct2cell(findRxnsFromGenes(model, expr_genes(de_indUP)));
regRxnsUp = vertcat(regRxnsUp{:});
regRxnsDown = struct2cell(findRxnsFromGenes(model, expr_genes(de_indDOWN)));
regRxnsDown = vertcat(regRxnsDown{:});

%% Create binary variables for all reactions
% Choose 1 of the following two statements

% Statement A : Relaxed conditions
% model_Bin = createBinaryUseVariable_enzyme(m_GPRT, epsilon, M_prime, regGenes, ones(numel(regGenes),1));

% Statement B : Stringent conditions (when using this condition, please
% make sure that generatio isnt out of bounds
  regGeneRatio(find(abs(regGeneRatio)<0.1)) = 0.1;
  regGeneRatio(find(abs(regGeneRatio)>10)) = 10;
  model_Bin = createBinaryUseVariable_enzyme(m_GPRT, epsilon, M_prime, regGenes, regGeneRatio);

%% Creating objectives for maximizing the consistency of DE reactions
% Choose 1 of the following two statements

% Statement A : Without weighting the objective function with gene
% expression changes
  [MILPstructure, model_binary] = createGeneConsistencyObj(model_Bin, expr_genes, de_indUP, de_indDOWN, Generatio, 0);

% Statement B : With weighting the objective function with gene
% expression changes
% Generatio(find(abs(Generatio)<0.1)) = 0.1;
% Generatio(find(abs(Generatio)>10)) = 10;
% [MILPstructure, model_binary] = createGeneConsistencyObj(model_Bin, expr_genes, de_indUP, de_indDOWN, Generatio, 1);

[~,NFind]=ismember(model.rxns,MILPstructure.varNames);

%% Add additional constraints - e.g. exchange fluxes
% Glucose exchange measurements - Include more data if desired
e_ids = find(~cellfun(@isempty,regexp(flux_list,'EX_')));
e_ids = e_ids(1);
[~,ids] = ismember(flux_list(e_ids), MILPstructure.varNames);
MILPstructure.lb(ids) = flux_data(e_ids,mut) - flux_data(e_ids,1);
MILPstructure.ub(ids) = flux_data(e_ids,mut) - flux_data(e_ids,1);

%% Solve for 1. Maximizing consistency in DE reactions
%Setup parameters for optimizer in gurobi
params.OutputFlag = 1;
params.DisplayInterval = 5;
params.TimeLimit = cobraParams.timeLimit;
params.IntFeasTol = 1e-09;
params.FeasibilityTol = cobraParams.feasTol;
params.OptimalityTol = cobraParams.optTol;

% Optimize for maximizing consistency
sol_Z_de = gurobi(MILPstructure, params);

%% Solve for 2. Minimizing the L2-norm of dFBA results
repMILP = MILPstructure;
[~,num_vars] = size(repMILP.A);
repMILP.lb(find(repMILP.obj)) = sol_Z_de.x(find(repMILP.obj));
repMILP.ub(find(repMILP.obj)) = sol_Z_de.x(find(repMILP.obj));

u_ind = [];
for b=1:numel(model.genes)
    k1 = find(~cellfun(@isempty,regexp(MILPstructure.varNames,strcat('u_',model.genes{b,1},'_\d'))));
    u_ind = union(k1,u_ind);
end

Qmat = sparse(num_vars,num_vars);
for i = 1:numel(u_ind)
    Qmat(u_ind(i), u_ind(i)) = 1;
end

repMILP = rmfield(repMILP,'obj');
repMILP.Q = Qmat;
repMILP.modelsense = 'min';

% Optimize for minimized L2norm with maximal consistency
sol_rep = gurobi(repMILP, params);


%% Compare Solution with measured flux data
Measured = flux_data(:,mut) - flux_data(:,wt);
[~,ids] = ismember(flux_list, repMILP.varNames);
Predicted = sol_rep.x(ids);
sel = setdiff(1:numel(Measured),union(find(Predicted==0), find(Measured==0)));
Cor = dot(Predicted(sel),Measured(sel))/(norm(Predicted(sel))*norm(Measured(sel)));
Res = CalcPerf(Measured, Predicted);
NRMSE = Res.RMSE;
sel2 = setdiff(1:numel(flux_list),e_ids);
sel = intersect(sel,sel2);
intCor = dot(Predicted(sel),Measured(sel))/(norm(Predicted(sel))*norm(Measured(sel)));
downRecall = numel(intersect(find(Predicted<0), find(Measured<0)))/numel(find(Measured<0));
upRecall = numel(intersect(find(Predicted>0), find(Measured>0)))/numel(find(Measured>0));
[~,u_ind2] = ismember(strcat('u_', regGenes),MILPstructure.varNames);
U_Predicted = sol_rep.x(u_ind2);
M_D = Measured;
P_D = Predicted;
P_D_F = sol_rep.x(NFind);
R_D = regGeneRatio;

%% Escher plots
T1 = table(flux_list, Measured);
T1.Properties.VariableNames = {'ID', 'Flux'};
writetable(T1,strcat('/morpheus1/users/ravis/D_FBA/Example/',conditions{1,mut},'_M_sGPR.csv'),'Delimiter',',')

T2 = table(flux_list, Predicted);
T2.Properties.VariableNames = {'ID', 'Flux'};
writetable(T2,strcat('/morpheus1/users/ravis/D_FBA/Example/',conditions{1,mut},'_P_sGPR.csv'),'Delimiter',',')

T4 = table(model.rxns, P_D_F);
T4.Properties.VariableNames = {'ID', 'Flux'};
writetable(T4,strcat('/morpheus1/users/ravis/D_FBA/Example/',conditions{1,mut},'_PF_sGPR.csv'),'Delimiter',',')

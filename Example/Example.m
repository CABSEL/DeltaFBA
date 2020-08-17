% Example to run D_FBA on a metabolic model. Here in this example, we show
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
[model_del, nochange_idx] = createDeltaModel(model, [], [],maxflux_val);

%% Set parameters
epsilon = 0.1;
M_prime =1e5;

%% Add NetFlux Variable
[model_delNet] = addNetDeltaFlux(model_del, model.rxns);


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

%% Calculate ratio and filter - Convert to reaction expression
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
cut_off = 1;
[~,de_indUP]=maxk(de_exprrxns,round(cut_off*size(find(de_exprrxns>1),1)));
[~,de_indDOWN]=mink(de_exprrxns,round(cut_off*size(find(de_exprrxns<1),1)));

regRxns=[de_rxns(de_indUP); de_rxns(de_indDOWN)];
regRxnsRatio = [de_exprrxns(de_indUP) ; de_exprrxns(de_indDOWN)];

%% Create binary variables that define use for all reactions
% Choose 1 of the following two statements

% Statement A : Relaxed conditions
% model_Bin = createBinaryUseVariable_leg(model_delNet, epsilon, M_prime, regRxns, ones(numel(regRxns),1));

% Statement B : Stringent conditions (when using this condition, please
% make sure that generatio isnt out of bounds
 regRxnsRatio(find(abs(regRxnsRatio)<0.1)) = 0.1;
 regRxnsRatio(find(abs(regRxnsRatio)>10)) = 10;
 model_Bin = createBinaryUseVariable_leg(model_delNet, epsilon, M_prime, regRxns, regRxnsRatio);

%% Gene Deletions - if necessary and present
ids = find(ismember(g2r(:,1), conditions(mut)));
if numel(ids)>0
    ids = find(modelIrrev.rxnGeneMat(:,find(ismember(model.genes,g2r{ids,2}))));
    if numel(ids)>0
        ids = ids(1);
        ids = find(ismember(model_Bin.varNames,strcat('Z2_', modelIrrev.rxns(ids))));
        model_Bin.varlb(ids) = 1;
    end
end

%% Creating objectives for maximizing the consistency with DE reactions
% Choose 1 of the following two statements

% Statement A : Without weighting the objective function with gene
% expression changes
[MILPstructure, model_binary] = createRxnConsistencyObj(model_Bin,de_rxns,de_indUP, de_indDOWN, nochange_idx, 1, de_exprrxns, 0);

% Statement B : With weighting the objective function with gene
% expression changes
% de_exprrxns(find(abs(de_exprrxns)<0.1)) = 0.1;
% de_exprrxns(find(abs(de_exprrxns)>10)) = 10;
% [MILPstructure, model_binary] = createRxnConsistencyObj(model_Bin,de_rxns,de_indUP, de_indDOWN, nochange_idx, 1, de_exprrxns, 1);

[~,NFind]=ismember(strcat('NF_',model.rxns),MILPstructure.varNames);

%% Add additional constraints - e.g. exchange fluxes
% Glucose exchange measurements - Include more data if desired
e_ids = find(~cellfun(@isempty,regexp(flux_list,'EX_')));
e_ids = e_ids(1);
[~,ids] = ismember(strcat('NF_',flux_list(e_ids)), MILPstructure.varNames);
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

Qmat = sparse(num_vars,num_vars);
for i = 1:numel(NFind)
    Qmat(NFind(i), NFind(i)) = 1;
end

repMILP = rmfield(repMILP,'obj');
repMILP.Q = Qmat;
repMILP.modelsense = 'min';

% Optimize for minimized L2norm with maximal consistency
sol_rep = gurobi(repMILP, params);


%% Compare Solution with measured flux data
Measured = flux_data(:,mut) - flux_data(:,wt);
[~,ids] = ismember(strcat('NF_',flux_list), repMILP.varNames);
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
M_D = Measured;
P_D = Predicted;
R_D = RxnExpr_mod;

%% Plots
figure('Units', 'pixels', 'Position', [100 100 800 675]);
s1 = plot(Measured, RxnExpr_mod);
s2 = xline(0,'--');
s3 = yline(1,'--');
set(s2,'Color',[0 0 .5]);
set(s3,'Color', [0 0 .5]);
set(s1,'Linestyle','none','Marker', 'o','MarkerSize', 6, 'MarkerEdgeColor' , 'none', 'MarkerFaceColor' , [.75 .75 1] );
hTitle  = title (strcat('Comparing Measured Flux and Expression Changes -',{' '}, conditions{1,mut}, ' vs Ref'));
hXLabel = xlabel('Flux Changes (mmol/g.DW.h^-^1)');
hYLabel = ylabel('Reaction Expression Ratio');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 12          );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1         );
set(gcf, 'PaperPositionMode', 'auto');

figure('Units', 'pixels', 'Position', [100 100 800 675]);
s1 = plot(Predicted, RxnExpr_mod);
s2 = xline(0,'--');
s3 = yline(1,'--');
set(s2,'Color',[0 0 .5]);
set(s3,'Color', [0 0 .5]);
set(s1,'Linestyle','none','Marker', 'o','MarkerSize', 6, 'MarkerEdgeColor' , 'none', 'MarkerFaceColor' , [.75 .75 1] );
hTitle  = title (strcat('Comparing Predicted Flux and Expression Changes -',{' '}, conditions{1,mut}, ' vs Ref'));
hXLabel = xlabel('Flux Changes (mmol/g.DW.h^-^1)');
hYLabel = ylabel('Reaction Expression Ratio');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 12          );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1         );
set(gcf, 'PaperPositionMode', 'auto');


figure('Units', 'pixels', 'Position', [100 100 800 675]);
s1 = plot(Measured, Predicted);
s2 = xline(0,'--');
s3 = yline(0,'--');
set(s2,'Color',[0 0 .5]);
set(s3,'Color', [0 0 .5]);
set(s1,'Linestyle','none','Marker', 'o','MarkerSize', 6, 'MarkerEdgeColor' , 'none', 'MarkerFaceColor' , [.75 .75 1] );
hTitle  = title (strcat('Comparing Flux Changes -',{' '}, conditions{1,mut}, ' vs Ref'));
hXLabel = xlabel('Measured Flux Changes (mmol/g.DW.h^-^1)');
hYLabel = ylabel('Predicted Flux Changes (mmol/g.DW.h^-^1)');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 12          );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1         );
set(gcf, 'PaperPositionMode', 'auto');


%% Escher plots
T1 = table(flux_list, Measured);
T1.Properties.VariableNames = {'ID', 'Flux'};
writetable(T1,strcat('/morpheus1/users/ravis/D_FBA/Example/',conditions{1,mut},'_M.csv'),'Delimiter',',')

T2 = table(flux_list, Predicted);
T2.Properties.VariableNames = {'ID', 'Flux'};
writetable(T2,strcat('/morpheus1/users/ravis/D_FBA/Example/',conditions{1,mut},'_P.csv'),'Delimiter',',')

T3 = table(flux_list, RxnExpr_mod);
T3.Properties.VariableNames = {'ID', 'Flux'};
writetable(T3,strcat('/morpheus1/users/ravis/D_FBA/Example/',conditions{1,mut},'_R.csv'),'Delimiter',',')

T4 = table(model.rxns, sol_rep.x(NFind));
T4.Properties.VariableNames = {'ID', 'Flux'};
writetable(T4,strcat('/morpheus1/users/ravis/D_FBA/Example/',conditions{1,mut},'_PF.csv'),'Delimiter',',')

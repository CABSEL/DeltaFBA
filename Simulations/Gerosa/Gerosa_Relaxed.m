%% Clear Workspace
clear all;
clc

%% Set paths, initiate cobra and solver parameters
addpath(genpath('D:\GitHub\D_FBA'));
addpath(genpath('C:\Users\ravis\Desktop\GPR Training\cobratoolbox'));
initCobraToolbox()
[cobraParams,solverParams] = parseSolverParameters('MILP');
cobraParams.timeLimit = 600;
cobraParams.printLevel = 1;

%% Load metabolic model
load('D:\Github\D_FBA\GEMs\iJO1366.mat');

%% Read model and optimze Wild type model
model = iJO1366;
model.csense = model.csense';
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
essential_idx = union(find(model.ub~=0 & model.ub<1000), find(model.lb~=0 & model.lb>-1000));
% Note = Check essential reactions and prune set

%% Check for reactions that have no flux flows
zeroflow_idx = find(model.lb==0 & model.ub==0); %zero fluxes by model definition

%% Create delta model by setting bounds of all reactions
[model_del, nochange_idx] = createDeltaModel(model, [], zeroflow_idx,maxflux_val);

%% Set parameters for Delta FBA
epsilon = 0.1;
M_prime =1e5;

%% Add NetFlux Variable
[model_delNet] = addNetDeltaFlux(model_del, model.rxns);

%% Load flux data (40 Measured flux reactions)
[num,txt,~] = xlsread('D:\GitHub\D_FBA\Data\Gerosa\Fluxomics.xlsx');
flux_data = num;
flux_list = txt(2:end,1);
checkFlux = (numel(find(ismember(flux_list, model.rxns)))==numel(flux_list));
conditions = txt(1,2:end);

%% Growth rate data
[gr_rate] = [0.2900, 0.4924, 0.1800, 0.6500, 0.4900, 0.5900, 0.3900, 0.5100];

%% Preinitialize Results vectors for solution stat
sol_stat = zeros(numel(conditions), numel(conditions));
sol_stat2 = zeros(numel(conditions), numel(conditions));

for w = 1:numel(conditions)
    
    %% Preinitialize Results vectors
    Cor = zeros(numel(conditions),1);
    intCor = zeros(numel(conditions),1);
    NRMSE = zeros(numel(conditions),1);
    downRecall = zeros(numel(conditions),1);
    upRecall = zeros(numel(conditions),1);
    M_D = zeros(numel(flux_list),numel(conditions));
    P_D = zeros(numel(flux_list),numel(conditions));
    R_D = zeros(numel(flux_list),numel(conditions));
    
    for j = 1:numel(conditions)
        if w~=j
            %% Identify gene expression for all metabolic genes in the condition
            wt=w;
            mut=j;
            [num] = readtable(strcat('D:\GitHub\D_FBA\Data\Gerosa\RawData\',conditions{1,mut},'-',conditions{1,wt},'.xlsx'));
            genes = table2cell(num(1:end,2));
            genes_data = table2array(num(1:end,7));
            
            geneLFC=nan(numel(model.genes),1);
            for i=1:numel(model.genes)
                [~,ind]=ismember(model.genes{i},genes);
                if ind>0
                    geneLFC(i)=genes_data(ind);
                end
            end
            
            %% Calculate ratio and filter
            Generatio=geneLFC;
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
            cut_off = .05;
            [~,de_indUP]=maxk(de_exprrxns,round(cut_off*size(find(de_exprrxns>0),1)));
            [~,de_indDOWN]=mink(de_exprrxns,round(cut_off*size(find(de_exprrxns<0),1)));
            
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
            % model_Bin = createBinaryUseVariable(model_delNet, epsilon, M_prime, regRxns, regRxnsRatio);
            
            
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
            
            %% Add additional constraints - e.g. exchange fluxes
            % Grwoth rate - Biomass 7-8
            
            MILPstructure.lb(NFind(8)) = gr_rate(mut) - gr_rate(wt);
            MILPstructure.ub(NFind(8)) = gr_rate(mut) - gr_rate(wt);
            
            % All exchange measurements
            e_ids = find(~cellfun(@isempty,regexp(flux_list,'EX_')));
            [~,ids] = ismember(strcat('NF_',flux_list(e_ids)), MILPstructure.varNames);
            MILPstructure.lb(ids) = flux_data(e_ids,mut) - flux_data(e_ids,wt);
            MILPstructure.ub(ids) = flux_data(e_ids,mut) - flux_data(e_ids,wt);
            
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
                sol_stat(wt, mut) = 1;
                sol_stat2(wt,mut) = sol_rep.mipgap;
            else
                sol_rep = sol_Z_de;
            end
            
            %% Compare Solution with measured flux data
            Measured = flux_data(:,mut) - flux_data(:,wt);
            [~,ids] = ismember(strcat('NF_',flux_list), repMILP.varNames);
            Predicted = sol_rep.x(ids);
            sel = setdiff(1:numel(Measured),union(find(Predicted==0), find(Measured==0)));
            Cor(j) = dot(Predicted(sel),Measured(sel))/(norm(Predicted(sel))*norm(Measured(sel)));
            Res = CalcPerf(Measured, Predicted);
            NRMSE(j) = Res.RMSE;
            sel2 = setdiff(1:numel(flux_list),e_ids);
            sel = intersect(sel,sel2);
            intCor(j) = dot(Predicted(sel),Measured(sel))/(norm(Predicted(sel))*norm(Measured(sel)));
            downRecall(j) = numel(intersect(find(Predicted<0), find(Measured<0)))/numel(find(Measured<0));
            upRecall(j) = numel(intersect(find(Predicted>0), find(Measured>0)))/numel(find(Measured>0));
            M_D(:,j) = Measured;
            P_D(:,j) = Predicted;
            R_D(:,j) = RxnExpr_mod;
            save(strcat('D:\GitHub\D_FBA\Simulations\Gerosa\MAT_Files\',conditions{1,j},'vs',conditions{1,wt},'.mat'))

        end
    end
    
    xlswrite(strcat('D:\GitHub\D_FBA\Simulations\Gerosa\Results',num2str(w),'.xlsx'),M_D,1);
    xlswrite(strcat('D:\GitHub\D_FBA\Simulations\Gerosa\Results',num2str(w),'.xlsx'),P_D,2);
    xlswrite(strcat('D:\GitHub\D_FBA\Simulations\Gerosa\Results',num2str(w),'.xlsx'),R_D,3);
    xlswrite(strcat('D:\GitHub\D_FBA\Simulations\Gerosa\Results',num2str(w),'.xlsx'),Cor,4);
    xlswrite(strcat('D:\GitHub\D_FBA\Simulations\Gerosa\Results',num2str(w),'.xlsx'),intCor,5);
    xlswrite(strcat('D:\GitHub\D_FBA\Simulations\Gerosa\Results',num2str(w),'.xlsx'),NRMSE,6);
    xlswrite(strcat('D:\GitHub\D_FBA\Simulations\Gerosa\Results',num2str(w),'.xlsx'),upRecall,7);
    xlswrite(strcat('D:\GitHub\D_FBA\Simulations\Gerosa\Results',num2str(w),'.xlsx'),downRecall,8);
end
xlswrite('D:\GitHub\D_FBA\Simulations\Gerosa\Solstats.xlsx',sol_stat,1);
xlswrite('D:\GitHub\D_FBA\Simulations\Gerosa\Solstats.xlsx',sol_stat2,2);
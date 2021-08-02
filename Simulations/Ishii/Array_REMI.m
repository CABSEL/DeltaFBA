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
cobraParams.timeLimit = 600;
cobraParams.printLevel = 1;

%% Follow Remi script
%% This data from Elux2 paper
load('OtherTools\REMI\remi\models\iJO1366.mat') % this model is taken from BiGG database
load('OtherTools\REMI\remi\data\test_expr.mat') % Data is used from this: E-Flux2 and SPOT: Validated Methods for Inferring Intracellular Metabolic Flux Distributions from Transcriptomic Data
genes=expr{1};
genes_data=expr{2};

%% add binary use variables
fmodel=iJO1366;
fmodel.rev=ones(numel(fmodel.rev),1);
fmodel=addUseVariablesDH(fmodel);
fmodel=addNetFluxVariablesNEW(fmodel);
% we want to add cellular growth which is atleat 80 % of the maximum growth

sol=solveTFBAmodel(fmodel, false, 'gurobi_direct');

fmodel.var_lb(find(fmodel.f))=sol.val*0.8;

%% initialize varibales
for mutiter = 1:numel(expr_cols)-1
    
    wt=1;
    mut=mutiter+1; % this is for the pgm model
    tarns_cut=0.1; % transcript cutoff
    up=.05; % this cutoff indicates 5% of all reaction
    down=.05;
    
    %% find up and down regulated genes
    geneWT=nan(numel(fmodel.genes),1);
    geneMUT=nan(numel(fmodel.genes),1);
    for i=1:numel(fmodel.genes)
        [~,ind]=ismember(fmodel.genes{i},genes);
        if ind>0
            geneWT(i)=genes_data(ind,wt);
            geneMUT(i)=genes_data(ind,mut);
        end
    end
    
    % calculate ratio (if transcrip label is too samll then we dont cal culate ratio)
    ind1=find(geneWT(geneWT>tarns_cut));
    ind2=find(geneMUT(geneMUT>tarns_cut));
    indC=intersect(ind1,ind2); % common index
    Generatio=geneMUT(indC)./geneWT(indC);
    
    %% find up and downregulated reactions usig 5% up and down criteria
    [rxns,ratio]=evaluateGPR(fmodel,fmodel.genes(indC), Generatio,@geomean,@mean);
    indUP=find(ratio>1);
    ratioUP=ratio(indUP);
    indDOWN=find(ratio<1);
    ratioDOWN=ratio(indDOWN);
    [s_ratioUP,Iup]=sort(ratioUP);
    [s_ratioDOWN,Idown]=sort(ratioDOWN);
    cut1=floor(numel(Iup)*(1-up));
    ind1=indUP(Iup(cut1:end));
    cut2=ceil(numel(Idown)*down);
    ind2=indDOWN(Idown(1:cut2));
    regRxnRatio=[ratio(ind1);ratio(ind2)]; % this is the combined up and down reaction ratios and reactions
    regIdx=[ind1;ind2];
    
     %% R-D
    [~,t2] = ismember(rxns,iJO1366.rxns);
    R_full = zeros(numel(iJO1366.rxns),1);
    R_full(t2) = ratio;
    R_full(isnan(R_full)) = 0;
    
    %% create a Gex model which only intgerate relative expression
    % add constraints for up and down regulated reactions
    %     tmpModel=addRelExpCons2Models(model1,model2,model1.rxns(regIdx),regRxnRatio); % this can be also used addRelExpCons2Models
    %     tmpModel.description='Regulated Model';
    gex=addRelConsExpression(fmodel,fmodel,fmodel.rxns(regIdx),regRxnRatio);
    gex.description='REMI-Gex';
    sol=solveTFBAmodel(gex, false, 'gurobi_direct');
    
    MCS=sol.val; %% This is the maximum consistency score
    
    
    %% Alternative analysis
    model = gex;
    index_updown = [gex.relExp.forB];
    
    % Number of maximum solutions
    numsol = 300;
    
    z_matrix=[];
    sol_matrix=[];
    
    % Find binary objective value and set it as constrain
    sol=solveTFBAmodel(model,false,'gurobi_direct');
    model.var_lb(find(model.f))=sol.val;
    
    % Store the result
    z_matrix=[z_matrix sol.x(index_updown)];
    sol_matrix=[sol_matrix sol.x];
    store_obj=[];
    
    % Find alternate solutions
    obj=sol.val;
    store_obj(end+1)=obj;
    [num_cons,num_vars] = size(model.A);
    
    i = 1;
    j = 1;
    while (i>=1 && j==1)
        
        mat=[zeros(1,num_vars)];
        nonzero_ind=find(sol.x(index_updown)>0.98);
        mat(index_updown(nonzero_ind))=1;
        model.A=[model.A;mat];
        model.constraintNames{end+1}=['integercut' num2str(i)];
        model.constraintType((end+1)) = {'<'};
        model.rhs(end+1)=numel(nonzero_ind)-1;
        
        sol=solveTFBAmodel(model,false,'gurobi_direct');
        if (numel(sol.x)>0)
            z_matrix=[z_matrix sol.x(index_updown)];
            store_obj(end+1)=sol.val;
            sol_matrix=[sol_matrix sol.x];
            i = i+1;
        else
            j = 2;
        end
    end
    
    
    %% Load measure fluxes data
    [num,txt,~] = xlsread('Data\Ishii\Fluxomics.xlsx');
    flux_data = num(:,[1,8,9,14,20,23,4,5]);
    flux_list = txt(2:end,1);
    [~,fids] = ismember(flux_list, model.rxns);
    M_Wt = repmat(flux_data(:,wt),1,size(z_matrix,2));
    M_Exp = repmat(flux_data(:,mut), 1, size(z_matrix,2));
    
    rxns = iJO1366.rxns;
    [~,NFind]=ismember(strcat('NF_',rxns),model.varNames);
    [~,PNFind]=ismember(strcat('PERTURB_NF_',rxns),model.varNames);
    R_D = R_full(fids);
    %% Parfor here (still for Pearson correlation and percentage errors)
    parfor eId=1:size(z_matrix,2)
        model_tmp=gex;
        
        % Minimize the fluxes
        objIndex=[model_tmp.relExp.forB];
        activeIdx=objIndex(find(z_matrix(:,eId)>0.99));
        model_tmp.var_lb(activeIdx)=1;
        f=zeros(numel(model_tmp.f),1);
        f(NFind)=-1;
        f(PNFind)=-1;
        model_tmp.f=f;
        
        sol=solveTFBAmodel(model_tmp,false,'gurobi_direct');
        if (numel(sol.x)>0)
            
            P_Wt_Full = sol.x(NFind);
            P_Exp_Full = sol.x(PNFind);
            
            P_Wt(:,eId) = P_Wt_Full(fids);
            P_Exp(:,eId) = P_Exp_Full(fids);
            Wt_Full(:,eId) = P_Wt_Full;
            Exp_Full(:,eId) = P_Exp_Full;
            sol_pfba(eId) = sol.val;
        end
    end
    xlswrite(strcat(root,'\Simulations\Ishii\Array_REMI/',expr_cols{1,mut}, '.xlsx'),M_Wt,1);
    xlswrite(strcat(root,'\Simulations\Ishii\Array_REMI/',expr_cols{1,mut}, '.xlsx'),M_Exp,2);
    xlswrite(strcat(root,'\Simulations\Ishii\Array_REMI/',expr_cols{1,mut}, '.xlsx'),P_Wt,3);
    xlswrite(strcat(root,'\Simulations\Ishii\Array_REMI/',expr_cols{1,mut}, '.xlsx'),P_Exp,4);
    xlswrite(strcat(root,'\Simulations\Ishii\Array_REMI/',expr_cols{1,mut}, '.xlsx'),R_D,5);
    xlswrite(strcat(root,'\Simulations\Ishii\Array_REMI/',expr_cols{1,mut}, '.xlsx'),Wt_Full,6);
    xlswrite(strcat(root,'\Simulations\Ishii\Array_REMI/',expr_cols{1,mut}, '.xlsx'),Exp_Full,7);
    xlswrite(strcat(root,'\Simulations\Ishii\Array_REMI/',expr_cols{1,mut}, '.xlsx'),R_full,8);
    xlswrite(strcat(root,'\Simulations\Ishii\Array_REMI/',expr_cols{1,mut}, '.xlsx'),sol_pfba,9);
    
    clear P_Wt P_Exp Wt_Full Exp_Full sol_pfba
end

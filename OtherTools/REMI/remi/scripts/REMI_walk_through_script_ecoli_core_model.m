%######
% with this script we will genarte Gex, GexM and M model 

%######


addpath(genpath('/Users/vikash/Desktop/BackupEPFL/CPLEX_Studio1251')); % this is for cplex solver folder and one need to adjust the path accordingly 
% addpath(genpath('/Users/vikash/Desktop/BackupEPFL/FBA_Toolboxes')); % this is the path if one want apply thermodynamic constraints  
    
% we want to build relative model for ecoli for two conditions

% load ecoli core model
load('/Users/vikash/Desktop/REMI/models/e_coli_core.mat'); 

% 1) we will add relative expression constraint 
   
% we provide one random expression data    
% based on gene ratios (between two conditions) we need evallue GPR for finding reaction ratio
% first load relative expression data
load('/Users/vikash/Desktop/REMI/data/RelExpData.mat')
up_cutoff=2; % more than 2 fold up regulated is known as up regulated genes (This paramter can change)
down_cutoff=0.5;% less than 1/2 fold  is known as down regulated genes 
% we need to identify genes which are up and down regulated
% input argument : 1) model
%                  2) gene names
%                  3) gene ratios
%                  4) and opertaor (evalute based on GPR)
%                  5) or operator
% length of argument 2 and 3 should be same
[rxns,ratio]=evaluateGPR(e_coli_core,e_coli_core.genes,expData.ratios,@geomean,@mean);
% find up- and down- regulated reactions
indUP=find(ratio>up_cutoff);
ratioUP=ratio(indUP);
indDOWN=find(ratio<down_cutoff);
ratioDOWN=ratio(indDOWN);
regRxnRatio=[ratioUP;ratioDOWN];
regRxns=[rxns(indUP)';rxns(indDOWN)'];

% avoid numerical error (more than 100 fold is taken only 100)
regRxnRatio(regRxnRatio>100)=100;
regRxnRatio(regRxnRatio<1/100)=1/100;

% if we want to add relative constraint into TFA model then we need to add
% net flux variable to the model using 'addNetFluxVariablesNEW'
% and if one want to add into FBA model  then evalute scripts given below
mTmp=e_coli_core;
mTmp.rev=ones(numel(mTmp.rev),1);
use_m=addUseVariablesDH(mTmp);
netM=addNetFluxVariablesNEW(use_m);

%% We are going to add constraints for only relative expression
% now we add relative expression constarint 
% input argument : 1) model1 represents condition 1
%                  2) model2 represents condition2 
%                  3) regulated rxns 
%                  4) regulated reaction rations

% now we are build gex Model:  which integeartes relative gene expression
[gex]=addRelConsExpression(netM,netM,regRxns,regRxnRatio)
sol_gex=solveTFBAmodel(gex)

% maximum consistency score (MCS) will be objective value of the solution.
MCS=sol_gex.val; % this is the maximum consistency score of GeX model 

%% ALTERNATIVE ANALYSIS on a given MCS
% enumerate alternatives on expression variabels or for each MCS score
% there can be many alternative set of constraints and if one want to
% enumerate those alternatives one should use follwing script


path_save='/Users/vikash/Desktop/REMI/data/AltSolNew.mat' % This path will save all alternative results
% all binary variables are stored as model.relExp.forB
%                  3) will be model.relExp.forB

% we want to find alternative at maximum consistency then we need to add
% maximum consistency to the lowerbound
coM2=gex;
coM2.var_lb(end)=MCS; % this will force maximum consistency can not be 
time=200;
% input argument : 1) numsol is number of alternatives
%                  2) model
%                  3) index for binary variables 
%                  4) path for saving the resultes
%                  5) time for solver
findAltCombi(2,coM2,coM2.relExp.forB,path_save,time);

%% ADD RELATIVE METABOLITE on M and Gex Model
%
% first load relative data
load('/Users/vikash/Desktop/REMI/data/relMetdata.mat')

mratio=metData.ratio;
indUP=find(mratio>up_cutoff);
indDOWN=find(mratio<down_cutoff);

regMetRatio=[mratio(indUP);mratio(indDOWN)];
regMets=[e_coli_core.mets(indUP);e_coli_core.mets(indDOWN)];
%% We are going to add constraints for only relative metabolites (M model)
% now we add relative expression constarint 
% input argument : 1) model1 represents condition 1
%                  2) model2 represents condition2 
%                  3) regulated mets 
%                  4) regulated metratios
%                  5) false (for only metaboliet), true(if expression is alreay intgerated)
[M]=addRelMetabolite(netM,netM,regMets,regMetRatio,false)
sol_M=solveTFBAmodel(M) 

% you can force consistency by putting MCS in the lower bound. 
% coM2.var_lb(end)=MCS; % this will force maximum consistency can not be 


% find alternatives
path_save='/Users/vikash/Desktop/REMI/data/AltSolOnlyMet.mat'

findAltCombi(2,M,M.metB,path_save,time);

%% We want to add relative expression and metabolites together (GexM model)

[gexM]=addRelMetabolite(coM1,coM1,regMets,regMetRatio,true);
sol=solveTFBAmodel(gexM)
% find alternatives
path_save='/Users/vikash/Desktop/REMI/data/AltSolExpMet.mat'
comIdx=[gexM.relExp.forB;gexM.metB];
findAltCombi(2,gexM,comIdx,path_save,time);

%% force model for one alternative solution

modelTMP=coMExpmet;
sol1=solveTFBAmodel(modelTMP)
index1=comIdx(find(sol.x(comIdx)>0.98));
index0=comIdx(find(sol.x(comIdx)<0.08));
modelTMP.var_lb(index1)=1;
modelTMP.var_ub(index0)=0;
sol2=solveTFBAmodel(modelTMP);





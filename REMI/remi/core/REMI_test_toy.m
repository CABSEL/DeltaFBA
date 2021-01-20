% addPaths
addpath(genpath('/Users/vikashpandey/git/fba_toolbox'))
addpath(genpath('/Users/vikashpandey/git/CPLEX_Studio1251'))

% create toy metabolic toyModel
toyModel.mets={'A';'B';'C';'D';'E';'F'};
toyModel.rxns={'1';'2';'3';'4';'5';'6';'7';'8';'9'};
toyModel.S=[1 0 -1 0 0 0 0 0 0;
            0 1 1 -1 -1 0 0 0 0;
            0 0 0 1 0 -1 0 0 0;
            0 0 0 0 1 1 -1 -1 0;
            0 0 0 0 0 0 1 0 -1;
            0 0 0 0 0 0 0 1 -1;
         ];
toyModel.c=zeros(numel(toyModel.rxns),1);
toyModel.lb=zeros(numel(toyModel.rxns),1);
toyModel.ub=2*ones(numel(toyModel.rxns),1);
toyModel.c(end)=1;
toyModel.description='toy'
% toyModel.lb(1)=1;
% toyModel.lb(2)=1;
% % toyModel.lb(4)=1;
% % toyModel.ub(1)=1;
% % toyModel.ub(2)=1;
toyModel.b=zeros(numel(toyModel.mets),1);
toyModel.rev=ones(numel(toyModel.rxns),1);
sol1=solveFBACplex(toyModel);
% add use variables
mjoin1=addUseVariablesDH(toyModel);
mjoin2=addNetFluxVariables(mjoin1);

%% add Relative constarint
% addpath(genpath('/Users/vikashpandey/Documents/MATLAB/RelativeFirstAttempt/'));
DErxns={'3';'6';'8'};
DEratio=[1.4;0.6;1.2];
coM=addRelConsExpression(mjoin2,mjoin2,DErxns,DEratio)
sol=solveTFBAmodel(coM);

on=coM.varNames(DEmodel.objIndex1B(find(sol.x(DEmodel.objIndex1B)>0.99)))
forward=strcat('F_',DErxns);
pforward=strcat('PERTURB_F_',DErxns);
[~,ind1]=ismember(forward,coM.varNames);
[~,ind2]=ismember(pforward,coM.varNames);
[nfind,~]=getAllVar(coM,{'NF'});
[pnfind,~]=getAllVar(coM,{'PERTURB_NF'});


[coM.varNames(ind1) num2cell(sol.x(ind1)) coM.varNames(ind2) num2cell(sol.x(ind2)) num2cell(DEratio) ]
[coM.rxns,num2cell(sol.x(nfind)) num2cell(sol.x(pnfind))]


%% add relative metabolite 
dmets={'A';'D';'E'}
dmratio=[0.25;4;0.25]
% dmets={'C'}
% dmratio=[2]
[coM2]= addRelMetabolite(mjoin2,mjoin2,dmets,dmratio,false)
sol=solveTFBAmodel(coM2);

on=coM.varNames(DEmodel.objIndex1B(find(sol.x(DEmodel.objIndex1B)>0.99)))
forward=strcat('F_',DErxns);
pforward=strcat('PERTURB_F_',DErxns);
[~,ind1]=ismember(forward,coM.varNames);
[~,ind2]=ismember(pforward,coM.varNames);
[nfind,~]=getAllVar(coM,{'NF'});
[pnfind,~]=getAllVar(coM,{'PERTURB_NF'});


[coM.varNames(ind1) num2cell(sol.x(ind1)) coM.varNames(ind2) num2cell(sol.x(ind2)) num2cell(DEratio) ]
[coM.rxns,num2cell(sol.x(nfind)) num2cell(sol.x(pnfind))]
[sol.x(coM2.metVars.cond1.P) sol.x(coM2.metVars.cond1.C) sol.x(coM2.metVars.cond2.P) sol.x(coM2.metVars.cond2.C)]

%% add relative metabolites on core model
load('/Users/vikashpandey/Documents/MATLAB/REMIWalkThrough/Model/e_coli_core.mat')
core=e_coli_core;
core.rev=ones(numel(core.rev),1);
mjoin1=addUseVariablesDH(core);
mjoin2=addNetFluxVariables(mjoin1);
idx=[1:10 20:25 30:35];
dmets=core.mets(idx);
dmratio=rand(numel(idx),1)*10;
% dmets={'C'}
% dmratio=[2]
[coM2]= addRelMetabolite(mjoin2,mjoin2,dmets,dmratio,false)
sol=solveTFBAmodel(coM2);
[sol.x(coM2.metVars.cond1.P) sol.x(coM2.metVars.cond1.C) sol.x(coM2.metVars.cond2.P) sol.x(coM2.metVars.cond2.C)];
 
% Find alternatives for metabolites
path_save='/Users/vikashpandey/Documents/MATLAB/REMIWalkThrough/Model/alt.mat';
findAltCombi(10,coM2,coM2.metB,path_save)

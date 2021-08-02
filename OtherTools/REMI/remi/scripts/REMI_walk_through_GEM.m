% ##################

% here we are showing how to build models without thermodynamics: REMI-Gex, REMI-GexM, REMI-M
% with thermodynamics REMI-TGex, REMI-TGexM, REMI-TM
%

%#####################

%% addPaths for cplex solver and FBA Toolboxes

% addpath(genpath('/Users/vikashpandey/git/CPLEX_Studio1251'));
% addpath(genpath('/Users/vikashpandey/git/FBA_Toolboxes'));
addpath(genpath('/Users/vikash/Desktop/REMI'));
% store data
load('/Users/vikash/Desktop/REMI/data/fmodel.mat')
fff=fmodel;

%% initialize varibales
wt=1;
mut=2; % this is for the pgm model
tarns_cut=0.1; % transcript cutoff
up=.05; % this cutoff indicates 5% of all reaction
down=.05;

%% This data from Elux2 paper
load('/Users/vikash/Desktop/REMI/models/iJO1366.mat') % this model is taken from BiGG database
load('/Users/vikash/Desktop/REMI/data/test_expr.mat') % Data is used from this: E-Flux2 and SPOT: Validated Methods for Inferring Intracellular Metabolic Flux Distributions from Transcriptomic Data
genes=expr{1};
genes_data=expr{2};

%% add binary use variables
fmodel=iJO1366;
fmodel.rev=ones(numel(fmodel.rev),1);
fmodel=addUseVariablesDH(fmodel);
fmodel=addNetFluxVariablesNEW(fmodel);
% we want to add cellular growth which is atleat 80 % of the maximum growth

% sol=solveTFBAmodel(fmodel);

fmodel.var_lb(find(fmodel.f))=sol.val*0.8;

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

%% create a Gex model which only intgerate relative expression
% add constraints for up and down regulated reactions
%     tmpModel=addRelExpCons2Models(model1,model2,model1.rxns(regIdx),regRxnRatio); % this can be also used addRelExpCons2Models
%     tmpModel.description='Regulated Model';
gex=addRelConsExpression(fmodel,fmodel,fmodel.rxns(regIdx),regRxnRatio);
gex.description='REMI-Gex';
% sol=solveTFBAmodel(gex);

MCS=sol.val %% This is the maximum consistency score

%% This piece of codes for adding relative metabolites data and creating GexM model
% add relative Expression into models


mutant=1;% concentartion data for mutant (pgm)
concWT=fff.Ishii.concData(:,8); % concentartion data for wild type
concMUT=fff.Ishii.concData(:,mutant);
metratio=concMUT./concWT;

% identify up and down regulated metabolites
indUP=find(metratio>1);
ratioUP=fff.mets(indUP);
indDOWN=find(metratio<1);
ratioDOWN=fff.mets(indDOWN);
[s_ratioUP,Iup]=sort(ratioUP);
[s_ratioDOWN,Idown]=sort(ratioDOWN);
cut1=floor(numel(Iup)*(1-up));
ind1=indUP(Iup(cut1:end));
cut2=ceil(numel(Idown)*down);
ind2=indDOWN(Idown(1:cut2));
regMetRatio=[metratio(ind1);metratio(ind2)];
regIdx=[ind1;ind2];

% build REMI-Gex model :add relative metabolomics data
gexM=addRelMetabolite(gex,gex,fff.mets(regIdx),regMetRatio,true) % true when we want to add metabolomics into a Gex Model
gexM.description='REMI-GexM';
% sol=solveTFBAmodel(gexM);
%% create only REMI-M model
% we need to use fmodel which is without intgerating relative
% expression
remiM=addRelMetabolite(fmodel,fmodel,fff.mets(regIdx),regMetRatio,false); % false when we want to add metabolomics only
% store models
remiM.description='REMI-M'; 
% sol=solveTFBAmodel(remiM);

%% we can also apply REMI on thermodynamic model 
load('/Users/vikash/Desktop/REMI/models/iJO1366_thermo.mat') % this is thermodynamic model
fmodel_thermo=iJO1366_thermo;
% sol=solveTFBAmodel(fmodel_thermo);
sol_val=0.8
fmodel_thermo.var_lb(find(fmodel_thermo.f))=sol.val*0.8;

%% create REMI-TGex model which intgerates relative expression in a thermodynamic consistent model
tgex=addRelConsExpression(fmodel_thermo,fmodel_thermo,fmodel_thermo.rxns(regIdx),regRxnRatio);
tgex.description='REMI-TGex';
% sol=solveTFBAmodel(gex);

%% create REMI-TGex model which intgerates relative expression, relative metabolomics into a thermodynamic consistent model
tgexM=addRelMetabolite(tgex,tgex,fff.mets(regIdx),regMetRatio,true) % true when we want to add metabolomics into a Gex Model
tgexM.description='REMI-tGexM';
% sol=solveTFBAmodel(gexM);

%% create REMI-TGex model which intgerates relative expression, relative metabolomics into a thermodynamic consistent model
tM=addRelMetabolite(fmodel_thermo,fmodel_thermo,fff.mets(regIdx),regMetRatio,false) % true when we want to add metabolomics into a Gex Model
tM.description='REMI-TM';
% sol=solveTFBAmodel(gexM);

%%  #### ALTERNATIVE ANALYSIS one should use findAltCombi function
% we shown in ecoli walk thorugh file
%time=500;
% input argument : 1) numsol is number of alternatives
%                  2) model
%                  3) index for binary variables 
%                  4) path for saving the resultes
%                  5) time for solver
% findAltCombi(2,coM2,coM2.relExp.forB,path_save,time);


% ####

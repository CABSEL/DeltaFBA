function [coM]= addSumPC(coM,metIdx,cond,para)
% In this function we want to add sum of production and some of consumption
% for each metaboliets
% coM is the combined relative model
% metIdx is metabolite index in the model
% cond is 1 or 2
S=coM.S;
pId=find(S(metIdx,:)>0);  % production reactions
cId=find(S(metIdx,:)<0); % consumption reactions
pcoeff=abs(S(metIdx,pId)); % coefficient of productions
scoeff=abs(S(metIdx,cId)); % coeffiecient of consumptions
% now we want to sum forward and backward reactions
% find NF indexes and PERTURB indexes
coM.metVars.C=[];
coM.metVars.P=[];
pcondF=[];
pcondR=[];
scondF=[];
scondR=[];
if cond==1
    if numel(pId)>0 && numel(cId)==0
        [~,pcondF]=ismember(strcat('F_',coM.rxns(pId)),coM.varNames);
        [~,pcondR]=ismember(strcat('R_',coM.rxns(pId)),coM.varNames);
    elseif numel(pId)==0 && numel(cId)>0
        [~,scondF]=ismember(strcat('F_',coM.rxns(cId)),coM.varNames);
        [~,scondR]=ismember(strcat('R_',coM.rxns(cId)),coM.varNames);
    elseif numel(pId)>0 && numel(cId)>0
        [~,pcondF]=ismember(strcat('F_',coM.rxns(pId)),coM.varNames);
        [~,pcondR]=ismember(strcat('R_',coM.rxns(pId)),coM.varNames);
        [~,scondF]=ismember(strcat('F_',coM.rxns(cId)),coM.varNames);
        [~,scondR]=ismember(strcat('R_',coM.rxns(cId)),coM.varNames);
    else
        disp('problem in production and consumption')
        return
    end
    
elseif cond==2
    if numel(pId)>0 && numel(cId)==0
        [~,pcondF]=ismember(strcat('PERTURB_F_',coM.rxns(pId)),coM.varNames);
        [~,pcondR]=ismember(strcat('PERTURB_R_',coM.rxns(pId)),coM.varNames);
    elseif numel(pId)==0 && numel(cId)>0
        [~,scondF]=ismember(strcat('PERTURB_F_',coM.rxns(cId)),coM.varNames);
        [~,scondR]=ismember(strcat('PERTURB_R_',coM.rxns(cId)),coM.varNames);
   elseif numel(pId)>0 && numel(cId)>0
        [~,pcondF]=ismember(strcat('PERTURB_F_',coM.rxns(pId)),coM.varNames);
        [~,pcondR]=ismember(strcat('PERTURB_R_',coM.rxns(pId)),coM.varNames);
        [~,scondF]=ismember(strcat('PERTURB_F_',coM.rxns(cId)),coM.varNames);
        [~,scondR]=ismember(strcat('PERTURB_R_',coM.rxns(cId)),coM.varNames);
    else
        disp('problem in production and consumption')
        return
    end
    
    
else
    disp ('condition can be 1 or 2')
    return
end

% create variable for sum of production in coditions 1
coM = addNewVariableInTFA(coM, strcat('Cond',num2str(cond),'SumFlux_P_', coM.mets{metIdx} ), 'C', [0 para.bigM]);
bIdx1=length(coM.varNames);
coM.metVars.P=[coM.metVars.P;bIdx1];
% add constarint for production of condition 1
CLHS.varIDs=[pcondF;scondR;bIdx1];

CLHS.varCoeffs=[pcoeff scoeff -1]';

coM=addNewConstraintInTFA(coM, strcat('Cond',num2str(cond),'CSumFlux_P_',coM.mets{metIdx}), '=', CLHS, 0);

% create variable for sum of production in coditions 1
coM = addNewVariableInTFA(coM, strcat('Cond',num2str(cond),'SumFlux_C_', coM.mets{metIdx} ), 'C', [0 para.bigM]);
bIdx2=length(coM.varNames);
%store indexes
coM.metVars.C=[coM.metVars.C;bIdx2];
% add constarint for production of condition 1
CLHS.varIDs=[scondF;pcondR;bIdx2];
CLHS.varCoeffs=[scoeff pcoeff -1]';
coM=addNewConstraintInTFA(coM, strcat('Cond',num2str(cond),'CSumFlux_C_',coM.mets{metIdx}), '=', CLHS, 0);
coM.PC={pcondF pcondR scondF scondR};
    
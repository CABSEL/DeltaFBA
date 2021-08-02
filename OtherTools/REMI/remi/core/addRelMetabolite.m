function [coM]= addRelMetabolite(model1,model2,regMets,regRatio,combined)
% if we want to integerate only metabolites then we need to
% add 2 models and make a relative models
% model in TFA format
% regMets: cells of regulated metabolites
% regMets: regualated ratios
% combined: test whether are you intgerating metabolite and expression together (True: we added relative expression and we want to add metabolites; False for adding only metaboliets ) 
%% testing whether model has expression integerated or not

if ~combined %do i wanna do thermo ?
    coM=combine2Models(model1,model2);
else
    coM=model1;
end

%% parameters
para.bigM=100000;
para.eps=0.00001;
para.basal=0.0001;
%% Here we want to store consumption and production variabels.
coM.metVars.cond1.P=[];
coM.metVars.cond1.C=[];
coM.metVars.cond2.P=[];
coM.metVars.cond2.C=[];


% find stoichiometry

metB=[];
[~,metIdx]=ismember(regMets,coM.mets);
metIdx=metIdx(metIdx>0);
metRatio=regRatio(metIdx>0);
%% now we want to add relative constraints for metabolites
for i=1:numel(metIdx)
    %% now we want to add for sum of production and consumption for all metabolites in condition 2
    [coM]= addSumPC(coM,metIdx(i),1,para);
    coM.metVars.cond1.P=[coM.metVars.cond1.P;coM.metVars.P];
    coM.metVars.cond1.C=[coM.metVars.cond1.C;coM.metVars.C];
    PC1=coM.PC;
    [coM]= addSumPC(coM,metIdx(i),2,para);
    coM.metVars.cond2.P=[coM.metVars.cond2.P;coM.metVars.P];
    coM.metVars.cond2.C=[coM.metVars.cond2.C;coM.metVars.C];
    PC2=coM.PC;
    %% now we want to add constraints for relative metabolites
    % res=[V ph1 V ph2 slack]
    coM = addNewVariableInTFA(coM, strcat('useB_', coM.mets{metIdx(i)}), 'B', [0 1]);
    useB=length(coM.varNames);
    [coM,resP]= addRelSumFlux(coM,coM.varNames{coM.metVars.cond1.P(i)},coM.varNames{coM.metVars.cond2.P(i)},useB,para);
    
    [coM,resC]= addRelSumFlux(coM,coM.varNames{coM.metVars.cond1.C(i)},coM.varNames{coM.metVars.cond2.C(i)},useB,para);
    
    % apply inherent flux constraint which really tess if a metabolite is up
    % regulated then each flux which are participating in production can
    % not be less than in condition 2 comapred to condition 1
   
    
    
    %%
    % add metaboliets
    
    % create binaray varaibales for the consistency
    coM = addNewVariableInTFA(coM, strcat('MET_', coM.mets{metIdx(i)}), 'B', [0 1]);
    bIdx=length(coM.varNames);
    metB=[metB;bIdx];
    % add basal flux to all sum of produsts and consumptions
    [coM]= addbasalSumFlux(coM,resP.ph1,resP.ph2,resC.ph1,resC.ph2,bIdx,para);
    % add slack variabels
    coM = addNewVariableInTFA(coM, strcat('SLACK_', coM.mets{metIdx(i)}), 'C', [0 para.bigM]);
    slack=length(coM.varNames);
    [coM,res]= addSlackB(coM,coM.varNames{slack},useB,para);
    slackV=res.v; % this is the new slack variable
    
%     coM = addNewVariableInTFA(coM, strcat('useB1_', coM.mets{metIdx(i)}), 'B', [0 1]);
%     useB1=length(coM.varNames);
    
    if metRatio(i)>=1
        % mian objective: ph2>r*ph1
        % we want more production
        
        
        
        % v2 -V1*r+slackV>0
        CLHS.varIDs=[resP.v2 resP.v1 slackV];
        CLHS.varCoeffs=[1 -metRatio(i) 1];
        coM=addNewConstraintInTFA(coM, strcat('Rel_P_up',coM.mets{metIdx(i)}), '>', CLHS, 0);
         
        % ph2 -ph1*(1/r)' -V2 + v1*(1/r)-slackV>0
        CLHS.varIDs=[resC.ph2 resC.ph1 resC.v2 resC.v1 slackV];
        CLHS.varCoeffs=[1 -1/metRatio(i) -1 1/metRatio(i) -1 ];
        coM=addNewConstraintInTFA(coM, strcat('Rel_C_up',coM.mets{metIdx(i)}), '<', CLHS, 0);
        coM=applyInherentCons(coM,PC1,PC2,coM.mets{metIdx(i)},'up',useB,para);
        
    else
        % v2 -V1*r-slackV<0
        CLHS.varIDs=[resP.v2 resP.v1 slackV];
        CLHS.varCoeffs=[1 -metRatio(i) -1];
        coM=addNewConstraintInTFA(coM, strcat('Rel_P_down',coM.mets{metIdx(i)}), '<', CLHS, 0);
        
        % ph2 -ph1*(1/r)' -V2 + v1*(1/r)+slackV>0
        CLHS.varIDs=[resC.ph2 resC.ph1 resC.v2 resC.v1 slackV];
        CLHS.varCoeffs=[1 -1/metRatio(i) -1 1/metRatio(i) 1 ];
        coM=addNewConstraintInTFA(coM, strcat('Rel_C_down',coM.mets{metIdx(i)}), '>', CLHS, 0);
        coM=applyInherentCons(coM,PC1,PC2,coM.mets{metIdx(i)},'down',useB,para);
    end
    % eps*(1-YBplus)<=slack<=eps+(1-YBplus)*bigM
    % -eps*yBplus-slack<=-eps and
    CLHS.varIDs=[bIdx slack];
    CLHS.varCoeffs=[-para.eps -1];
    coM=addNewConstraintInTFA(coM, strcat('SlackMet1_',coM.mets{metIdx(i)}), '<', CLHS, -para.eps);
    
    %slack+yBplus*bigM<=eps+bigM
    CLHS.varIDs=[bIdx slack];
    CLHS.varCoeffs=[para.bigM 1];
    coM=addNewConstraintInTFA(coM, strcat('SlackMet2_',coM.mets{metIdx(i)}), '<', CLHS, para.eps+para.bigM);
    
    
    
end

%% add sum of metabolites
if ~combined %do i wanna do thermo ?
    coM = addNewVariableInTFA(coM, 'Sum_MREG', 'I', [0 Inf]);
    bIdx1=length(coM.varNames);
    % add sum of constarints
    CLHS.varIDs=[metB' bIdx1];
    CLHS.varCoeffs=[ones(1,numel(metB))*-1 1];
    coM=addNewConstraintInTFA(coM, 'REGSumM', '=', CLHS, 0);
    coM.metB=metB;
    f=zeros(numel(coM.varNames),1);
    f(end)=1;
    
    coM.objtype = -1; %-1; % 1 - minimize, -1 - maximize
    coM.f = f;
else
    coM = addNewVariableInTFA(coM, 'Sum_MREG', 'I', [0 Inf]);
    bIdx1=length(coM.varNames);
    % add sum of constarints
    if isfield(coM,'objIndex1B')
        CLHS.varIDs=[coM.objIndex1B' metB' bIdx1];
        CLHS.varCoeffs=[ones(1,numel(coM.objIndex1B))*-1 ones(1,numel(metB))*-1 1];
    elseif isfield(coM,'relExp')
         CLHS.varIDs=[coM.relExp.forB' metB' bIdx1];
         CLHS.varCoeffs=[ones(1,numel(coM.relExp.forB'))*-1 ones(1,numel(metB))*-1 1];
    else
        disp ('there is error in combining')
    end
    coM=addNewConstraintInTFA(coM, 'REGSumM', '=', CLHS, 0);
    coM.metB=metB;
    f=zeros(numel(coM.varNames),1);
    f(end)=1;
    
    coM.objtype = -1; %-1; % 1 - minimize, -1 - maximize
    coM.f = f;
end

function [coM]=addRelConsExpression(model1,model2,regRxns,regRatio)
    
    % model1: in condition1 cobra Or TFA format
    % model2: in condition1 cobra Or TFA format
    % regRxns: reaction for which expression is pertubed
    % regRatio: regulated reaction ratios
    
    %% parameters 
    para.bigM=100000;
    para.eps=0.00001;
    para.basal=0.0025;%0.0001;
    %%
    % now we are going to add two models whereas one represent condition1
    % and another represents condition2 
    coM=combine2Models(model1,model2);
    
    % indexes for expression analysis analysis  
    coM.relExp.forC=[]; % this is the index for forward and backward reactions (continuous varaiable)
    coM.relExp.forB=[]; % this is the index for forward and backward reactions (Binary varaiable)

    
    for i=1:numel(regRxns)
        [~,Find]=ismember(['F_' regRxns{i}],coM.varNames);
        [~,PFind]=ismember(['PERTURB_F_' regRxns{i}],coM.varNames);
        
        [coM]=addConsEachDir(coM,Find,PFind,'F',regRxns{i},regRatio(i),para);
        
        [~,Rind]=ismember(['R_' regRxns{i}],coM.varNames);
        [~,PRind]=ismember(['PERTURB_R_' regRxns{i}],coM.varNames);
        [coM]=addConsEachDir(coM,Rind,PRind,'R',regRxns{i},regRatio(i),para);
    end
                
    % add sum of the variabels
    coM = addNewVariableInTFA(coM, 'Sum_REG', 'I', [0 Inf]);
    bIdx1=length(coM.varNames);
    % add sum of constarints
    CLHS.varIDs=[coM.relExp.forB' bIdx1];
    CLHS.varCoeffs=[ones(1,numel(coM.relExp.forB))*-1 1];
    coM=addNewConstraintInTFA(coM, 'REGSum', '=', CLHS, 0);
    
    f=zeros(numel(coM.varNames),1);
    f(end)=1;
    
    coM.objtype = -1; %-1; % 1 - minimize, -1 - maximize
    coM.f = f;
    
    
end



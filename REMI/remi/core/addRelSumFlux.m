function [coM,reqVarIdx]= addRelSumFlux(coM,cond1Var,cond2Var,useB,para)

% coM is the combined relative model
% cond1Var and cond2Var for each metabolite 
% para is parameters
%% now we want to add basal flux

% we need to generate varaibels
% V=sumFlux*useB, %% V, useB, slack,
[~,ph1]=ismember(cond1Var,coM.varNames); % sum of flux
[~,ph2]=ismember(cond2Var,coM.varNames);
% for condition1
coM = addNewVariableInTFA(coM, strcat(cond1Var, '_contV' ), 'C', [0 para.bigM]);
contV1=length(coM.varNames);
% coM = addNewVariableInTFA(coM, strcat(cond1Var, '_useB' ), 'B', [0 1]);
% useB1=length(coM.varNames);

%
% for condition2
coM = addNewVariableInTFA(coM, strcat(cond2Var, '_contV' ), 'C', [0 para.bigM]);
contV2=length(coM.varNames);
% coM = addNewVariableInTFA(coM, strcat(cond2Var, '_useB' ), 'B', [0 1]);
% useB2=length(coM.varNames);



% [coM]= add3ConsBinCont(coM,contV1,ph1, useB1,para) % V=C*B V: contV1 C=ph1 and B=useB1
% [coM]= add3ConsBinCont(coM,contV2,ph2, useB2,para) % V=C*B V: contV1 C=ph1 and B=useB1
[coM]= add3ConsBinCont(coM,contV1,ph1, useB,para); % V=C*B V: contV1 C=ph1 and B=useB1
[coM]= add3ConsBinCont(coM,contV2,ph2, useB,para); % V=C*B V: contV1 C=ph1 and B=useB1



reqVarIdx.v1=contV1; 
reqVarIdx.ph1=ph1; 
reqVarIdx.v2=contV2; 
reqVarIdx.ph2=ph2; 


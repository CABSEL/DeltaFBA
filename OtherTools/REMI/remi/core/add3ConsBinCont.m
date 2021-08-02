function [coM]= add3ConsBinCont(coM,tIdx,cIdx, bIdx,para)
%%
% tIdx is the index of tarnsformed variable (V=C*B)
% cIdx is the index of continuous variables (C)
%% bIdx is the binary variable
% V-MB<=0
CLHS.varIDs=[tIdx bIdx];
CLHS.varCoeffs=[1 -para.bigM];
coM=addNewConstraintInTFA(coM, strcat('c1',coM.varNames{tIdx}  ), '<', CLHS, 0);
% V-C<=0
CLHS.varIDs=[tIdx cIdx];
CLHS.varCoeffs=[1 -1];
coM=addNewConstraintInTFA(coM, strcat('c2',coM.varNames{tIdx}  ), '<', CLHS, 0);

% V-C-MB>=-M
CLHS.varIDs=[tIdx cIdx bIdx];
CLHS.varCoeffs=[1 -1 -para.bigM];
coM=addNewConstraintInTFA(coM, strcat('c3',coM.varNames{tIdx} ), '>', CLHS, -para.bigM);


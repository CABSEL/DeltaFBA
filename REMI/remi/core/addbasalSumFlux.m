function [coM]= addbasalSumFlux(coM,pcond1Var,pcond2Var,scond1Var,scond2Var,bIdx,para)

% add basal reaction
% for products
% ph1-basalB>=0
CLHS.varIDs=[pcond1Var bIdx];
CLHS.varCoeffs=[1 -para.basal];
coM=addNewConstraintInTFA(coM, strcat('basal',coM.varNames{pcond1Var} ), '>', CLHS, 0);

% ph2-basalB>=0
CLHS.varIDs=[pcond2Var bIdx];
CLHS.varCoeffs=[1 -para.basal];
coM=addNewConstraintInTFA(coM, strcat('basal',coM.varNames{pcond2Var} ), '>', CLHS, 0);


% for substarates
% ph1-basalB>=0
CLHS.varIDs=[scond1Var bIdx];
CLHS.varCoeffs=[1 -para.basal];
coM=addNewConstraintInTFA(coM, strcat('basal',coM.varNames{scond1Var} ), '>', CLHS, 0);

% ph2-basalB>=0
CLHS.varIDs=[scond2Var bIdx];
CLHS.varCoeffs=[1 -para.basal];
coM=addNewConstraintInTFA(coM, strcat('basal',coM.varNames{scond2Var} ), '>', CLHS, 0);

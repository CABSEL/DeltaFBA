function [newM]=combine2Models(model1,model2)
% model1 represents  condition1
% model2 represents condition2
% newM is the combined model
newM=model1;
if ~isfield(newM,'var_lb')
    disp ('model first need to make in TFA format')
%     newM.var_ub = [model1.ub;model2.ub];
%     newM.var_lb = [model1.lb;model2.lb];
%     newM.varNames = [strcat('NF_', model1.rxns);strcat('PERTURB_NF_', model2.rxns)];
%     newM.vartypes = [repmat({'C'},numel(model1.rxns),1);repmat({'C'},numel(model2.rxns),1)];
%     newM.constraintNames=[strcat('M_', model1.mets);strcat('M_', model2.mets)];
%     newM.constraintType=[repmat({'='},numel(model1.mets),1);repmat({'='},numel(model2.mets),1)];
%     newM.A=sparse(numel(constraintNames), numel(varNames));
%     newM.A(1:numel(model1.mets),1:numel(model1.rxns))=model1.S;
%     newM.A(numel(model1.mets)+1:numel(constraintNames),numel(model1.rxns)+1:numel(varNames))=model2.S;
%     newM.rhs=[zeros(numel(model1.mets),1);zeros(numel(model2.mets),1)];
    
else
    
    % initialize model variables for control and perturbed conditions
    newM.var_ub = [model1.var_ub;model2.var_ub];
    newM.var_lb = [model1.var_lb;model2.var_lb];
    newM.varNames =[model1.varNames;strcat('PERTURB_',model2.varNames)];
    newM.vartypes = [model1.vartypes;model2.vartypes];
    newM.constraintNames=[model1.constraintNames;strcat('PERTURB_',model2.constraintNames)];
    newM.constraintType=[model1.constraintType;model2.constraintType];
    newM.A=sparse(numel(newM.constraintNames), numel(newM.varNames));
    newM.A(1:numel(model1.constraintNames),1:numel(model1.varNames))=model1.A;
    newM.A(numel(model1.constraintNames)+1:numel(newM.constraintNames),numel(model1.varNames)+1:numel(newM.varNames))=model2.A;
    newM.rhs=[model1.rhs;model2.rhs];
end
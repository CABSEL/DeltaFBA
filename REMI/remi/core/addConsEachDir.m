function [coM]=addConsEachDir(coM,Find,PFind,var,rxn,ratio,para)

% This function add relative constariants as well as basal fluxes
% coM is the combined relative model
% Find is the index for the flux in condition 1
% PFind is the index for the flux in condition 2
% var can be 'F' or  'R' 
% rxn is the reaction name
% para is parameter
if (Find>0) & (PFind>0)
    
    % frist create variabels for Forward
    coM = addNewVariableInTFA(coM, strcat('Yplus_', var, '_',rxn ), 'C', [0 para.bigM]);
    bIdx1=length(coM.varNames);
    % create bianry variable
    coM = addNewVariableInTFA(coM, strcat('YBplus_', var, '_',rxn ), 'B', [0 1]);
    bIdx1B=length(coM.varNames);
    % now we want to comapre 2 expressions
    if ratio>=1
        % F2<F1*ratio + Yplus
        % F2-F1*ratio-Yplus <0
        
%         CLHS.varIDs=[Find PFind bIdx1];
%         CLHS.varCoeffs=[-ratio 1 -1];
%         coM=addNewConstraintInTFA(coM, strcat('Cplus_', var, '_',rxn), '<', CLHS, 0);
                % F2 >F1*ratio-Yplus
                % F2-F1*raio+Yplus>0
        CLHS.varIDs=[Find PFind bIdx1];
        CLHS.varCoeffs=[-ratio 1 1];
        coM=addNewConstraintInTFA(coM, strcat('Cplus_', var, '_',rxn), '>', CLHS, 0);
    else
        % F2> F1*ratio-Yplus
        % F2-F1*raio+ Yplus>0
%         CLHS.varIDs=[Find PFind bIdx1];
%         CLHS.varCoeffs=[-ratio 1 1];
%         coM=addNewConstraintInTFA(coM, strcat('Cplus_', var, '_',rxn), '>', CLHS, 0);
                   % F2< F1*ratio+Yplus
                    % F2-Yplus-F1*raio<0
        CLHS.varIDs=[Find PFind bIdx1];
        CLHS.varCoeffs=[-ratio 1 -1];
        coM=addNewConstraintInTFA(coM, strcat('Cplus_', var, '_',rxn), '<', CLHS, 0);
    end
    end
    
    % eps*(1-YBplus)<=yPlus<=eps+(1-YBplus)*bigM
    % -eps*yBplus-yPlus<=-eps and
    CLHS.varIDs=[bIdx1B bIdx1];
    CLHS.varCoeffs=[-para.eps -1];
    coM=addNewConstraintInTFA(coM, strcat('C1plus_', var, '_',rxn), '<', CLHS, -para.eps);
    
    %yplus+yBplus*bigM<=eps+bigM
    CLHS.varIDs=[bIdx1B bIdx1];
    CLHS.varCoeffs=[para.bigM 1];
    coM=addNewConstraintInTFA(coM, strcat('C2plus_', var, '_',rxn), '<', CLHS, para.eps+para.bigM);
    % here we want to store the reactions
    coM.relExp.forC=[coM.relExp.forC;bIdx1];
    coM.relExp.forB=[coM.relExp.forB;bIdx1B];
    % mow we are adding basal flux for each fluxes
    % F1>delta*YBplus_F
    % F1-delta*YBplus_F >0
    CLHS.varIDs=[Find bIdx1B];
    CLHS.varCoeffs=[1 -para.basal];
    coM=addNewConstraintInTFA(coM, strcat('basal1_',var, '_',rxn), '>', CLHS, 0);
    % F2-delta*YBplus_F >0
    CLHS.varIDs=[PFind bIdx1B];
    CLHS.varCoeffs=[1 -para.basal];
    coM=addNewConstraintInTFA(coM, strcat('basal2_', var, '_',rxn), '>', CLHS, 0);
    
end
function [coM]=applyInherentCons(coM,fcond1,fcond2,metVar,sign,useB,para)
% here we want to add add constraints
% coM is the combined model
% flux index of condition 1
% flux index of condition 2
% sign can be < or >
P=[1 3];
C=[2 4];

if isequal(sign,'up')
    
    for i=1:numel(P)
        % for production
        % cond2-cond1>=0
        
        for j=1:numel(fcond2{P(i)})
            
            % add new variables 
            coM = addNewVariableInTFA(coM, strcat(metVar,'_inherP_',num2str(i),'_',num2str(j),'cond1' ), 'C', [0 para.bigM]);
            contV1=length(coM.varNames);
            % nevar=fluxvar*B
            [coM]= add3ConsBinCont(coM,contV1,fcond1{P(i)}(j), useB,para);
            % add new variables 
            coM = addNewVariableInTFA(coM, strcat(metVar,'_inherP_',num2str(i),'_',num2str(j),'cond2' ), 'C', [0 para.bigM]);
            contV2=length(coM.varNames);
            % nevar=fluxvar*B
            [coM]= add3ConsBinCont(coM,contV2,fcond2{P(i)}(j), useB,para);
            % add new variables 
            CLHS.varIDs=[contV2 contV1];
            CLHS.varCoeffs=[1 -1];
            coM=addNewConstraintInTFA(coM, strcat(metVar,'_inherP_',num2str(i),'_',num2str(j) ), '>', CLHS, 0);
            
        end
    end
    % for consumption
    % cond2-cond1>=0
    for i=1:numel(C)
        for j=1:numel(fcond2{C(i)})
             % add new variables 
            coM = addNewVariableInTFA(coM, strcat(metVar,'_inherC_',num2str(i),'_',num2str(j),'cond1' ), 'C', [0 para.bigM]);
            contV1=length(coM.varNames);
            % nevar=fluxvar*B
            [coM]= add3ConsBinCont(coM,contV1,fcond1{C(i)}(j), useB,para);
            % add new variables 
            coM = addNewVariableInTFA(coM, strcat(metVar,'_inherC_',num2str(i),'_',num2str(j),'cond2' ), 'C', [0 para.bigM]);
            contV2=length(coM.varNames);
            % nevar=fluxvar*B
            [coM]= add3ConsBinCont(coM,contV2,fcond2{C(i)}(j), useB,para);
            % add new variables 
            CLHS.varIDs=[contV2 contV1 fcond2{C(i)}(j) fcond1{C(i)}(j)];
            CLHS.varCoeffs=[-1 1 1 -1];
            coM=addNewConstraintInTFA(coM, strcat(metVar,'_inherC_',num2str(i),'_',num2str(j) ), '<', CLHS, 0);
        end
    end
    
    
elseif isequal(sign,'down')
    
    for i=1:numel(P)
        % for production
        % cond2-cond1>=0
        for j=1:numel(fcond2{P(i)})
             % add new variables 
            coM = addNewVariableInTFA(coM, strcat(metVar,'_inherP_',num2str(i),'_',num2str(j),'cond1' ), 'C', [0 para.bigM]);
            contV1=length(coM.varNames);
            % nevar=fluxvar*B
            [coM]= add3ConsBinCont(coM,contV1,fcond1{P(i)}(j), useB,para);
            % add new variables 
            coM = addNewVariableInTFA(coM, strcat(metVar,'_inherP_',num2str(i),'_',num2str(j),'cond2' ), 'C', [0 para.bigM]);
            contV2=length(coM.varNames);
            % nevar=fluxvar*B
            [coM]= add3ConsBinCont(coM,contV2,fcond2{P(i)}(j), useB,para);
            % add new variables 
            CLHS.varIDs=[contV2 contV1];
            CLHS.varCoeffs=[1 -1];
            coM=addNewConstraintInTFA(coM, strcat(metVar,'_inherP_',num2str(i),'_',num2str(j) ), '<', CLHS, 0);
        end
    end
    % for consumption
    % cond2-cond1>=0
    for i=1:numel(C)
        for j=1:numel(fcond2{C(i)})
             % add new variables 
            coM = addNewVariableInTFA(coM, strcat(metVar,'_inherC_',num2str(i),'_',num2str(j),'cond1' ), 'C', [0 para.bigM]);
            contV1=length(coM.varNames);
            % nevar=fluxvar*B
            [coM]= add3ConsBinCont(coM,contV1,fcond1{C(i)}(j), useB,para);
            % add new variables 
            coM = addNewVariableInTFA(coM, strcat(metVar,'_inherC_',num2str(i),'_',num2str(j),'cond2' ), 'C', [0 para.bigM]);
            contV2=length(coM.varNames);
            % nevar=fluxvar*B
            [coM]= add3ConsBinCont(coM,contV2,fcond2{C(i)}(j), useB,para);
            % add new variables 
            CLHS.varIDs=[contV2 contV1 fcond2{C(i)}(j) fcond1{C(i)}(j)];
            CLHS.varCoeffs=[-1 1 1 -1];
            coM=addNewConstraintInTFA(coM, strcat(metVar,'_inherC_',num2str(i),'_',num2str(j) ), '>', CLHS, 0);
        end
    end
    
    
    
else
    
    disp('please put sign up and down')
    
end


function [model]=addRelExpCons2Models(model1,model2,regRxns,regRatio)
    
    % model1: in condition1 cobra Or TFA format
    % model2: in condition1 cobra Or TFA format
    % regRxns: reaction for which expression is pertubed
    % regRatio: regulated reaction ratios

    if ~isfield(model1,'var_lb')
        var_ub = [model1.ub;model2.ub];
        var_lb = [model1.lb;model2.lb];
        varNames = [strcat('NF_', model1.rxns);strcat('PERTURB_NF_', model2.rxns)];
        vartypes = [repmat({'C'},numel(model1.rxns),1);repmat({'C'},numel(model2.rxns),1)];
        constraintNames=[strcat('M_', model1.mets);strcat('M_', model2.mets)];
        constraintType=[repmat({'='},numel(model1.mets),1);repmat({'='},numel(model2.mets),1)];
        A=sparse(numel(constraintNames), numel(varNames));
        A(1:numel(model1.mets),1:numel(model1.rxns))=model1.S;
        A(numel(model1.mets)+1:numel(constraintNames),numel(model1.rxns)+1:numel(varNames))=model2.S;
        rhs=[zeros(numel(model1.mets),1);zeros(numel(model2.mets),1)];
        
    else
        % initialize model variables for control and perturbed conditions
        var_ub = [model1.var_ub;model2.var_ub];
        var_lb = [model1.var_lb;model2.var_lb];
        varNames =[model1.varNames;strcat('PERTURB_',model2.varNames)];
        vartypes = [model1.vartypes;model2.vartypes];
        constraintNames=[model1.constraintNames;strcat('PERTURB_',model2.constraintNames)];
        constraintType=[model1.constraintType;model2.constraintType];
        A=sparse(numel(constraintNames), numel(varNames));
        A(1:numel(model1.constraintNames),1:numel(model1.varNames))=model1.A;
        A(numel(model1.constraintNames)+1:numel(constraintNames),numel(model1.varNames)+1:numel(varNames))=model2.A;
        rhs=[model1.rhs;model2.rhs];
    end
    weight=[];
    
    bigM=100000;
    eps=0.00001;
    basal=0.0001;
    objIndex1=[]; % this is for storing results
    objIndex1B=[];
    objIndex2=[];
    objIndex2B=[];
    objBasal=[];
    for i=1:numel(regRxns)
        [~,Find]=ismember(['F_' regRxns{i}],varNames);
        [~,PFind]=ismember(['PERTURB_F_' regRxns{i}],varNames);
        
        if (Find>0) & (PFind>0) 
            %now make free bound for perturbed state
            LB=var_lb(PFind);
            UB=var_ub(PFind);
            [PLB,PUB,B]=makefreeLBUB(LB,UB,1000);
            var_lb(PFind)=PLB;
            var_ub(PFind)=PUB;
            % frist create variabels for Forward
            [num_constr,num_vars] = size(A);
            varNames{length(varNames)+1,1} = ['Yplus_'  'F' regRxns{i}];
            var_ub(length(varNames),1) = bigM;
            var_lb(length(varNames),1) = 0;
            vartypes{length(varNames),1} = 'C';
            bIdx1=length(varNames);
            [num_constr,num_vars] = size(A);
            varNames{length(varNames)+1,1} = ['YBplus_' 'F' regRxns{i}];
            var_ub(length(varNames),1) = 1;
            var_lb(length(varNames),1) = 0;
            vartypes{length(varNames),1} = 'B';
            bIdx1B=length(varNames);

            if regRatio(i)>=1
                
                % F2 >F1*ratio-Yplus
                % F2-F1*raio+Yplus>0
                [num_constr,num_vars] = size(A);
                rhs(num_constr+1,1) = 0;
                constraintNames{num_constr+1,1} = ['Cplus_' 'F' regRxns{i}];
                constraintType{num_constr+1,1} = '>';
                A(num_constr+1,Find) = -regRatio(i);
                A(num_constr+1,PFind) =1;
                A(num_constr+1,bIdx1) = 1; % for yplus
%                 A(num_constr+1,bIdx2) = 1; % for Yconflict
            else
                % F2< F1*ratio+Yplus
                % F2-Yplus-F1*raio<0
                [num_constr,num_vars] = size(A);
                rhs(num_constr+1,1) = 0;
                constraintNames{num_constr+1,1} = ['Cplus_' 'F' regRxns{i}];
                constraintType{num_constr+1,1} = '<';
                A(num_constr+1,Find) = -regRatio(i);
                A(num_constr+1,PFind) =1;
                A(num_constr+1,bIdx1) = -1; % for y plus
%                 A(num_constr+1,bIdx2) = -1; % for y conflict
            end
          
            % eps*(1-YBplus)<=yPlus<=eps+(1-YBplus)*bigM
            % -eps*yBplus-yPlus<=-eps and
            [num_constr,num_vars] = size(A);
            rhs(num_constr+1,1) = -eps;
            constraintNames{num_constr+1,1} = ['C1plus_' 'F' regRxns{i}];
            constraintType{num_constr+1,1} = '<';
            A(num_constr+1,bIdx1B) = -eps; % for YBplus
            A(num_constr+1,bIdx1) = -1; % for Yplus
            %yplus+yBplus*bigM<=eps+bigM
            [num_constr,num_vars] = size(A);
            rhs(num_constr+1,1) = eps+bigM;
            constraintNames{num_constr+1,1} = ['C2plus_' 'F' regRxns{i}];
            constraintType{num_constr+1,1} = '<';
            A(num_constr+1,bIdx1B) = bigM; % for YBplus
            A(num_constr+1,bIdx1) = 1; % for Yplus

            objIndex1=[objIndex1;bIdx1];
%             objIndex2=[objIndex2;bIdx2];
            objIndex1B=[objIndex1B;bIdx1B];
%             objIndex2B=[objIndex2B;bIdx2B];
        end
        
        [~,Rind]=ismember(['R_' regRxns{i}],varNames);
        [~,PRind]=ismember(['PERTURB_R_' regRxns{i}],varNames);
        if (Rind>0) & (PRind>0) 
            %now make free bound for perturbed state
            LB=var_lb(PRind);
            UB=var_ub(PRind);
            [PLB,PUB,B]=makefreeLBUB(LB,UB,1000);
            var_lb(PRind)=PLB;
            var_ub(PRind)=PUB;
            % frist create variabels for backward
            [num_constr,num_vars] = size(A);
            varNames{length(varNames)+1,1} = ['Yplus_'  'R' regRxns{i}];
            var_ub(length(varNames),1) = bigM;
            var_lb(length(varNames),1) = 0;
            vartypes{length(varNames),1} = 'C';
            bIdx1=length(varNames);
            [num_constr,num_vars] = size(A);
            varNames{length(varNames)+1,1} = ['YBplus_' 'R' regRxns{i}];
            var_ub(length(varNames),1) = 1;
            var_lb(length(varNames),1) = 0;
            vartypes{length(varNames),1} = 'B';
            bIdx1B=length(varNames);

            if regRatio(i)>=1
                
                % NF2 > NF1*ratio-Yplus
                % NF2-NF1*raio+Yplus>0
                [num_constr,num_vars] = size(A);
                rhs(num_constr+1,1) = 0;
                constraintNames{num_constr+1,1} = ['Cplus_' 'R' regRxns{i}];
                constraintType{num_constr+1,1} = '>';
                A(num_constr+1,Rind) = -regRatio(i);
                A(num_constr+1,PRind) =1;
                A(num_constr+1,bIdx1) = 1; % for yplus
%                 A(num_constr+1,bIdx2) = 1; % for Yconflict
            else
                % NF2< NF1*ratio+Yplus
                % NF2-Yplus-NF1*raio<0
                [num_constr,num_vars] = size(A);
                rhs(num_constr+1,1) = 0;
                constraintNames{num_constr+1,1} = ['Cplus_' 'R' regRxns{i}];
                constraintType{num_constr+1,1} = '<';
                A(num_constr+1,Rind) = -regRatio(i);
                A(num_constr+1,PRind) =1;
                A(num_constr+1,bIdx1) = -1; % for y plus
%                 A(num_constr+1,bIdx2) = -1; % for y conflict
            end
             % now add five constarint
            % eps*(1-YBplus)<=yPlus<=eps+(1-YBplus)*bigM
            % -eps*yBplus-yPlus<=-eps and
            [num_constr,num_vars] = size(A);
            rhs(num_constr+1,1) = -eps;
            constraintNames{num_constr+1,1} = ['C1plus_' 'R' regRxns{i}];
            constraintType{num_constr+1,1} = '<';
            A(num_constr+1,bIdx1B) = -eps; % for YBplus
            A(num_constr+1,bIdx1) = -1; % for Yplus
            %yplus+yBplus*bigM<=eps+bigM
            [num_constr,num_vars] = size(A);
            rhs(num_constr+1,1) = eps+bigM;
            constraintNames{num_constr+1,1} = ['C2plus_' 'R' regRxns{i}];
            constraintType{num_constr+1,1} = '<';
            A(num_constr+1,bIdx1B) = bigM; % for YBplus
            A(num_constr+1,bIdx1) = 1; % for Yplus 
            
            objIndex1=[objIndex1;bIdx1];

            objIndex1B=[objIndex1B;bIdx1B];

        end
        %% add basal flux in control conditions
        if (Find>0) & (Rind>0) & (PFind>0) & (PRind>0)
            
             %%
            
            [~,plusF]=ismember(['YBplus_' 'F' regRxns{i}],varNames);
            % F1>delta*YBplus_F 
            [num_constr,num_vars] = size(A);
            rhs(num_constr+1,1) = 0;
            constraintNames{num_constr+1,1} = ['basalF1_'  regRxns{i}];
            constraintType{num_constr+1,1} = '>';
            A(num_constr+1,Find) = 1; 
            A(num_constr+1,plusF) = -basal; 
            % F2>delta*YBplus_F 
            [num_constr,num_vars] = size(A);
            rhs(num_constr+1,1) = 0;
            constraintNames{num_constr+1,1} = ['basalF2_'  regRxns{i}];
            constraintType{num_constr+1,1} = '>';
            A(num_constr+1,PFind) = 1; 
            A(num_constr+1,plusF) = -basal; 
            %
            [~,plusR]=ismember(['YBplus_' 'R' regRxns{i}],varNames);
            
            % F1>delta*YBplus_R 
            [num_constr,num_vars] = size(A);
            rhs(num_constr+1,1) = 0;
            constraintNames{num_constr+1,1} = ['basalR1_'  regRxns{i}];
            constraintType{num_constr+1,1} = '>';
            A(num_constr+1,Rind) = 1; 
            A(num_constr+1,plusR) = -basal; 
            % F2>delta*YBplus_R 
            [num_constr,num_vars] = size(A);
            rhs(num_constr+1,1) = 0;
            constraintNames{num_constr+1,1} = ['basalR2_'  regRxns{i}];
            constraintType{num_constr+1,1} = '>';
            A(num_constr+1,PRind) = 1; 
            A(num_constr+1,plusR) = -basal; 
            %
            

        end
        
        
    end
                
                
               %%
    [num_constr,num_vars] = size(A);
    varNames{length(varNames)+1,1} = strcat('Sum_REG');
    var_ub(length(varNames),1) = Inf;
    var_lb(length(varNames),1) = 0;
    vartypes{length(varNames),1} = 'I';
    %% make weighted objective coefficient
    rhs(num_constr+1,1) =  0;
    constraintNames{num_constr+1,1} = strcat('REGSum');
    constraintType{num_constr+1,1} = '=';
    A(num_constr+1,[objIndex1B]) =-1; %-weight; %-1; %% weigted objective coefficient
%     A(num_constr+1,[objBasal]) =-1;
    A(num_constr+1,length(varNames)) = 1;
    %%finished
    
    
    model.A = A;
    model.varNames = varNames;
    model.vartypes = vartypes;
    f=zeros(numel(varNames),1);
    f(end)=1;
    model.var_lb = var_lb;
    model.var_ub = var_ub;
    model.objtype = -1; %-1; % 1 - minimize, -1 - maximize
    model.f = f; % objective vector for TFBA problem
    model.constraintNames = constraintNames;
    model.constraintType = constraintType;
    model.rhs = rhs;
    model.objIndex1=objIndex1;
    model.objIndex2=objIndex2;
    model.objIndex1B=objIndex1B;
    model.objIndex2B=objIndex2B;
    model.objBasal=objBasal;
end


function [LB,UB,B]=makefreeLBUB(LB,UB,bound)
    if LB<0 & UB>0
        LB=-bound;
        UB=bound;
        B=0;
    elseif LB<0 & UB<=0
        LB=-bound;
        UB=0;
        B=0;
    elseif LB>=0 & UB>0
        LB=0;
        UB=bound;
        B=0;
    else
        B=1;
    end
end

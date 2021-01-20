function [model]=addRelSumCons(model,regRxns,regRatio)
     % model: in cobra Or TFA format
    % regRxns: metabolites for which concentartion is pertubed
    % regRatio: regulated reaction ratios
    if ~isfield(model,'var_lb')
         model.var_ub = model.ub;
         model.var_lb = model.lb;
         model.varNames = strcat('NF_', model.rxns);
         model.vartypes = repmat({'C'},numel(model.rxns),1);
         model.constraintNames=strcat('M_', model.mets);
         model.constraintType=repmat({'='},numel(model.mets),1);
         model.A=sparse(model.S);
         model.rhs=zeros(numel(model.mets),1);
    else
        var_ub = model.var_ub;
        var_lb = model.var_lb;
        varNames =model.varNames;
        vartypes =model.vartypes;
        constraintNames=model.constraintNames;
        constraintType=model.constraintType;
        A=model.A;
        rhs=model.rhs;
        
    end
    
    bigM=100000;
    eps=0.00001;
    basal=0.001;
    objIndexM1=[]; % this is for storing results
    objIndexM1B=[];
  
    for i=1:numel(regRxns)
        [~,Find]=ismember(['SumFlux_' regRxns{i}],varNames);
        [~,PFind]=ismember(['Pertb_SumFlux_' regRxns{i}],varNames);
        
        if (Find>0) & (PFind>0) 
            %now make free bound for perturbed state
            LB=var_lb(PFind);
            UB=var_ub(PFind);
            [PLB,PUB,B]=makefreeLBUB(LB,UB,1000);
            var_lb(PFind)=PLB;
            var_ub(PFind)=PUB;
            % frist create variabels for Forward
            [num_constr,num_vars] = size(A);
            varNames{length(varNames)+1,1} = ['Yplus_'   regRxns{i}];
            var_ub(length(varNames),1) = bigM;
            var_lb(length(varNames),1) = 0;
            vartypes{length(varNames),1} = 'C';
            bIdx1=length(varNames);
            [num_constr,num_vars] = size(A);
            varNames{length(varNames)+1,1} = ['YBplus_'  regRxns{i}];
            var_ub(length(varNames),1) = 1;
            var_lb(length(varNames),1) = 0;
            vartypes{length(varNames),1} = 'B';
            bIdx1B=length(varNames);

            if regRatio(i)>=1
                
                % F2 >F1*ratio-Yplus
                % F2-F1*raio+Yplus>0
                [num_constr,num_vars] = size(A);
                rhs(num_constr+1,1) = 0;
                constraintNames{num_constr+1,1} = ['Cplus_'  regRxns{i}];
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
                constraintNames{num_constr+1,1} = ['Cplus_'  regRxns{i}];
                constraintType{num_constr+1,1} = '<';
                A(num_constr+1,Find) = -regRatio(i);
                A(num_constr+1,PFind) =1;
                A(num_constr+1,bIdx1) = -1; % for y plus
%                 A(num_constr+1,bIdx2) = -1; % for y conflict
            end
            % now add five constarint
            % eps*(1-YBplus)<=yPlus<=eps+(1-YBplus)*bigM
            % -eps*yBplus-yPlus<=-eps and
            [num_constr,num_vars] = size(A);
            rhs(num_constr+1,1) = -eps;
            constraintNames{num_constr+1,1} = ['C1plus_'  regRxns{i}];
            constraintType{num_constr+1,1} = '<';
            A(num_constr+1,bIdx1B) = -eps; % for YBplus
            A(num_constr+1,bIdx1) = -1; % for Yplus
            %yplus+yBplus*bigM<=eps+bigM
            [num_constr,num_vars] = size(A);
            rhs(num_constr+1,1) = eps+bigM;
            constraintNames{num_constr+1,1} = ['C2plus_'  regRxns{i}];
            constraintType{num_constr+1,1} = '<';
            A(num_constr+1,bIdx1B) = bigM; % for YBplus
            A(num_constr+1,bIdx1) = 1; % for Yplus

            objIndexM1=[objIndexM1;bIdx1];

            objIndexM1B=[objIndexM1B;bIdx1B];

        end
        
        
        %% add basal flux in control conditions
        if (Find>0)& (PFind>0)
            
             %%
            
            [~,plusF]=ismember(['YBplus_' regRxns{i}],varNames);
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
           
            

        end
        
        
    end
                
                
               %%
    [num_constr,num_vars] = size(A);
    varNames{length(varNames)+1,1} = strcat('Sum_REGEXPMET');
    var_ub(length(varNames),1) = Inf;
    var_lb(length(varNames),1) = 0;
    vartypes{length(varNames),1} = 'I';
    %% make weighted objective coefficient
    rhs(num_constr+1,1) =  0;
    constraintNames{num_constr+1,1} = strcat('Sum_REGEXPMET_C');
    constraintType{num_constr+1,1} = '=';
    A(num_constr+1,[objIndexM1B]) =-1; %-weight; %-1; %% weigted objective coefficient
    A(num_constr+1,[model.objIndex1B]) =-1;
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
    model.objIndexM1=objIndexM1;
    model.objIndexM1B=objIndexM1B;
  
   
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

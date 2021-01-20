function [model]=addSumVars(model,regMets)
    % model in TFA format
    % regMets: cells of regulated metabolites
    TM=[]; % turnover rate for metabolites in cond1 (normal)
    PERT_TM=[]; % turnover rate for metabolites in cond2 (disease)
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
    S=model.S;
    [~,metIdx]=ismember(regMets,model.mets);
    metIdx=metIdx(metIdx>0);
    for i=1:numel(metIdx)
        rId=find(S(metIdx(i),:));
        coeff=abs(S(metIdx(i),rId))*0.5;
     
        % find NF indexes and PERTURB indexes
        [~,nfIdx]=ismember(strcat('NF_',model.rxns(rId)),model.varNames);
  
        [~,pnfIdx]=ismember(strcat('PERTURB_NF_',model.rxns(rId)),model.varNames);
        
        if numel(nfIdx) & numel(pnfIdx)
            bigM=sum(var_ub(nfIdx))/2.0;
            PbigM=sum(var_ub(pnfIdx))/2.0;
            % crate variable for SumFlux
            
            [num_constr,num_vars] = size(A);
            varNames{length(varNames)+1,1} = ['SumFlux_'   model.mets{metIdx(i)}];
            var_ub(length(varNames),1) = bigM;
            var_lb(length(varNames),1) = 0;
            vartypes{length(varNames),1} = 'C';
            bIdx1=length(varNames);
            % crate variable for Perturb SumFlux
            [num_constr,num_vars] = size(A);
            varNames{length(varNames)+1,1} = ['Pertb_SumFlux_'   model.mets{metIdx(i)}];
            var_ub(length(varNames),1) = PbigM;
            var_lb(length(varNames),1) = 0;
            vartypes{length(varNames),1} = 'C';
            bIdx2=length(varNames);
            % add constarint sum for condition 1
            [num_constr,num_vars] = size(A);
            rhs(num_constr+1,1) =  0;
            constraintNames{num_constr+1,1} = ['CSumFlux_'   model.mets{metIdx(i)}];
            constraintType{num_constr+1,1} = '=';
            A(num_constr+1,nfIdx) =-coeff; 
            A(num_constr+1,bIdx1) = 1;
            
            % add constarint sum for condition 2
            [num_constr,num_vars] = size(A);
            rhs(num_constr+1,1) =  0;
            constraintNames{num_constr+1,1} = ['Pertb_CSumFlux_'   model.mets{metIdx(i)}];
            constraintType{num_constr+1,1} = '=';
            A(num_constr+1,pnfIdx) =-coeff; 
            A(num_constr+1,bIdx2) = 1;
            %
            TM=[TM;bIdx1];
            PERT_TM=[PERT_TM;bIdx2];
        end
        
    end
    
    model.A = A;
    model.varNames = varNames;
    model.vartypes = vartypes;
    f=zeros(numel(varNames),1);
    model.var_lb = var_lb;
    model.var_ub = var_ub;
    model.objtype = -1; %-1; % 1 - minimize, -1 - maximize
    model.f = f; % objective vector for TFBA problem
    model.constraintNames = constraintNames;
    model.constraintType = constraintType;
    model.rhs = rhs;
    model.TM=TM;
    model.PERT_TM=PERT_TM;
end
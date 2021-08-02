function [tmpModel]=createMetEnumModel(model,regMets,regRatio,up,down)
    % this function will gives us combined model with condition1 and
    % conditon 2
    %%%Input  model1  is TFA model in condition1 (cond1)
    %         model2 is TFA model in condition2 (cond2)
    %         regGenes  is cellarray of gene names
    %         regRatio is ratio for the rxns (cond2/cond1)
    %         up is upregulation cutoff
    %         down is downregulation cutoff
    
    if nargin<5
        up=.1;
    end
    if nargin<4
        down=0.1;
    end
    
    % add relative Expression into models
    
    indUP=find(regRatio>1);
    ratioUP=regRatio(indUP);
    indDOWN=find(regRatio<1);
    ratioDOWN=regRatio(indDOWN);
    [s_ratioUP,Iup]=sort(ratioUP);
    [s_ratioDOWN,Idown]=sort(ratioDOWN);
    cut1=floor(numel(Iup)*(1-up));
    ind1=indUP(Iup(cut1:end));
    cut2=ceil(numel(Idown)*down);
    ind2=indDOWN(Idown(1:cut2));
    
   
    regMetRatio=[regRatio(ind1);regRatio(ind2)];
    regIdx=[ind1;ind2];
    % find index which are not regulated
    
%     % write function for adding relative metabolites
%     s_model=addSumVars(model,regMets(regIdx));
%     tmpModel=addRelSumCons(s_model,regMets(regIdx),regMetRatio);
    tmpModel=addRelMetabolite(model,model,regMets(regIdx),regMetRatio,true)
 
    tmpModel.description='Regulated Model';
    
end
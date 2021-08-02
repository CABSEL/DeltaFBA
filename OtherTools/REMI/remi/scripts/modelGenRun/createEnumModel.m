function [tmpModel]=createEnumModel(model1,model2,regGenes,regRatio,up,down)
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
    [rxns,ratio]=evaluateGPR(model1,regGenes,regRatio,@geomean,@mean);
    indUP=find(ratio>1);
    ratioUP=ratio(indUP);
    indDOWN=find(ratio<1);
    ratioDOWN=ratio(indDOWN);
    [s_ratioUP,Iup]=sort(ratioUP);
    [s_ratioDOWN,Idown]=sort(ratioDOWN);
    cut1=floor(numel(Iup)*(1-up));
    ind1=indUP(Iup(cut1:end));
    cut2=ceil(numel(Idown)*down);
    ind2=indDOWN(Idown(1:cut2));
    
   
    regRxnRatio=[ratio(ind1);ratio(ind2)];
    regIdx=[ind1;ind2];
    % find index which are not regulated
    tmpId=setdiff([1:numel(model1.rxns)],regIdx);
    
    tmpModel=addRelExpCons2Models(model1,model2,model1.rxns(regIdx),regRxnRatio);
    tmpModel.description='Regulated Model';
    
end
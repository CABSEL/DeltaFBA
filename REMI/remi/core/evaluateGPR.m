function [rxns,rxnExp]=evaluateGPR(model,genes,geneRatio,f_and,f_or)
    % MUmodel is TFA fomat model
    % geneIdx: gene indexes in the model
    % geneRatio: coreesponding gene ratios
    % output
         % rxnIdx= reaction idexes for coreesponding gene indexes
         % ratios=reaction ratios for corresponding ratios
              ratio= geneRatio;%pow2(geneRatio);
              %% crate map from gene index to ration
              geneRatio_map=containers.Map(genes,ratio);
               %[rxnRatio]=geneToreaction_ratios(MUtmodel,geneIdx,ratio,@geomean,@mean)
              [GPRrules,store] = parsingGPRs(model);
        %     load('/Users/vikashpandey/Documents/MATLAB/MetGene/parsedGPR_MUmodel.mat')
              [rxns,rxnExp]=getRxnvalues(model,GPRrules,geneRatio_map,f_and,f_or);
        
end


function [rxns,rxExp]=getRxnvalues(model,parseGPRs,geneRatio_map,f_and,f_or)
    % model
    % parseGPRs: reactions genes GPR
    % geneRatio_map : geneRatio_map
    rxns=parseGPRs.keys;
    rxExp=nan(numel(rxns),1);
    parseRules=parseGPRs.values;
    
    for r=1:numel(rxns)
       
        rxExp(r)=eval_eachRxn(model,rxns{r},parseRules{r},geneRatio_map,f_and,f_or);
    end
end

function  final_val=eval_eachRxn(model,rxn,ruleArray,geneRatio,f_and,f_or)
    % ruleArray is number of cells that are parsed with 'or' rules
    % and inside each cell there is another cell array than contains genes
    % which are seperated by and rule.
   
    
    if numel(ruleArray)==1 & numel(ruleArray{1})==1
        % this is for those genes which are catalized by one gene
        if isKey(geneRatio,ruleArray{1}{1})
            final_val=geneRatio(ruleArray{1}{1});
        else
            final_val=NaN;
        end
       
    elseif numel(ruleArray)==1 & numel(ruleArray{1})>1 %% this is the case when they are catlized by and only
        
        andVal=[];
        for j=1:numel(ruleArray{1})
            % get value and apply operation
            if isKey(geneRatio,ruleArray{1}{j})
                andVal=[andVal;geneRatio(ruleArray{1}{j})];
            end
            
        end
        % this script is for ignoring nan values
        andVal=andVal(~isnan(andVal));
        if ~isempty(andVal)
            final_val=f_and(andVal);
            if numel(find(isnan(andVal)))>0
                disp('please check andVal contains nan (can not use geomean)')
            end
        else
            final_val=NaN;
        end
    else
        % this is for those catlized by or as well as and
        onVal=[];
        for i=1:numel(ruleArray)
            andVal=[];
            for j=1:numel(ruleArray{i})
                %%get value and apply operation
                if isKey(geneRatio,ruleArray{i}{j})
                    andVal=[andVal;geneRatio(ruleArray{i}{j})];
                end
            end
            % this script is for ignoring nan values
            andVal=andVal(~isnan(andVal));
            if ~isempty(andVal)
                % apply function of andRule
                
                onVal=[onVal;f_and(andVal)];
                if numel(find(isnan(andVal)))>0
                    disp('please check andVal contains nan (can not use geomean)')
                end
            end
        end
         % this script is for ignoring nan values
        onVal=onVal(~isnan(onVal));
        if ~isempty(onVal)
            
            
            final_val=f_or(onVal);
        else
            final_val=NaN;
        end
        
    end
end


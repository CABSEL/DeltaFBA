function reaction_levels = gene_to_reaction_levels( model, genes, levels, f_and, f_or )
% Convert gene expression levels to reaction levels using GPR associations.
% Level is NaN if there is no GPR for the reaction or no measured genes.
%
% INPUTS
%       model - cobra model
%       genes - gene names
%       levels - gene expression levels
%       f_and - function to replace AND
%       f_or - function to replace OR
%
% OUTPUTS
%       reaction_levels - reaction expression levels
%
% Author: Daniel Machado, 2013

    reaction_levels = zeros(length(model.rxns), 1);

    for i = 1:length(model.rxns)
        level = eval_gpr(model.grRules{i}, genes, levels, f_and, f_or);
        reaction_levels(i) = level;
    end

end

function [result, status] = eval_gpr(rule, genes, levels, f_and, f_or)
% Evaluate the expression level for a single reaction using the GPRs.
% Note: Computes the expression level even if there are missing measured
% values for the given rule. This implementation is a modified version of
% an implementation provided in [Lee et al, BMC Sys Biol, 2012]

    EVAL_OK = 1;
    PARTIAL_MEASUREMENTS = 0;
    NO_GPR_ERROR = -1;
    NO_MEASUREMENTS = -2;
    MAX_EVALS_EXCEEDED = -3;

    MAX_EVALS = 1000;
    NONETYPE = 'NaN';

    NUMBER = '[0-9\.\-e]+';
    MAYBE_NUMBER = [NUMBER '|' NONETYPE];

    expression = rule;
    result = NaN;
    status = EVAL_OK;

    if isempty(expression)
        status = NO_GPR_ERROR;
    else
        rule_genes = setdiff(regexp(expression,'\<(\w|\-)+\>','match'), {'and', 'or'});
        
        total_measured = 0;
        
        for i = 1:length(rule_genes)
            j = find(strcmp(rule_genes{i}, genes));
            if isempty(j)
                level = NONETYPE;
            else
                level = num2str(levels(j));
                total_measured = total_measured + 1;
            end
            expression = regexprep(expression, ['\<', rule_genes{i}, '\>'], level );
        end
        
        
        if total_measured == 0
            status = NO_MEASUREMENTS;
        else
            if total_measured < length(rule_genes)
                status = PARTIAL_MEASUREMENTS;
            end
            
            maybe_and = @(a,b)maybe_functor(f_and, a, b);
            maybe_or = @(a,b)maybe_functor(f_or, a, b); 
            str_wrapper = @(f, a, b)num2str(f(str2double(a), str2double(b)));

            counter = 0;
            
            while isnan(result)

                counter = counter + 1;
                if counter > MAX_EVALS
                    status = MAX_EVALS_EXCEEDED;
                    break
                end

                try 
                    result = eval(expression);            
                catch e   
                    paren_expr = ['\(\s*(', MAYBE_NUMBER,')\s*\)'];
                    and_expr = ['(',MAYBE_NUMBER,')\s+and\s+(',MAYBE_NUMBER,')'];
                    or_expr = ['(',MAYBE_NUMBER,')\s+or\s+(',MAYBE_NUMBER,')'];

                    expression = regexprep(expression, paren_expr, '$1');
                    expression = regexprep(expression, and_expr, '${str_wrapper(maybe_and, $1, $2)}');
                    expression = regexprep(expression, or_expr, '${str_wrapper(maybe_or, $1, $2)}');
                end
            end
            
        end
    end

end


function c = maybe_functor(f, a, b)
    
    if isnan(a) && isnan(b)
        c = nan;
    elseif ~isnan(a) && isnan(b)
        c = a;
    elseif isnan(a) && ~isnan(b)
        c = b;
    else 
        c = f(a,b);
    end
end

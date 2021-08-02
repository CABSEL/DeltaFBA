function fluxes = call_MADE(model, gene_names, gene_exp, gene_exp_ref, obj_frac, elf_model, bounds_ref )
% Calls the implementation of MADE provided in [Jensen and Papin, Bioinformatics, 2011].
% Please install the code directly from: http://www.bme.virginia.edu/csbl/downloads-made.php
% 
% Note: This call wrapper only considers two experimental conditions.
%
% INPUTS
%       model - cobra model
%       gene_names - genes ids
%       gene_exp - gene expression
%       gene_exp_ref - gene expression for reference condition
%       obj_frac - minimum fraction of the biological objective
%       elf_model - ELF model (optional)
%       bounds_ref - flux bounds for reference condition
%
% OUTPUTS
%       fluxes - flux distribution
%
% Author: Daniel Machado, 2013 
    
    fold_change = gene_exp ./ gene_exp_ref;
    pvals = zeros(size(fold_change));
    
    if isempty(elf_model)
        elf_model = assert_elf_model(model);
    end
    
    %lb and ub must be same size of elf model
    bounds.lb = [model.lb; elf_model.lb(length(model.rxns)+1:end)];
    bounds.ub = [model.ub; elf_model.ub(length(model.rxns)+1:end)];
    bounds_ref.lb = [bounds_ref.lb; elf_model.lb(length(model.rxns)+1:end)];
    bounds_ref.ub = [bounds_ref.ub; elf_model.ub(length(model.rxns)+1:end)];
    
    [~,~,sol,~] = made(elf_model, fold_change, pvals, obj_frac, ...
        'gene_names', gene_names, 'weighting', 'unit', 'bounds', {bounds_ref, bounds});
    
    if isempty(sol.x)
        fluxes = [];
    else
        fluxes = sol.x(length(elf_model.c)+1:length(elf_model.c)+length(model.rxns));
    end
                           
end


function fluxes = call_EFlux(model, gene_names, gene_exp, scale_rxn, scale_value)
% Implements E-Flux as defined in [Colijn et al, PLoS Comp Bio, 2009].
% 
% Note: E-Flux ignores all previously defined flux bounds.
%
% INPUTS
%       model - cobra model
%       gene_names - gene ids
%       gene_exp - gene expression
%       scale_rxn - reaction to scale
%       scale_value - flux of scale reaction
%
% OUTPUTS
%       fluxes - flux distribution
%
% Author: Daniel Machado, 2013 
    
    levels = gene_to_reaction_levels(model, gene_names, gene_exp, @min, @(x,y)(x+y));
    levels = levels / max(levels);
    levels(isnan(levels)) = 1;
    
    blocked_lb = model.lb >= 0;
    blocked_ub = model.ub <= 0;

    model.lb = -levels;
    model.lb(blocked_lb) = 0;
    model.ub = levels;
    model.ub(blocked_ub) = 0;

    sol = optimizeCbModel(model);
    
    if isempty(sol.x)
        fluxes = [];
    else
        scale_idx = strcmp(scale_rxn, model.rxns);
        if sol.x(scale_idx) == 0
            fluxes = [];
        else
            fluxes = sol.x * abs(scale_value / sol.x(scale_idx));
        end
    end
end

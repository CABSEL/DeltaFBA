function [model, fitted_rates] = constrain_model_to_data(model, reactions, rates)
% Constrain model to given flux rates.
%
% Author: Daniel Machado, 2013

    TOL = 1e-9;
    fitted_rates = fit_fluxes_to_model(model, reactions, rates);
    [~, reaction_idxs] = ismember(reactions, model.rxns);
    model.lb(reaction_idxs) = fitted_rates - TOL * abs(fitted_rates);
    model.ub(reaction_idxs) = fitted_rates + TOL * abs(fitted_rates);
end
function [ v_new ] = remove_futile_cycles(model, v)
% Remove futile cyles from a given flux distribution
%
% INPUTS
%       model - cobra model
%       v - flux distribution
%
% OUTPUTS
%       v_new - new flux distribution
%
% Author: Daniel Machado, 2013

    MAX_FLUX = 100; %potential futile cycle threshold
    TOL = 1e-6;
    
    valid_fluxes = find(abs(v) <= MAX_FLUX);
    
    if length(valid_fluxes) < length(v)

        model.lb(valid_fluxes) = v(valid_fluxes) - TOL;
        model.ub(valid_fluxes) = v(valid_fluxes) + TOL;

        sol = optimizeCbModel(model, 'max', 'one');
        
        if ~isempty(sol.x)
            v_new = sol.x;
        else
            v_new = v;
            disp('Failed to remove futile cycles.');
        end
    else
        v_new = v;
    end
end


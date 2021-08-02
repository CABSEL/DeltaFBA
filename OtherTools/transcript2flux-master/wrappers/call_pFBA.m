function fluxes = call_pFBA(model, knockouts)
% Call pFBA
%
% INPUTS
%       model - cobra model
%       knockouts - genes to delete (if any)
%
% OUTPUTS
%       fluxes - flux distribution
%
% Author: Daniel Machado, 2013    

    if ~isempty(knockouts)
        model = deleteModelGenes(model, knockouts);
    end
    
    result = optimizeCbModel(model, 'max', 'one');
    fluxes = result.x;    
end


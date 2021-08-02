function data = precompute_GXFBA( model_ref )
% Pre-compute flux variability for GX-FBA
%
% Author: Daniel Machado, 2013 

    data.wt_sol = optimizeCbModel(model_ref,'max');
    [wt_minf, wt_maxf] = fluxVariability(model_ref,0,'max');
    data.wt_minf = wt_minf;
    data.wt_maxf = wt_maxf;

end


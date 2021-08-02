function [fluxes, status, runtime] = call_method(method, args, clean_futile_cycles)
% Generic wrapper for calling any method
%
% INPUTS
%       method - method name
%       args - structure with all method specific arguments
%       clean_futile_cycles - remove futile cycles after computation (bool)
%
% OUTPUTS
%       fluxes - flux distribution
%       status - exit status
%       runtime - runtime of method (seconds)
%
% Author: Daniel Machado, 2013

    tstart = tic();
    
    switch method
        case 'pFBA'
            fluxes = call_pFBA(args.model, args.FBA.knockouts);
            
        case 'GIMME'
            fluxes = call_GIMME(args.model, args.gene_names, args.gene_exp, ...
                args.obj_frac, args.gimme.threshold);
                            
        case 'iMAT'
            fluxes = call_iMAT(args.model, args.gene_names, args.gene_exp, ...
                args.iMAT.low_threshold, args.iMAT.up_threshold, args.iMAT.eps);
            
        case 'MADE'
            fluxes = call_MADE(args.model, args.gene_names, args.gene_exp, ...
                args.gene_exp_ref, args.obj_frac, args.MADE.elf_model, args.MADE.bounds_ref);
            
        case 'E-Flux'
            fluxes = call_EFlux(args.model, args.gene_names, args.gene_exp, ...
                args.flux_scale_rxn, args.scale_value);

        case 'Lee-12'
            fluxes = call_Lee12(args.model, args.gene_names, args.gene_exp, ...
                args.gene_exp_sd, args.gene_scale_rxn, args.flux_scale_rxn, args.scale_value);
            
        case 'RELATCH'
            fluxes = call_RELATCH(args.model, args.gene_names, args.gene_exp, ...
                args.external_rxns, args.external_rates);
            
        case 'GX-FBA'
            fluxes = call_GXFBA(args.model, args.gene_names, args.gene_exp, ...
                args.gene_exp_ref, args.GXFBA);
            
        otherwise
            fluxes = [];
            fprintf(1, 'Unknown method: %s\n', method);
    end
    
    runtime = toc(tstart);
    status = ~isempty(fluxes);
    
    if status && clean_futile_cycles
        fluxes = remove_futile_cycles(args.model, fluxes);
    end
end


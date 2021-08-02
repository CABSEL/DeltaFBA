function result = evaluate_method(model, dataset, method, condition, ref_condition, options)
% Evaluate a method for a specific condition from a given dataset.
%
% INPUTS
%       model - cobra model
%       dataset - dataset structure (created by load_dataset)
%       method - method name
%       condition - condition id
%       ref_condition - reference condition id
%       options - simulation options
%
% OUTPUTS
%       result - evaluation result
%
% Author: Daniel Machado, 2013
    
    OBJ_FRAC = 0.9;
    GIMME_LOWER_QUANTILE = 0.25;
    IMAT_LOWER_QUANTILE = 0.25;
    IMAT_UPPER_QUANTILE = 0.75;
    IMAT_EPS = 1;
    
    if isfield(options, 'eval_str')
        eval(options.eval_str)
    end
    
    switch options.experiment_type
        case 'sim_all'
            options.set_growth_rate = false;
            options.set_secretion_rates = false;
            options.ignore_internal_fluxes = false;
        case 'sim_intra'
            options.set_growth_rate = true;
            options.set_secretion_rates = true;
            options.ignore_internal_fluxes = false;
        case 'sim_secr'
            options.set_growth_rate = false;
            options.set_secretion_rates = false;
            options.ignore_internal_fluxes = true;
    end
    
    cond_idx = strcmp(condition, dataset.conditions);
    ref_cond_idx = strcmp(ref_condition, dataset.conditions);

    gene_exp = dataset.transcriptomics(:,cond_idx);
    gene_exp_ref = dataset.transcriptomics(:,ref_cond_idx);
    
    fluxes_exp = dataset.fluxomics(:,cond_idx);
    fluxes_exp_ref = dataset.fluxomics(:,ref_cond_idx);
        
    simulated = dataset.reactions;
    
    model_ref = model;
    
    biomass = find(model.c);
    biomass_rxn = model.rxns(biomass);

    if options.set_growth_rate
        growth_idx = find(strcmp(biomass_rxn, dataset.reactions));
        if isempty(growth_idx)
            disp('Unable to set growth rate.')
        else
            model.lb(biomass) = fluxes_exp(growth_idx);
            model.ub(biomass) = fluxes_exp(growth_idx);
            model_ref.lb(biomass) = fluxes_exp_ref(growth_idx);
            model_ref.ub(biomass) = fluxes_exp_ref(growth_idx);
            simulated = setdiff(simulated, biomass_rxn);
            OBJ_FRAC = 1;
        end
    end
    
    fluxes_exp_old = fluxes_exp;
            
    if options.reestimate_data
        fluxes_exp = fit_fluxes_to_model(model, dataset.reactions, fluxes_exp_old);
        fluxes_exp_ref = fit_fluxes_to_model(model_ref, dataset.reactions, fluxes_exp_ref);
    end
    
    if options.set_secretion_rates
        exchange_rxns = intersect(dataset.reactions, [model.uptk_rxns; model.secr_rxns]);
    else
        exchange_rxns = intersect(dataset.reactions, model.uptk_rxns);
    end

    [~, exchange_idxs_data] = ismember(exchange_rxns, dataset.reactions);
    exchange_rates = fluxes_exp(exchange_idxs_data);
    model  = constrain_model_to_data(model, exchange_rxns, exchange_rates);
    exchange_rates_ref = fluxes_exp_ref(exchange_idxs_data);
    model_ref = constrain_model_to_data(model_ref, exchange_rxns, exchange_rates_ref);
    
    if options.ignore_internal_fluxes
        simulated = intersect(dataset.reactions, [biomass_rxn; model.secr_rxns]);
    else
        simulated = setdiff(simulated, exchange_rxns);
    end

    [~, simulated_idx_data] = ismember(simulated, dataset.reactions);
    [~, simulated_idx_model] = ismember(simulated, model.rxns);
    [~, exchange_idxs_model] = ismember(exchange_rxns, model.rxns);
    
    fluxes_exp_sim = fluxes_exp(simulated_idx_data);

    uptake_rate = abs(fluxes_exp(strcmp(model.subtrate_uptake_rxn, dataset.reactions)));

    args.model = model;
    args.model_ref = model_ref;
    args.gene_names = dataset.genes;
    args.gene_exp = gene_exp;
    args.gene_exp_sd = [];
    args.gene_exp_ref = gene_exp_ref;
    args.external_rxns = exchange_rxns;
    args.external_rates = exchange_rates;
    args.external_rates_ref = exchange_rates_ref;
    args.obj_frac = OBJ_FRAC;
    args.gene_scale_rxn = model.subtrate_transport_rxn;
    args.flux_scale_rxn = model.subtrate_uptake_rxn;
    args.scale_value = uptake_rate;
    args.gimme.threshold = quantile(gene_exp, GIMME_LOWER_QUANTILE);
    args.FBA.knockouts = dataset.knockouts{cond_idx};
    args.iMAT.low_threshold = quantile(gene_exp, IMAT_LOWER_QUANTILE);
    args.iMAT.up_threshold = quantile(gene_exp, IMAT_UPPER_QUANTILE);
    args.iMAT.eps = IMAT_EPS;
    
    if strcmp('MADE', method)
        if isfield(options.precomputed,'MADE')
            args.MADE = options.precomputed.MADE;
        else
            args.MADE = precompute_MADE(model_ref);
            result.precomputed.MADE = args.MADE;
        end
    end
    
    if strcmp('GX-FBA', method)
        if isfield(options.precomputed,'GXFBA')
            args.GXFBA = options.precomputed.GXFBA;
        else
            args.GXFBA = precompute_GXFBA(model_ref);
            result.precomputed.GXFBA = args.GXFBA;
        end
    end
    
    if isfield(dataset, 'transcriptomics_std')
        args.gene_exp_sd = dataset.transcriptomics_std(:,cond_idx);
    end
        
    [fluxes, status, runtime] = call_method(method, args, true);
    
    result.status = status;
    result.simulated = simulated;
    result.fluxes = fluxes;
    result.fluxes_exp = fluxes_exp;
    result.runtime = runtime;
    result.fluxes_exp_sim = fluxes_exp_sim;
    result.exchange_rxns = exchange_rxns;

    if status
        result.fluxes_sim = fluxes(simulated_idx_model);
        result.error = norm(result.fluxes_sim - fluxes_exp_sim) / norm(fluxes_exp_sim);
        result.error_exchange = norm(fluxes(exchange_idxs_model) - exchange_rates);
    else
        result.fluxes_sim = [];
        result.error = NaN;
        result.error_exchange = NaN;
    end
            
end

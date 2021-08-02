function experiment = benchmark_method(method, model, dataset, options)
% Test a method using all conditions from a given dataset.
%
% INPUTS
%       method - method name
%       model - cobra model
%       dataset - dataset structure (created by load_dataset)
%       options - simulation options
%
% OUTPUTS
%       experiment - experiment results
%
% Author: Daniel Machado, 2013

    N = length(dataset.conditions);
    
    status_all = zeros(1,N);
    error_all = zeros(1,N);
    runtime_all = zeros(1,N);
    fluxes_exp_all = cell(1,N);
    fluxes_sim_all = cell(1,N);
    
    options.precomputed = [];
    
    h = waitbar(0, sprintf('testing %s for %s dataset\n', method, dataset.name));
    set(findall(h,'type','text'),'Interpreter','none');

    for i=1:N
        condition = dataset.conditions{i};
        ref_condition = dataset.conditions{1};
        result = evaluate_method(model, dataset, method, condition, ref_condition, options);
        
        status_all(i) = result.status;
        if result.status
            fluxes_sim_all{i} = result.fluxes_sim;
            fluxes_exp_all{i} = result.fluxes_exp_sim;
            error_all(i) = result.error; 
            runtime_all(i) = result.runtime;
        end
        
        if isfield(result, 'precomputed')
            options.precomputed = result.precomputed;
        end
        waitbar(i/N);
    end
    
    experiment.name = sprintf('%s_%s_%s', method, dataset.name, options.experiment_type);
    experiment.conditions = dataset.conditions;
    experiment.status_all = status_all;
    experiment.error_all = error_all;
    experiment.runtime_all = runtime_all;
    experiment.sim_reactions = result.simulated;
    experiment.fluxes_exp_all = fluxes_exp_all;
    experiment.fluxes_sim_all = fluxes_sim_all;
    
    save(['results/' experiment.name '.mat'], 'experiment');
    close(h);
end


function analysis = sensitivity_analysis(model, dataset, method, parameter, min_val, max_val, points, scale, options)
% Analyse the sensitivity of a method to its parameters.
%
% INPUTS
%       model - cobra model
%       dataset - dataset structure (created by load_dataset)
%       method - method name
%       parameter - parameter id
%       min_val - minimum value
%       max_val - maximum value
%       points - number of sampling points
%       scale - 'lin' or 'log'
%       options - simulation options
%
% OUTPUTS
%       analysis - analysis result
%
% Author: Daniel Machado, 2013
        
    N = length(dataset.conditions);
    
    switch scale
        case 'lin'
            values = linspace(min_val, max_val, points);
        case 'log'
            values = logspace(min_val, max_val, points);
    end

    error_data = cell(points,1);    
    options.precomputed = [];
        
    h = waitbar(0, sprintf('sensitivity analysis of %s to %s for %s\n', method, parameter, dataset.name));
    set(findall(h,'type','text'),'Interpreter','none');
    
    for j = 1:points
        options.eval_str = sprintf('%s = %f;', parameter, values(j));
    
        status_all = zeros(1,N);
        error_all = zeros(1,N);

        for i=1:N
                        
            condition = dataset.conditions{i};
            ref_condition = dataset.conditions{1};
            result = evaluate_method(model, dataset, method, condition, ref_condition, options);

            status_all(i) = result.status;
            error_all(i) = result.error;
            
            if isfield(result, 'precomputed')
                options.precomputed = result.precomputed;
            end
            
            waitbar(((j-1)*N + i) / (points*N));
        end
        
        error_data{j} = error_all(status_all > 0);
    end
    
    analysis.parameter = parameter;
    analysis.method = method;
    analysis.param_values = values;
    analysis.scale = scale;
    analysis.error_data = error_data;
    save(sprintf('results/sensitivity_%s_%s_%s.mat', method, parameter, dataset.name), 'analysis');
    close(h)

end


function analysis = robustness_analysis(model, dataset, method, condition, ref_condition, steps, points, options)
% Analyse the robustness of a method to noisy data.
%
% INPUTS
%       model - cobra model
%       dataset - dataset structure (created by load_dataset)
%       method - method name
%       condition - condition id
%       ref_condition - reference condition id
%       steps - number of steps
%       points - number of replicates per step
%       options - simulation options
%
% OUTPUTS
%       analysis - analysis result
%
% Author: Daniel Machado, 2013

    options.precomputed = [];
        
    alpha = linspace(0,1,steps);
    cond_idx = strcmp(condition, dataset.conditions);
    
    measured = dataset.transcriptomics(:,cond_idx);
    if isfield(dataset, 'transcriptomics_std')
        measured_std = dataset.transcriptomics_std(:,cond_idx);
    end
    
    error_data = cell(steps,1);
    
    h = waitbar(0, sprintf('robustness analysis of %s for %s (%s)\n', method, dataset.name, condition));
    set(findall(h,'type','text'),'Interpreter','none');
    
    for i = 1:steps
        
        status_all = zeros(1,points);
        error_all = zeros(1,points);
        
        for j = 1:points
                        
            permutation = randperm(length(measured));
            scrambled = measured(permutation);
            noisy = measured + alpha(i) * (scrambled - measured);
            dataset.transcriptomics(:,cond_idx) = noisy;

            if isfield(dataset, 'transcriptomics_std')
                scrambled_std = measured_std(permutation);
                noisy_std = measured_std + alpha(i) * (scrambled_std - measured_std);
                dataset.transcriptomics_std(:,cond_idx) = noisy_std;
            end
            
            result = evaluate_method(model, dataset, method, condition, ref_condition, options);
            status_all(j) = result.status;
            error_all(j) = result.error;
            
            if isfield(result, 'precomputed')
                options.precomputed = result.precomputed;
            end
            
            if alpha(i) == 0
                break
            end
            
            waitbar(((i-1)*points + j) / (points*steps));
        end
        
        error_data{i} = error_all(status_all > 0);
    end
        
    analysis.method = method;
    analysis.dataset = dataset.name;
    analysis.condition = condition;
    analysis.alpha = alpha;
    analysis.error_data = error_data;

    save(sprintf('results/robustness_%s_%s_%s.mat', method, dataset.name, condition), 'analysis');
    close(h)

end


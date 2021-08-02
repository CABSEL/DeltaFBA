function heatmap(experiment, show_labels, clim)
% Build heat map.
%
% Author: Daniel Machado, 2013

    n_rows = length(experiment.sim_reactions);
    n_cols = length(experiment.conditions);
    data = nan*ones(n_rows, n_cols);
    
    for i = 1:n_cols
        if experiment.status_all(i)
            data(:,i) = experiment.fluxes_sim_all{i} - experiment.fluxes_exp_all{i};
        end
    end
    
    values = reshape(data,1,[]);
    values = values(~isnan(values));
    std_data = std(values);
    
    if isempty(clim)
        clim = 3*[-std_data std_data];
    end
    
    h = imagesc(data, clim);
    set(h, 'AlphaData', ~isnan(data));

    colormap(blue_yellow())
    colorbar('SouthOutside')
    set(gca, 'XTick', []);
    if show_labels
        set(gca, 'YTick', 1:n_rows, 'YTickLabel', experiment.sim_reactions, 'FontWeight', 'bold');
    else
        set(gca, 'YTick', 1:n_rows, 'YTickLabel', []);
    end
    set(gca, 'LineWidth', 1);
    
end

function map = blue_yellow()
 
    N = 64;
    X = 0.5: -1/(N-1):-0.5;
    X = abs(X).*2;

    R = [zeros(1,N/2) X(:,(N/2 + 1):N)]';
    B = [X(:,1:N/2) zeros(1,N/2)]';
    G = [zeros(1,N/2) X(:,(N/2 + 1):N)]';

    map = [R G B];
end
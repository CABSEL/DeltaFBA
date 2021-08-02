function plot_error( error )
% Plot scaled error levels.
%
% Author: Daniel Machado, 2013
    
    n = length(error);
    plot(error, '.', 'MarkerSize', 12);
    xlim([0.25 n+0.75])
    min_y = min(error); max_y = max(error);
    dif = max_y - min_y;
    if dif == 0
        dif = max_y;
    end
    ylim([min_y max_y] + [-0.05 0.05]*dif);
    xlabel('experiment', 'FontWeight', 'bold')
    ylabel('error', 'FontWeight', 'bold')
    set( gca, 'LineWidth', 1, 'FontWeight', 'bold' )
    set( gca, 'Xtick', [] )
    set( gca, 'Ytick', [] )
    box on
    camva('manual')
end


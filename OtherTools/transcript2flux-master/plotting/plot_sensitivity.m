function plot_sensitivity(values, data, scale, ymax, x_label, y_label )
% Plot sensitivity/robustness analysis.
%
% Author: Daniel Machado, 2013
    
    n = length(data);
    
    data_avg = zeros(1, n);
    data_std = zeros(1, n);
    
    for i = 1:n
        data_avg(i) = mean(data{i});
        data_std(i) = std(data{i});
    end

    errorbar(values, data_avg, data_std, '.-', 'LineWidth', 1, 'MarkerSize', 12);
    xlabel(x_label);
    
    
    switch scale
        case 'lin'
            set(gca,'xscale','lin')
            xlim([min(values) max(values)] + [-0.05 0.05])
        case 'log'
            set(gca,'xscale','log')
            xlim([min(values) max(values)] .* [0.5 2])
    end
    
    if ymax
        ylim([0 ymax])
        set(gca, 'YTick', [0 ymax])
        ylabel(y_label);
    else
        y_lim = [min(data_avg - data_std), max(data_avg + data_std)];
        y_lim = y_lim + [-1 1]*0.05*norm(y_lim(2) - y_lim(1));
        ylim(y_lim)
        set(gca, 'YTick', [])        
    end
    
    set(gca, 'XTick', [min(values) max(values)])    
    set( gca, 'LineWidth', 1, 'FontWeight', 'bold' )

end


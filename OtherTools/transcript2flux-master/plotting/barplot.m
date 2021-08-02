function barplot(values, reference, plot_title, labels)
% Build bar plot.
%
% Author: Daniel Machado, 2013

    n = length(values);
    hold on
    box on
    bar(1, reference, 'r')
    bar(2:n+1, values)
    
    set(gca, 'XTick', 1:n+1, 'XTickLabel', ['Measured' labels]);
    set(gca, 'LineWidth', 1, 'FontWeight', 'bold')
    
    xlim([0.5, n+1.5]);

    y_lim = [min([values 0 reference]), max([values 0 reference 1e-1])];
    y_lim = y_lim + [-0.1 0.1]*norm(y_lim);
    if norm(y_lim) > 0
        ylim(y_lim);
    end
    
    title(plot_title);
    rotateXLabels(gca, 90);
    hold off
end


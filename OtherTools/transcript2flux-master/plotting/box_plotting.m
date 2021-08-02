function box_plotting( data, labels , ymax)
% Build box plot.
%
% Author: Daniel Machado, 2013

    values = [];
    groups = [];
    max_points = 0;
    
    for i = 1:length(data)
        n = length(data{i});
        max_points = max(max_points, n);
        values = [values data{i}];
        groups = [groups repmat(i, 1, n)];
    end
    
    h = boxplot(values, groups, 'labels', labels, 'labelorientation', 'inline', 'symbol', '.');
    
    set( h, 'LineWidth', 1 )
    ylim([0 ymax])
    set( gca, 'LineWidth', 1, 'FontWeight', 'bold' )
    set(findobj(gca,'Type','text'),'FontWeight', 'bold')
end


[fluxId]=textread('/Users/vikash/Desktop/REMI/data/fluxid.txt','%s');

% We are making scatter plot for reviewer 

exps={ 'pgm', 'pgi', 'gapC', 'zwf', 'rpe', 'wt5', 'wt7', 'NOX',	'ATPase'};
%% TGexM

folderIn='/Users/vikash/Desktop/REMI/simData/ModelsSolutions/TGexM/store_';
for i=1 %1:numel(exps)
    try
        load([folderIn exps{i} '.mat'])
        XX=storeXX(:,:,1);
        X=XX(:,1);
        Y=XX(:,2);
        plot(X,Y,'o','MarkerSize',10,...
        'MarkerEdgeColor','k')
%         dx = 0.1; dy = 0.1;
%         text(X+dx, Y+dy, fluxId, 'Fontsize', 10);
%         scatter(X,Y,50,'k','o');
         xlabel('Experimental flux (unit: mmol/g DCW/h)')
         ylabel('Experimental flux (unit: mmol/g DCW/h)')
    end
end

data=[fluxId num2cell([X Y])];

% after arranging the data 

[~,~,data]=xlsread('/Users/vikash/Desktop/REMI/simData/scatter_plot/TGexM.xlsx','zwf');
ll=data(:,1);
% please load manually from excel file
xx=cell2mat(x);
X=xx(:,1);
Y=xx(:,2);
plot(X,Y,'o','MarkerSize',12,...
    'MarkerEdgeColor','k')
dx = 0.03; dy = 0.03;
text(X+dx, Y+dy, ll, 'Fontsize', 12);

xlabel('Experimental flux (mmol/g DW/h)')
ylabel('Predicted flux (mmol/g DW/h)')
title ('zwf')

set(gca, 'Fontsize', 15)



%%

% after arranging the data 

[~,~,data]=xlsread('/Users/vikash/Desktop/REMI/simData/scatter_plot/TGexM.xlsx','pgm');
ll=data(:,1);
% please load manually from excel file
xx=cell2mat(x);
X=xx(:,1);
Y=xx(:,2);
plot(X,Y,'o','MarkerSize',12,...
    'MarkerEdgeColor','k')
dx = 0.03; dy = 0.03;
text(X+dx, Y+dy, ll, 'Fontsize', 12);

xlabel('Experimental flux (mmol/g DW/h)')
ylabel('Predicted flux (mmol/g DW/h)')
title ('pgm')

set(gca, 'Fontsize', 15)





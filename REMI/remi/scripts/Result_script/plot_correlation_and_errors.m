% With this script we are plotting correlation between experiment flux and
% simulated flux. We did it for for all REMI models and GX-FBA 

%% read REMI-MET
exps={ 'pgm', 'pgi', 'gapC', 'zwf', 'rpe', 'wt5', 'wt7', 'NOX',	'ATPase'}; % this is the type of mutants
storeGXFBA1=nan(numel(exps),4); % this to store mean correlations and percentage error
% load the results which are generated and saved in the folder :
% ModelsSolutions
load('/Users/vikash/Desktop/REMI/simData/ModelsSolutions/GXFBA/storeBIshiiNew_gxfba.mat') % Ishii data set
k=1;
for i=1:numel(storeB)
    storeGXFBA1(k,:)=storeB{i};
    k=k+1;
end
% clear storeB

load('/Users/vikash/Desktop/REMI/simData/ModelsSolutions/GXFBA/storeBHolm_gxfba.mat') % holm data  data set

for i=1:numel(storeB)
    storeGXFBA1(k,:)=storeB{i};
    k=k+1;
end

%% store results for Gex model
storeGXFBA2=nan(numel(exps),4); % this to store mean correlations and percentage error
storeSTD=nan(numel(exps),4); % this to store standard deviation
folderIn='/Users/vikash/Desktop/REMI/simData/ModelsSolutions/Gex/store_'

for i=1:numel(exps)
    try
        load([folderIn exps{i} '.mat'])
        m=mean(storeRes,2);
        s=std(storeRes,0,2);
        storeGXFBA2(i,:)=m';
        storeSTD(i,:)=s';
    catch
        disp('eroor in loding')
        exps{i}
    end
end
%% store results for TGex model
storeGXFBA2_thermo=nan(numel(exps),4);
storeSTD_thermo=nan(numel(exps),4);
folderIn='/Users/vikash/Desktop/REMI/simData/ModelsSolutions/TGex/store_'

for i=1:numel(exps)
    try
        load([folderIn exps{i} '.mat'])
        m=mean(storeRes,2);
        s=std(storeRes,0,2);
        storeGXFBA2_thermo(i,:)=m';
        storeSTD_thermo(i,:)=s';
    catch
        disp('eroor in loding')
        exps{i}
    end
end

%% store results for TGexM model

storeGXFBA2_tgexm=nan(numel(exps),4);
storeSTD_tgexm=nan(numel(exps),4);
%folderIn='/Users/vikashpandey/Documents/MATLAB/EcoliRelative/GXFBA2_MET/OUTRes/store_'
folderIn='/Users/vikash/Desktop/REMI/simData/ModelsSolutions/TGexM/store_';
for i=1:numel(exps)
    try
        load([folderIn exps{i} '.mat'])
        m=mean(storeRes,2);
        s=std(storeRes,0,2);
        storeGXFBA2_tgexm(i,:)=m';
        storeSTD_tgexm(i,:)=s';
    catch
        disp('eroor in loding')
        exps{i}
    end
end

%% store results for GexM model

storeGXFBA2_gexm=nan(numel(exps),4);
storeSTD_gexm=nan(numel(exps),4);
%folderIn='/Users/vikashpandey/Documents/MATLAB/EcoliRelative/GXFBA2_MET/OUTRes/store_'
folderIn='/Users/vikash/Desktop/REMI/simData/ModelsSolutions/GexM/store_';
for i=1:numel(exps)
    try
        load([folderIn exps{i} '.mat'])
        m=mean(storeRes,2);
        s=std(storeRes,0,2);
        storeGXFBA2_gexm(i,:)=m';
        storeSTD_gexm(i,:)=s';
    catch
        disp('eroor in loding')
        exps{i}
    end
end


%% store results for TM model

storeGXFBA2_tm=nan(numel(exps),4);
storeSTD_tm=nan(numel(exps),4);
%folderIn='/Users/vikashpandey/Documents/MATLAB/EcoliRelative/GXFBA2_MET/OUTRes/store_'
folderIn='/Users/vikash/Desktop/REMI/simData/ModelsSolutions/TM/store_';
for i=1:numel(exps)
    try
        load([folderIn exps{i} '.mat'])
        m=mean(storeRes,2);
        s=std(storeRes,0,2);
        storeGXFBA2_tm(i,:)=m';
        storeSTD_tm(i,:)=s';
    catch
        disp('eroor in loding')
        exps{i}
    end
end

%% store results for M model

storeGXFBA2_m=nan(numel(exps),4);
storeSTD_m=nan(numel(exps),4);
%folderIn='/Users/vikashpandey/Documents/MATLAB/EcoliRelative/GXFBA2_MET/OUTRes/store_'
folderIn='/Users/vikash/Desktop/REMI/simData/ModelsSolutions/M/store_';
for i=1:numel(exps)
    try
        load([folderIn exps{i} '.mat'])
        m=mean(storeRes,2);
        s=std(storeRes,0,2);
        storeGXFBA2_m(i,:)=m';
        storeSTD_m(i,:)=s';
    catch
        disp('eroor in loding')
        exps{i}
    end
end


%% for ploting 
for k=1:numel(exps)-2
      y=storeGXFBA1(k,[1 2 4]);
    err=[0 0 0];
    y=[y;storeGXFBA2(k,[1 2 4])];
    err=[err;storeSTD(k,[1 2 4])];
    
   
        
    y=[y;storeGXFBA2_thermo(k,[1 2 4])];
    err=[err;storeSTD_thermo(k,[1 2 4])];
    
    if k<8
            % tgexM
            y=[y;storeGXFBA2_tgexm(k,[1 2 4])];
            err=[err;storeSTD_tgexm(k,[1 2 4])];
            % gexM
            y=[y;storeGXFBA2_gexm(k,[1 2 4])];
            err=[err;storeSTD_gexm(k,[1 2 4])];
            % tM
            y=[y;storeGXFBA2_tm(k,[1 2 4])];
            err=[err;storeSTD_tm(k,[1 2 4])];
                % M
            y=[y;storeGXFBA2_m(k,[1 2 4])];
            err=[err;storeSTD_m(k,[1 2 4])];
        end
    
     subplot(4,2,k)
     
    h = bar(y,'grouped');

    set(h,'BarWidth',1);    % The bars will now touch each other
    
   


    set(gca,'YGrid','on')

    set(gca,'GridLineStyle','-')
    
     set(gca,'XTicklabel',{'GX-FBA','Gex','TGex','TGexM','GexM','TM','M'})
%       xticks([1 2 3])
%      xticklabels({'GX-FBA','REMI-Gex','REMI-TGex'})
    title([exps{k} ' vs. Ref'], 'FontSize', 15);
    set(gca,'fontsize',15)

%     
    
    hold on
    
    numgroups = size(y, 1); 
    numbars = size(y, 2); 
    groupwidth = min(0.8, numbars/(numbars+1.5));
    
    for i = 1:numbars
        if k<8;
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        errorbar(x(2), y(2,i), err(2,i), 'k', 'linestyle', 'none'); % 2 because only we can add error bar for GXFBA2
        errorbar(x(3), y(3,i), err(3,i), 'k', 'linestyle', 'none');
        errorbar(x(4), y(4,i), err(4,i), 'k', 'linestyle', 'none');
        errorbar(x(5), y(5,i), err(5,i), 'k', 'linestyle', 'none');
        errorbar(x(6), y(6,i), err(6,i), 'k', 'linestyle', 'none');
        errorbar(x(7), y(7,i), err(7,i), 'k', 'linestyle', 'none');
        else
            x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        errorbar(x(2), y(2,i), err(2,i), 'k', 'linestyle', 'none'); % 2 because only we can add error bar for GXFBA2
        errorbar(x(3), y(3,i), err(3,i), 'k', 'linestyle', 'none');
%         errorbar(x(4), y(4,i), err(4,i), 'k', 'linestyle', 'none');
        end
    end
end

rect = [0.25, 0.5, 0.25, 0.25];

lh = legend('flux corr (cond1)','flux corr (cond2)','percent error');
set(lh, 'Position', rect, 'FontSize',15)

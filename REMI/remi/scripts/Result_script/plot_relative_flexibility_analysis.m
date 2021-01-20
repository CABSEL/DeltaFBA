% relative fluxebility analysis with and without thermodynamics
% we used flux variability analysis to compute flux range and then computed
% relative flexibility as described in main manuscript. 


%% load minMax ranges from flux variability analysis 
load('/Users/vikash/Desktop/REMI/simData/FVAMM/Combine/MM_wild_with_obj.mat')
Refmm=[refmm;refmm];
out='/Users/vikash/Desktop/REMI/simData/FVAMM/Combine/';
gexM='/Users/vikash/Desktop/REMI/simData/FVAMM/GexM/TMM';
onlyM='/Users/vikash/Desktop/REMI/simData/FVAMM/M/TMM';
experiments={'pgm vs Ref'; 'pgi vs Ref'; 'gapC vs Ref'; 'zwf vs Ref'; 'rpe vs Ref'; 'wt5 vs Ref'; 'wt7 vs Ref'; 'NOX vs Ref'; 'ATPase vs Ref'};
tGexM='/Users/vikash/Desktop/REMI/simData/FVAMM/TGexM/TMM';
onlyTM='/Users/vikash/Desktop/REMI/simData/FVAMM/TM/TMM'
cut=1e-7;
for i= 1:7
    
    try
        % tfa
        subplot(3,3,i)
        
        hold on
        load(strcat(tGexM,num2str(i),'.mat'))
        tgexM_mm=mm;
        [flexDrop1,res1]=RF(Refmm,tgexM_mm,cut);
        
        sortedVector = sort(flexDrop1);
        
        relativeCounts = (1:length(sortedVector))/length(sortedVector);
        
         plot(sortedVector,relativeCounts,'r-','LineWidth',3)
        
%           plot([mean(sortedVector) mean(sortedVector)], [0,1],'b--','LineWidth', 2)
% %         hold on
        
        
        
        %semilogx(sortedVector,relativeCounts,'b-','LineWidth',3)
        %hold on
        plot([mean(sortedVector) mean(sortedVector)], [0,1],'r--','LineWidth', 2)
        m1=mean(sortedVector);
        load(strcat(out,'outT_',num2str(i),'.mat'))
        t_mm=[[mini maxi];mm];
        
        [flexDrop1,res1]=RF(Refmm,t_mm,cut);
        
        sortedVector = sort(flexDrop1);
        
        relativeCounts = (1:length(sortedVector))/length(sortedVector);
        
         plot(sortedVector,relativeCounts,'g-','LineWidth',3)
        
        %semilogx(sortedVector,relativeCounts,'g-','LineWidth',3)
        
        plot([mean(sortedVector) mean(sortedVector)], [0,1],'g--','LineWidth', 2)
            m2=mean(sortedVector);
        load(strcat(out,num2str(i),'MM.mat'))
        f_mm=[mini maxi];
        
        [flexDrop1,res1]=RF(Refmm,f_mm,cut);
        
        sortedVector = sort(flexDrop1);
        
        relativeCounts = (1:length(sortedVector))/length(sortedVector);
        
         plot(sortedVector,relativeCounts,'b-','LineWidth',3)
%         hold on
        % TGEXM
         %semilogx(sortedVector,relativeCounts,'g-','LineWidth',3)
%         hold on
        plot([mean(sortedVector) mean(sortedVector)], [0,1],'b--','LineWidth', 2)
            m3=mean(sortedVector);
        
        % load GexM minmax 
        load(strcat(gexM,num2str(i),'.mat'))
        gexM_mm=mm;
        
        [flexDrop1,res1]=RF(Refmm,gexM_mm,cut);
        
        sortedVector = sort(flexDrop1);
        
        relativeCounts = (1:length(sortedVector))/length(sortedVector);
        
         plot(sortedVector,relativeCounts,'c-','LineWidth',3)
%         hold on
        % TGEXM
         %semilogx(sortedVector,relativeCounts,'g-','LineWidth',3)
%         hold on
        plot([mean(sortedVector) mean(sortedVector)], [0,1],'c--','LineWidth', 2)
            m4=mean(sortedVector);
          % load onlTM_minMax minmax 
        load(strcat(onlyTM,num2str(i),'.mat'))
        onlyTM_mm=mm;
        
        [flexDrop1,res1]=RF(Refmm,onlyTM_mm,cut);
        
        sortedVector = sort(flexDrop1);
        
        relativeCounts = (1:length(sortedVector))/length(sortedVector);
        
         plot(sortedVector,relativeCounts,'m-','LineWidth',3)
%         hold on
        % TGEXM
         %semilogx(sortedVector,relativeCounts,'g-','LineWidth',3)
%         hold on
        plot([mean(sortedVector) mean(sortedVector)], [0,1],'m--','LineWidth', 2)
            m5=mean(sortedVector);    
            
                        
          % load onlM_minMax minmax 
        load(strcat(onlyM,num2str(i),'.mat'))
        onlyM_mm=mm;
        
        [flexDrop1,res1]=RF(Refmm,onlyM_mm,cut);
        
        sortedVector = sort(flexDrop1);
        
        relativeCounts = (1:length(sortedVector))/length(sortedVector);
        
         plot(sortedVector,relativeCounts,'k-','LineWidth',3)
%         hold on
        % TGEXM
         %semilogx(sortedVector,relativeCounts,'g-','LineWidth',3)
%         hold on
        plot([mean(sortedVector) mean(sortedVector)], [0,1],'k--','LineWidth', 2)
            m6=mean(sortedVector);  
        
        xlabel(strcat('Relative flexibility (',experiments{i},')'))
        ylabel('Cumulative distribution')
        %legend('with thermo','with thermo avg','without thermo','without thermo avg')
        
%legend('TFA','XTFA','XFBA')
        set(gca,'fontsize',20,'linewidth',2)
        
        
    catch
        disp('error in loading')
    end
    
end
lgd=legend('TGexM','m(TGexM)','TGex','m(TGex)','Gex','m(Gex)','GexM','m(GexM)','TM','m(TM)','M','m(M)')
% 
lgd.FontSize = 10;

%% Relative flexibility of data set2

cut=1e-7;
for i=8:9
    
    try
        % tfa
        subplot(1,2,i-7)
        
        hold on
%         load(strcat(tGexM,num2str(i),'.mat'))
%         tgexM_mm=mm;
%         [flexDrop1,res1]=RF(Refmm,tgexM_mm,cut);
%         
%         sortedVector = sort(flexDrop1);
%         
%         relativeCounts = (1:length(sortedVector))/length(sortedVector);
%         
%          plot(sortedVector,relativeCounts,'r-','LineWidth',3)
%         
%           plot([mean(sortedVector) mean(sortedVector)], [0,1],'b--','LineWidth', 2)
% %         hold on
        
        
        
        %semilogx(sortedVector,relativeCounts,'b-','LineWidth',3)
        %hold on
%         plot([mean(sortedVector) mean(sortedVector)], [0,1],'r--','LineWidth', 2)
%         m1=mean(sortedVector);
        load(strcat(out,'outT_',num2str(i),'.mat'))
        t_mm=[[mini maxi];mm];
        
        [flexDrop1,res1]=RF(Refmm,t_mm,cut);
        
        sortedVector = sort(flexDrop1);
        
        relativeCounts = (1:length(sortedVector))/length(sortedVector);
        
         plot(sortedVector,relativeCounts,'g-','LineWidth',3)
        
        %semilogx(sortedVector,relativeCounts,'g-','LineWidth',3)
        
        plot([mean(sortedVector) mean(sortedVector)], [0,1],'g--','LineWidth', 2)
            m2=mean(sortedVector);
        load(strcat(out,num2str(i),'MM.mat'))
        f_mm=[mini maxi];
        
        [flexDrop1,res1]=RF(Refmm,f_mm,cut);
        
        sortedVector = sort(flexDrop1);
        
        relativeCounts = (1:length(sortedVector))/length(sortedVector);
        
         plot(sortedVector,relativeCounts,'b-','LineWidth',3)
%         hold on
        % TGEXM
         %semilogx(sortedVector,relativeCounts,'g-','LineWidth',3)
%         hold on
        plot([mean(sortedVector) mean(sortedVector)], [0,1],'b--','LineWidth', 2)
            m3=mean(sortedVector);
        
        
        
        xlabel(strcat('Relative flexibility (',experiments{i},')'))
        ylabel('Cumulative distribution')
        %legend('with thermo','with thermo avg','without thermo','without thermo avg')
        
%legend('TFA','XTFA','XFBA')
        set(gca,'fontsize',20,'linewidth',2)
        
        
    catch
        disp('error in loading')
    end
    
end
lgd=legend('TGex','mean(TGex)','Gex','mean(Gex)')





    
    





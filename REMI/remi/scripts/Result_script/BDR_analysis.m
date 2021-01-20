% With this script we identified bidirectional and uni directional reactions 


%% take BDR form TGEX and GEX model :fo this we analyzed BDR using IdentifyBDR function
load('/Users/vikash/Desktop/REMI/simData/FVAMM/Combine/BDRAnalysis.mat')
cut=1e-6;
%% load minmax analysis of GexM model
mmAll={};
load('/Users/vikash/Desktop/REMI/simData/FVAMM/GexM/TMM1.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/GexM/TMM2.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/GexM/TMM3.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/GexM/TMM4.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/GexM/TMM5.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/GexM/TMM6.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/GexM/TMM7.mat')
mmAll{end+1}=mm;
% combine BDR
 bdr_GexM={};
 for i=1:7
        t_mm=mmAll{i};
        t_bdr=IdentifyBDR(t_mm,cut);
        bdr_GexM{i}=t_bdr;
 end

%% load minmax analysis of TGexM model
mmAll={};
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TGexM/TMM1.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TGexM/TMM2.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TGexM/TMM3.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TGexM/TMM4.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TGexM/TMM5.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TGexM/TMM6.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TGexM/TMM7.mat')
mmAll{end+1}=mm;
% combine BDR
 bdr_TGexM={};
 for i=1:7
        t_mm=mmAll{i};
        t_bdr=IdentifyBDR(t_mm,cut);
        bdr_TGexM{i}=t_bdr;
 end
 
% load minmax analysis of TM model
mmAll={};
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TM/TMM1.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TM/TMM2.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TM/TMM3.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TM/TMM4.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TM/TMM5.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TM/TMM6.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/TM/TMM7.mat')
mmAll{end+1}=mm;
% combine BDR
 bdr_TM={};
 for i=1:7
        t_mm=mmAll{i};
        t_bdr=IdentifyBDR(t_mm,cut);
        bdr_TM{i}=t_bdr;
 end


% load minmax analysis of M model
mmAll={};
load('/Users/vikash/Desktop/REMI/simData/FVAMM/M/TMM1.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/M/TMM2.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/M/TMM3.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/M/TMM4.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/M/TMM5.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/M/TMM6.mat')
mmAll{end+1}=mm;
load('/Users/vikash/Desktop/REMI/simData/FVAMM/M/TMM7.mat')
mmAll{end+1}=mm;
% combine BDR
 bdr_M={};
 for i=1:7
        t_mm=mmAll{i};
        t_bdr=IdentifyBDR(t_mm,cut);
        bdr_M{i}=t_bdr;
 end

%% BDR analysis

combo=[];
for i=1:7
    Gex=[numel(bdr_fba{i}.bi) numel(find(bdr_fba{i}.bi<=2583)) numel(find(bdr_fba{i}.bi>2583))];
    TGex=[numel(bdr_tfa{i}.bi) numel(find(bdr_tfa{i}.bi<=2583)) numel(find(bdr_tfa{i}.bi>2583))];
    TGexM=[numel(bdr_TGexM{i}.bi) numel(find(bdr_TGexM{i}.bi<=2583)) numel(find(bdr_TGexM{i}.bi>2583))];
    GexM=[numel(bdr_GexM{i}.bi) numel(find(bdr_GexM{i}.bi<=2583)) numel(find(bdr_GexM{i}.bi>2583))];
    onlyTM=[numel(bdr_TM{i}.bi) numel(find(bdr_TM{i}.bi<=2583)) numel(find(bdr_TM{i}.bi>2583))];
    onlyM=[numel(bdr_M{i}.bi) numel(find(bdr_M{i}.bi<=2583)) numel(find(bdr_M{i}.bi>2583))];
    combo=[combo;Gex TGex TGexM GexM onlyTM onlyM];
end
%% with this script we will generate reltaive Relative flexibility for the each subsystems of the metabolic model


% load flux ranges from flux variability analysis

load('/Users/vikash/Desktop/REMI/simData/FVAMM/Combine/MM_wild_with_obj.mat')
Refmm=[refmm;refmm];
out='/Users/vikash/Desktop/REMI/simData/FVAMM/Combine/';
gexM='/Users/vikash/Desktop/REMI/simData/FVAMM/GexM/TMM';
onlyM='/Users/vikash/Desktop/REMI/simData/FVAMM/M/TMM';
experiments={'pgm vs Ref'; 'pgi vs Ref'; 'gapC vs Ref'; 'zwf vs Ref'; 'rpe vs Ref'; 'wt5 vs Ref'; 'wt7 vs Ref'; 'NOX vs Ref'; 'ATPase vs Ref'};
tGexM='/Users/vikash/Desktop/REMI/simData/FVAMM/TGexM/TMM';
onlyTM='/Users/vikash/Desktop/REMI/simData/FVAMM/TM/TMM'
cut=1e-7;

%% initialize variables
tmm_all=[];
fmm_all=[];
tGexM_all=[];
GexM_all=[];
M_all=[];
TM_all=[];

% combine min max from flux varaibility analysis
for i=1:9
        try
        load(strcat(out,'outT_',num2str(i),'.mat'))
        t_mm=[[mini maxi];mm];
        
        load(strcat(out,num2str(i),'MM.mat'))
        f_mm=[mini maxi];
        
        tmm_all=[tmm_all t_mm];
        fmm_all=[fmm_all f_mm];
        
        load(strcat(tGexM,num2str(i),'.mat'))
        tGexM_all=[tGexM_all mm];
        
        load(strcat(gexM,num2str(i),'.mat'))
        GexM_all=[GexM_all mm];
        
        load(strcat(onlyTM,num2str(i),'.mat'))
        TM_all=[TM_all mm];
        load(strcat(onlyM,num2str(i),'.mat'))
        M_all=[M_all mm];
        catch
            disp('error in loading')
            disp(i)
        end
 end
 
 %% Do subsystem analysis
n=7; % number of conditions
load('/Users/vikash/Desktop/REMI/data//iJO1366.mat')
usubsys=unique(iJO1366.subSystems);
mapSubSys=containers.Map();
for i=1:numel(usubsys)
    mapSubSys(usubsys{i})=find(ismember(iJO1366.subSystems,usubsys{i}));
end

pathways=mapSubSys.keys;
pathways=pathways(2:end);

avgFluxRed=nan(numel(pathways)-1,6*n);
tol=1e-6;
for k=1:numel(pathways)-1
    k
    vIdx=mapSubSys(pathways{k});
    vIdx=[vIdx;vIdx+2583];
    vv=[];
    for i=1:2:2*n
    
        
        FR_TGex=RFSubSys(Refmm(vIdx,:),tmm_all(vIdx,i:i+1),tol);
        
        FR_Gex=RFSubSys(Refmm(vIdx,:),fmm_all(vIdx,i:i+1),tol);
        
        FR_TGexM=RFSubSys(Refmm(vIdx,:),tGexM_all(vIdx,i:i+1),tol);
        FR_GexM=RFSubSys(Refmm(vIdx,:),GexM_all(vIdx,i:i+1),tol);
        FR_M=RFSubSys(Refmm(vIdx,:),M_all(vIdx,i:i+1),tol);
        FR_TM=RFSubSys(Refmm(vIdx,:),TM_all(vIdx,i:i+1),tol);
        vv=[vv [FR_TGexM.avg FR_TGex.avg FR_GexM.avg FR_Gex.avg FR_TM.avg FR_M.avg]];
       
    end
    avgFluxRed(k,:)=vv;
end
%
combo=[pathways(1:end-1)' num2cell(avgFluxRed)] % this is reduction of RF according to each subsystems
%% we want to analyze perconditions
subsys=pathways(1:end-1)';
comboAll={};

for i=1:6:6*n
    
    tmp=avgFluxRed(:,i:i+5);
    x=tmp(:,1);
    [~,I]=sort(x);
    
    comboAll=[ comboAll [subsys(I) num2cell(tmp(I,:))]];
    
end
% 
% filePh = fopen('/Users/vikashpandey/Documents/MATLAB/EcoliRelative/RepeatMET/ThermoMet/features.txt','w');
% [rows,cols]=size(comboAll);
%  for r=1:rows
%      for c=1:cols
%         fprintf(filePh,'%s\t',comboAll{r,c});
%      end
%      fprintf(filePh,'\n');
%  end

 
%% subsystem wise relative flexibility of second case


n=9; % number of conditions
load('/Users/vikash/Desktop/REMI/data//iJO1366.mat')
usubsys=unique(iJO1366.subSystems);
mapSubSys=containers.Map();
for i=1:numel(usubsys)
    mapSubSys(usubsys{i})=find(ismember(iJO1366.subSystems,usubsys{i}));
end

pathways=mapSubSys.keys;
pathways=pathways(2:end);

avgFluxRed=nan(numel(pathways)-1,2*n);
tol=1e-6;
for k=1:numel(pathways)-1
    k
    vIdx=mapSubSys(pathways{k});
    vIdx=[vIdx;vIdx+2583];
    vv=[];
    for i=1:2:2*n
    
        
        FR_TGex=RFSubSys(Refmm(vIdx,:),tmm_all(vIdx,i:i+1),tol);
        
        FR_Gex=RFSubSys(Refmm(vIdx,:),fmm_all(vIdx,i:i+1),tol);
        
       
        vv=[vv [ FR_TGex.avg FR_Gex.avg ]];
       
    end
    avgFluxRed(k,:)=vv;
end

% [crap,I]=sort(avgFluxRed(:,1));
% sortedFluxRed=avgFluxRed(I,:)
combo=[pathways(1:end-1)' num2cell(avgFluxRed)]

%% we want to analyze perconditions
subsys=pathways(1:end-1)';
comboAll={};

for i=1:2:2*n
    
    tmp=avgFluxRed(:,i:i+1);
    x=tmp(:,1);
    [~,I]=sort(x);
    
    comboAll=[ comboAll [subsys(I) num2cell(tmp(I,:))]];
    
end
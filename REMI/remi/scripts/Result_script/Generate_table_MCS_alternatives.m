%% This script is for generating Table for maximum consistency score and their alternatives
% all alternatives for both (Ishii and holm data set) are generated through scripts through "modelGenRun" folder
% for example to genrate alternatives and modles for 'Gex and TGex'we used script: create_Gex_TGeX.m'  


 load('/Users/vikash/Desktop/REMI/simData/AlternativeMCS/GeX/storeExp.mat')
 load('/Users/vikash/Desktop/REMI/simData/AlternativeMCS/TGex/storeThermoExp_modified_zwf.mat')
 % Now new anlysis is generated. 
 load('/Users/vikash/Desktop/REMI/simData/AlternativeMCS/TGexM/storeTGexM.mat')
 storeTGexM=storeTMX;
 load('/Users/vikash/Desktop/REMI/simData/AlternativeMCS/GexM/storeGexM.mat')
 %
 storeGexM=storeTMX;
 combo1=[];
% load('/Users/vikashpandey/Documents/MATLAB/EcoliRelative/AltAnalysis/EXPMetThermo/EnumTXM/storeTMX.mat')
common=storeExp{1}.onVar;
for i=1:7
    resGex=[storeExp{i}.sol.val size(storeExp{i}.zvarMat,2) numel(storeExp{i}.onVar)]; % model gene
    resTGex=[storeTExp{i}.sol.val size(storeTExp{i}.zvarMat,2) numel(storeTExp{i}.onVar)]; % model gene thermo 
    resTGexM=[storeTGexM{i}.sol.val size(storeTGexM{i}.zvarMat,2) numel(storeTGexM{i}.onVar)]; % model gene thermo metabolite
    resGexM=[storeGexM{i}.sol.val size(storeGexM{i}.zvarMat,2) numel(storeGexM{i}.onVar)]; % model gene thermo metabolite
    combo1=[combo1;[resGex resTGex resTGexM resGexM]];
    common=intersect(common,storeExp{i}.onVar);
end

% this is the table 

% With this script once can do do flux variability analysis for the REMI
% model

    
fpgm % this is the model 
cpx=changeToCPLEX(fpgm);
inds=getAllVar(fpgm,{'NF','PERTURB_NF'}); % this is to find 'NF' and 'NF' variables
path=% this is the path to save the results  
d_num % this will divide in subgroups and save results form each subgroup 
% even par for breaks after some iteration still we can save the results
[mini,maxi] = parfor_runTMMBig(cpx,inds,120,path,d_num=10)


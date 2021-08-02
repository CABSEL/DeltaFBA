load('/Users/vikashpandey/Documents/MATLAB/EcoliRelative/RepeatMET/ThermoMet/onlyMET/storeMet.mat')
addpath(genpath('/Users/vikashpandey/git/CPLEX_Studio1251'));
addpath(genpath('/Users/vikashpandey/git/FBA_Toolboxes'));
for i=2:7
    
    fpgm=storeTMX{i}.model;
    cpx=changeToCPLEX(fpgm);
    inds=getAllVar(fpgm,{'NF','PERTURB_NF'});
    path=strcat('/Users/vikashpandey/Documents/MATLAB/EcoliRelative/RepeatMET/ThermoMet/onlyMET/MM/TMM',num2str(i),'index1.mat')
    [mini,maxi] = parfor_runTMMBig(cpx,inds,120,path,10);
     mm=[mini maxi];
    save(strcat('/Users/vikashpandey/Documents/MATLAB/EcoliRelative/RepeatMET/ThermoMet/onlyMET/MM/TMM',num2str(i),'.mat'),'mm')
end



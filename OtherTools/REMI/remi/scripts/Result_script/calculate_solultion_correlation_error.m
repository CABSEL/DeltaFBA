% With this script we would are calculating solutions by minimizing sum of fluxes
% we then compared solutions with experimental data which includes some
% fluxomics for core carbon metabolism

%% we are taking an example to show TGexM model 
% this analysis can be done with all diffrent kind of models 


exps=textread('/Users/vikash/Desktop/REMI/data/expType.txt','%s'); % this is mutant and wild types
exps=strrep(exps,'.mat',''); % remove the .mat
% set Path for cplex solver
% addpath(genpath('/Users/vikashpandey/git/CPLEX_Studio1251')); 
% addpath(genpath('/Users/vikashpandey/git/FBA_Toolboxes'));
path='/Users/vikash/Desktop/REMI/simData/AlternativeMCS/TGexM/Enum'; % this path contains number of alternative solutions for
                                                                      % each consistency score  
load('/Users/vikash/Desktop/REMI/data/iJO1366.mat') % load genome scale model

result_save= % SET PATH TO STORE ALL RESULTS: CORRLEATION, SOLUTIONS, PERCENTAGE EROOR/store_'
% We already gave saved results in 'ModelsSolutions' folder
rxns=iJO1366.rxns;
ishii=1:7;
holm=8:9;
for i=1:numel(exps)
    try
        load([path exps{i} '.mat'])  % loded enumerated solutions
        clear model % delete model 
        % load model for each condition and form specific type: We are
        % taking TGexM
        load(['/Users/vikash/Desktop/REMI/simData/ModelsSolutions/TGexM/' exps{i} '.mat']);
        fg1=find(ishii==i);
        fg2=find(holm==i);
        if numel(fg1)>0
            folderIn='/Users/vikash/Desktop/REMI/simData/FluxData/' % data is taken from E-Flux2 and SPOT: Validated Methods for Inferring Intracellular Metabolic Flux Distributions from Transcriptomic Data
            FLUX_LIST1 = [folderIn 'fluxmatched.txt'];
            MEASURED_FLUX1 = [folderIn 'measured_flux_rf'  '.txt'];
            FLUX_LIST2 = [folderIn 'fluxmatched.txt'];
            MEASURED_FLUX2 = [folderIn 'measured_flux_' exps{i}  '.txt'];
            % model is the specific model for each condition 
            
            [storeRes1,storeRes2,storeRes3,storeRes4,storeXX]=getCorrErorrPar(model,rxns,z_matrix,FLUX_LIST1,MEASURED_FLUX1,FLUX_LIST2,MEASURED_FLUX2);
            storeRes=[storeRes1;storeRes2;storeRes3;storeRes4];
%             save([result_save exps{i} '.mat'], 'storeRes','storeXX'); %
%             This results are already saved in ModelsSolutions folder
        elseif numel(fg2)>0
            folderIn='/Users/vikash/Desktop/REMI/simData/FluxData/'
            FLUX_LIST1 = [folderIn 'holm_match_flux.txt'];
            MEASURED_FLUX1 = [folderIn 'measured_flux_Ref'  '.txt'];
            FLUX_LIST2 = [folderIn 'holm_match_flux.txt'];
            MEASURED_FLUX2 = [folderIn 'measured_flux_' exps{i}  '.txt'];
            
           [storeRes1,storeRes2,storeRes3,storeRes4,storeXX]=getCorrErorrPar(model,rxns,z_matrix,FLUX_LIST1,MEASURED_FLUX1,FLUX_LIST2,MEASURED_FLUX2);
            storeRes=[storeRes1;storeRes2;storeRes3;storeRes4];
%             save([result_save exps{i} '.mat'], 'storeRes','storeXX');
            % This results are already saved in ModelsSolutions folder
        else
            disp('index is not from holm or ishii')
        end
    catch
        disp('error in loding')
        exps{i}
    end
end
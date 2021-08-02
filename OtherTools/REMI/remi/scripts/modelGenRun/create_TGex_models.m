%% this script for creating models for ishii data set
%% addPaths

% addpath(genpath('/Users/vikashpandey/git/CPLEX_Studio1251'));
% addpath(genpath('/Users/vikashpandey/git/FBA_Toolboxes'));


storeModels={};
for mut_exp=2:8
    %% initialize varibales
    wt=1;
    mut=mut_exp;
    tarns_cut=0.1; % transcript cutoff
    up=.05;
    down=.05;
    
    %% This data from Elux2 paper
    load('/Users/vikash/Desktop/REMI/data/iJO1366.mat')
    load('/Users/vikash/Desktop/REMI/data/test_expr.mat')
    genes=expr{1};
    genes_data=expr{2};
    
    %% check FBA model
    % fmodel=iJO1366;
    % fmodel.rev=ones(numel(fmodel.rev),1);
    % fmodel=addUseVariablesDH(fmodel);
    % fmodel=addNetFluxVariablesNEW(fmodel);
    % sol=solveTFBAmodel(fmodel);
    % initialize FBAmodel with maximum flux
    load('/Users/vikash/Desktop/REMI/data/iJO1366_thermo.mat')
    fmodel=thermoM_true;
    m_fmodel=fmodel;
    sol=solveTFBAmodel(fmodel);
    m_fmodel.var_lb(find(fmodel.f))=sol.val*0.8;
    
    % find reaction and reaction ratios between two conditions
    geneWT=nan(numel(fmodel.genes),1);
    geneMUT=nan(numel(fmodel.genes),1);
    for i=1:numel(fmodel.genes)
        [~,ind]=ismember(fmodel.genes{i},genes);
        if ind>0
            geneWT(i)=genes_data(ind,wt);
            geneMUT(i)=genes_data(ind,mut);
        end
    end
    
    % calculate ratio (if transcrip label is too samll then we dont cal culate ratio)
    ind1=find(geneWT(geneWT>tarns_cut));
    ind2=find(geneMUT(geneMUT>tarns_cut));
    indC=intersect(ind1,ind2); % common index
    ratio=geneMUT(indC)./geneWT(indC);
    [tmpModel]=createEnumModel(m_fmodel,m_fmodel,m_fmodel.genes(indC),ratio,up,down); % genes and ratio should be same length
    storeModels{end+1}=tmpModel;
end


%% models for holm data set


% addpath(genpath('/Users/vikashpandey/git/CPLEX_Studio1251'));
% addpath(genpath('/Users/vikashpandey/git/FBA_Toolboxes'));
% % read transcriptomics data 

EXPR_DATA_FNAME = '/Users/vikash/Desktop/REMI/simData/Expression_data/holm_transcriptomics.csv'; % txt file name where gene expresison data is stored

fid = fopen(EXPR_DATA_FNAME);
C = textscan(fid, '%s%f%f%f', 'Headerlines', 1, 'Delimiter', '\t'); % first column : gene names used in the model, the other columns : measured gene expression data

fclose(fid);
gene_names = C{1};
expr_data = cell2mat(C(2:end));

expr = cell(1, 2);
expr{1} = gene_names;
expr{2} = expr_data;


storeB={};
for mut_exp=2:3
    %% initialize varibales
    wt=1;
    mut=mut_exp;
    tarns_cut=0.1; % transcript cutoff
    up=.05;
    down=.05;
    
    %% This data from Elux2 paper
 load('/Users/vikash/Desktop/REMI/data/iJO1366.mat')
    
    genes=expr{1};
    genes_data=expr{2};
    
    %% check FBA model
    % fmodel=iJO1366;
    % fmodel.rev=ones(numel(fmodel.rev),1);
    % fmodel=addUseVariablesDH(fmodel);
    % fmodel=addNetFluxVariablesNEW(fmodel);
    % sol=solveTFBAmodel(fmodel);
    % initialize FBAmodel with maximum flux
    load('/Users/vikash/Desktop/REMI/data/iJO1366_thermo.mat')
    fmodel=thermoM_true;
    m_fmodel=fmodel;
    sol=solveTFBAmodel(fmodel);
    m_fmodel.var_lb(find(fmodel.f))=sol.val*0.8;
    
    % find reaction and reaction ratios between two conditions
    geneWT=nan(numel(fmodel.genes),1);
    geneMUT=nan(numel(fmodel.genes),1);
    for i=1:numel(fmodel.genes)
        [~,ind]=ismember(fmodel.genes{i},genes);
        if ind>0
            geneWT(i)=genes_data(ind,wt);
            geneMUT(i)=genes_data(ind,mut);
        end
    end
    
    % calculate ratio (if transcrip label is too samll then we dont cal culate ratio)
    ind1=find(geneWT(geneWT>tarns_cut));
    ind2=find(geneMUT(geneMUT>tarns_cut));
    indC=intersect(ind1,ind2); % common index
    ratio=geneMUT(indC)./geneWT(indC);
    [tmpModel]=createEnumModel(m_fmodel,m_fmodel,m_fmodel.genes(indC),ratio,up,down); % genes and ratio should be same length
    storeModels{end+1}=tmpModel;
end
% save('iShii_and_holm_Model.mat','storeModels')
% 
% exp={ 'pgm', 'pgi', 'gapC', 'zwf', 'rpe', 'wt5', 'wt7', 'NOX',	'ATPase'};
% fid=fopen('expType.txt','w');
% for i=1:numel(exp)
%     fprintf(fid,'%s\n',[exp{i} '.mat']);
%     model=storeModels{i};
%     sol=solveTFBAmodel(model);
%     model.var_lb(end)=sol.val;
%     save([exp{i} '.mat'],'model')
% end
% exps={ 'pgm', 'pgi', 'gapC', 'zwf', 'rpe', 'wt5', 'wt7', 'NOX',	'ATPase'};
% for i=1:numel(exps)
%     path_save=['/Users/vikashpandey/Documents/MATLAB/EcoliRelative/ThermoGXFBA2/EnumModels/enumT' exps{i} '.mat'];
%     clear model
%     load(['/Users/vikashpandey/Documents/MATLAB/EcoliRelative/ThermoGXFBA2/EnumModels/' exps{i} '.mat'])
%     findAltCombi(1000,model,model.objIndex1B,path_save)
% end
% 
% exps={  'NOX',	'ATPase'};
% for i=1:numel(exps)
%     path_save=['/Users/vikashpandey/Documents/MATLAB/EcoliRelative/ThermoGXFBA2/Res/enumT' exps{i} '.mat'];
%     clear model
%     load(['/Users/vikashpandey/Documents/MATLAB/EcoliRelative/ThermoGXFBA2/EnumModels/' exps{i} '.mat'])
%     findAltCombi(1000,model,model.objIndex1B,path_save)
% end
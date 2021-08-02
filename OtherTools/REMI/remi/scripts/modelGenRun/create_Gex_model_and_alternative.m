%% this script for creating Gex model form models for ishii and holm data data set
% we also used this script for latwrrnative anaysis
% addpath(genpath('/Users/vikashpandey/git/CPLEX_Studio1251'));
% addpath(genpath('/Users/vikashpandey/git/FBA_Toolboxes'));    add Path
% for cplex solver

% build a model for pgm vs wild type

storeModels={};
%% we are creating model form ishii data set initialize varibales
for mut_exp=2:8
 wt=1;
 mut=mut_exp;% this is for pgm 
 tarns_cut=0.1; % transcript cutoff
 up=.05; % this is the cutoff for slecting up regulated 5 %
 down=.05; % this is for the downregulated  5 %
    
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
    load('/Users/vikash/Desktop/REMI/data/fmodel.mat')
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
    tmpModel=createEnumModel(m_fmodel,m_fmodel,m_fmodel.genes(indC),ratio,up,down); % genes and ratio should be same length
    storeModels{end+1}=tmpModel;
end

%% we are creating models from holm data set
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
    fmodel=iJO1366;
    fmodel.rev=ones(numel(fmodel.rev),1);
    fmodel=addUseVariablesDH(fmodel);
    fmodel=addNetFluxVariablesNEW(fmodel);
    % sol=solveTFBAmodel(fmodel);
    % initialize FBAmodel with maximum flux
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

%% This is for alternative analysis
exp={ 'pgm', 'pgi', 'gapC', 'zwf', 'rpe', 'wt5', 'wt7'};
fid=fopen('expType.txt','w');
for i=1:numel(exp)
    fprintf(fid,'%s\n',[exp{i} '.mat']);
    model=storeModels{i};
    index=[model.objIndex1B;model.metB];
%     f=zeros(numel(model.f),1);
%     f(index)=1;
%     model.f=f;
    sol=solveTFBAmodel(model);
    model.var_lb(end)=sol.val; % set maximum consistency
%     save(['TGexM/' exp{i} '.mat'],'model')
      %  need to save this models we saved in Folder:
      %  simData/ModelsSolutions/GeX
end

%%
exp={ 'pgm', 'pgi', 'gapC', 'zwf', 'rpe', 'wt5', 'wt7'};
for i=1:numel(exp)
    path_save=['simData/AlternativeMCS/Gex/Enum' exp{i} '.mat']; % once need to adjust Path 
    clear model
    load(['simData/ModelsSolutions/Gex/' exp{i} '.mat']); % this is the path of the model which one we created above
    index=[model.objIndex1B;model.metB];
%     f=zeros(numel(model.f),1);
%     f(index)=1;
%     model.f=f;
   
    findAltCombi(1000,model,index,path_save,300) % this script find alternative solutions at maximum consistency score
end
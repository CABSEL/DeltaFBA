% addpath(genpath('/Users/vikashpandey/git/CPLEX_Studio1251'));
% addpath(genpath('/Users/vikashpandey/git/FBA_Toolboxes'));
% read transcriptomics data 

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
    % fmodel.rev=ones(numel(fmodel.rev),1);
    % fmodel=addUseVariablesDH(fmodel);
    % fmodel=addNetFluxVariablesNEW(fmodel);
    % sol=solveTFBAmodel(fmodel);
    % initialize FBAmodel with maximum flux
%     load('fmodel.mat')
    m_fmodel=fmodel;
%     sol=solveFBACplex(fmodel);
%     m_fmodel.lb(find(fmodel.f))=sol.obj_val;
    
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
    changeCobraSolver('gurobi')
    [gxFBAsolution,wflux] = gxFBA(m_fmodel,m_fmodel,m_fmodel.genes(indC),ratio,1)
    
    
    measured_condition=1;
    eflux2_flux = wflux.x;
    folderIn='/Users/vikash/Desktop/REMI/simData/FluxData/'
    FLUX_LIST = [folderIn 'holm_match_flux.txt'];
    condition = {'Ref',	'NOX',	'ATPase'};
    MEASURED_FLUX = [folderIn 'measured_flux_' condition{measured_condition} '.txt'];
    
    [correl_m1,m1,p1]=getMeasuredPrediction(eflux2_flux,FLUX_LIST,MEASURED_FLUX,measured_condition);
    
    measured_condition=mut;
    eflux2_flux = gxFBAsolution.x;
    folderIn='/Users/vikash/Desktop/REMI/simData/FluxData/'
    FLUX_LIST = [folderIn 'holm_match_flux.txt'];
    ccondition = {'Ref',	'NOX',	'ATPase'};
    MEASURED_FLUX = [folderIn 'measured_flux_' condition{measured_condition} '.txt'];
    
    
    [correl_m2,m2,p2]=getMeasuredPrediction(eflux2_flux,FLUX_LIST,MEASURED_FLUX,measured_condition);
    xx=[m1 m2 p1 p2];
    
    % rel1=m2./m1;
    % rel2=p2./p1;
    %
    % % see the binary on
    %
    % active=tmpModel.varNames(find(sol.x(tmpModel.objIndex1B)<0.1));
    % active=strrep(active,'YBplus_F','');
    % active=strrep(active,'YBplus_R','');
    % active=unique(active);
    % % find net flux
    % [~,a_nfi]=ismember(strcat('NF_',active),tmpModel.varNames);
    % [~,a_pnfi]=ismember(strcat('PERTURB_NF_',active),tmpModel.varNames);
    
    % calculate relative difference
    %% skip nan zero rows
    xx(any(isnan(xx),2),:) = [];
    %%
    flux_cut=1e-05;
    xx(xx<flux_cut & xx>-flux_cut)=0;
    exp=xx(:,[1 2]);
    sim=xx(:,[3 4]);
    d_exp=(exp(:,2)-exp(:,1))./(abs(exp(:,2))+abs(exp(:,1)));
    % find zero rows
    d_exp(all(exp==0,2))=0;
    
    d_sim=(sim(:,2)-sim(:,1))./(abs(sim(:,2))+abs(sim(:,1)));
    % find zero rows
    d_sim(all(sim==0,2))=0;
    % find correlation
    correl_m = dot(d_sim,d_exp)/(norm(d_sim)*norm(d_exp));
    percentErorr=mean(abs(d_sim-d_exp));
    storeB{end+1}=[correl_m1 correl_m2 correl_m percentErorr];
    
end

% we run two times 1) to generate and saved the results and
% ModelsSolutuons/GXFBA
% save('storeBHolm_gxfbaMin.mat','storeB')


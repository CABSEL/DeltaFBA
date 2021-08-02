function fluxes = call_RELATCH(model, gene_names, gene_exp, external_rxns, external_rates)
% Calls the implementation of RELATCH provided in [Kim and Reed, Genome Biol, 2012].
%
% Note: Only the RELATCH_Reference.m method is necessary.
%
% INPUTS
%       model - cobra model
%       gene_names - gene ids
%       gene_exp - gene expression
%       external_rxns - measured exchange reactions
%       external_rates - measured rates
%
% OUTPUTS
%       fluxes - flux distribution
%
% Author: Daniel Machado, 2013 

    
    % Full deletions must be removed to avoid division by zero
    relatch_gene_data.genes = gene_names(gene_exp > 0);
    gene_exp = gene_exp(gene_exp > 0);
    
    %RELATCH is sensitive to the scaling of the gene expression
    %so it's better to scale it using the mean
    relatch_gene_data.val = log(gene_exp / mean(gene_exp));
    
    external_flux.rxns = external_rxns;
    external_flux.val = external_rates;
    external_flux.err = zeros(size(external_flux.val));
    
    
    solutionRef = RELATCH_Reference(model, relatch_gene_data, external_flux);
    
    if solutionRef.stat == 1
        fluxes = solutionRef.w;
    else
        fluxes = [];
    end
end


function solutionRef = RELATCH_Reference(model,Gene_Expression,External_Flux,MFA_Flux,solver)
%RELATCH_Reference estimates the flux distribution and enzyme contributions
%in a reference state
%
%   solutionRef = RELATCH_Reference(model,Gene_Expression,External_Flux,MFA_Flux,solver)
%
%INPUTS
% model                     Reference metabolic model (CobraToolBox model with GPR)
% Gene_Expression
%   Gene_Expression.genes   Name of genes
%   Gene_Expression.val     Log2 normalized signal intensity
% External_Flux
%   External_Flux.rxns      Name of external fluxes
%   External_Flux.val       Measured flux values
%   External_Flux.err       Standard deviations
% MFA_Flux
%   MFA_Flux.rxns           Name of MFA fluxes
%   MFA_Flux.val            MFA flux values
%   MFA_Flux.err            MFA confidence intervals
% solver                    QP solver ('CobraQP' (default) or 'cplex_direct')
%
%OUTPUTS
% solutionRef
%   solutionRef.w           Reference flux distribution
%   solutionRef.W           Reference enzyme contributions
%   solutionRef.stat        Reference solution status
%
%Notes:
% RELATCH_Reference requires CobraToolBox 2.0 and a LP/QP solver
%
% Joonhoon Kim and Jennifer L. Reed 9/13/2012

if (nargin < 2)
    Gene_Expression = [];
end
if (nargin < 3)
    External_Flux = [];
end
if (nargin < 4)
    MFA_Flux = [];
end
if (nargin < 5)
    solver = 'CobraQP';
end

fprintf('Solving a FBA problem..\n');

% Solve FBA problem
solutionFBA = optimizeCbModel(model);

if solutionFBA.stat < 1
    error('Reference FBA problem is infeasible or unbounded');
end

% Generate GPR mapping
fprintf('Generating GPR mapping..\n');

[nMets,nRxns] = size(model.S);
% Find reactions with GPR
GPR_ID=find(~cellfun(@isempty,model.rules));
% Exclude Porin transport reactions and spontaneous reactions
Porin_ID=find(~cellfun(@isempty,strfind(model.subSystems,'Porin')));
Spont_ID=find(strcmp(model.grRules,'s0001'));
Excl_Rxn_ID=union(Porin_ID,Spont_ID);
GPR_ID=setdiff(GPR_ID,Excl_Rxn_ID);
% Find reversible reactions with GPR
[GPRrev_ID, GPRrev_ID2]=intersect(GPR_ID,find(model.rev));
nRxnsGPR=size(GPR_ID,1);
nRxnsGPRrev=size(GPRrev_ID,1);

% Generate GPR matrix
GPR=regexp(model.rules(GPR_ID),'\|','split');
GPR=cellfun(@(c) regexp(c(:),'x\((\d+)\)','tokens'), GPR, 'UniformOutput', 0);
Enz=cellfun(@(c) str2double([c{:}]), vertcat(GPR{:}), 'UniformOutput',0);
nEnzs=size(Enz,1);
nGenes=size(model.genes,1);

numEnz=zeros(nRxnsGPR,1);
for i=1:nRxnsGPR; numEnz(i)=size(GPR{i},1); end;
Rxn2Enz=sparse(nRxnsGPR,nEnzs);
for i=1:nRxnsGPR; Rxn2Enz(i,sum(numEnz(1:i-1))+1:sum(numEnz(1:i)))=ones(1,numEnz(i)); end;
numSub=zeros(nEnzs,1);
for i=1:nEnzs; numSub(i)=numel(Enz{i}); end;
Enz2Gene=sparse(nEnzs,nGenes);
for i=1:nEnzs; Enz2Gene(i,Enz{i})=1; end;

% Setting up the problem
% Variables in the following problem are
% x = [w;W]
% where w = flux distribution in the reference state
%       W = Enzyme contribution in the reference state
fprintf('Setting up the problem..\n');

F1 = sparse(nRxns,nRxns);
c1 = zeros(nRxns,1);

if ~isempty(External_Flux)
    [dummy, Ext_ID]=ismember(External_Flux.rxns,model.rxns);
    if find(Ext_ID==0)
        error('Exchange reactions not found in the model');
    else
        External_Flux.err(External_Flux.err==0)=1e-3;
        F1(Ext_ID,Ext_ID)=diag(External_Flux.err.^(-2));
        c1(Ext_ID)=-External_Flux.val./(External_Flux.err.^2);
        Nouse_Exch_ID=find(findExcRxns(model) & (solutionFBA.x==0));
        model.lb(Nouse_Exch_ID)=0;
        model.ub(Nouse_Exch_ID)=0;
        model.lb(Ext_ID) = External_Flux.val - External_Flux.err;
        model.ub(Ext_ID) = External_Flux.val + External_Flux.err;
    end
end

if ~isempty(MFA_Flux)
    [dummy, MFA_ID]=ismember(MFA_Flux.rxns,model.rxns);
    if find(MFA_ID==0)
        error('MFA reactions not found in the model');
    else
        MFA_Flux.err(MFA_Flux.err==0)=1e-3;
        F1(MFA_ID,MFA_ID)=diag(MFA_Flux.err.^(-2));
        c1(MFA_ID)=-MFA_Flux.val./(MFA_Flux.err.^2);
    end
end

if ~isempty(Gene_Expression)
    [Gene_ID, Exp_ID]=ismember(Gene_Expression.genes,model.genes);
    Exp_ID(Exp_ID==0)=[];
    Gene_Expression.val(~Gene_ID)=[];
    Exp_val=zeros(nGenes,1);
    Exp_val(Exp_ID)=exp(Gene_Expression.val);
    Exp_val(Exp_val==0)=exp(mean(Gene_Expression.val));
    F2 = diag(1./(Enz2Gene*Exp_val));
else
    F2 = sparse(1:nEnzs,1:nEnzs,1e-3);
end

if isempty(find(F1))
    model.lb(model.c==1)=solutionFBA.f;
    model.ub(model.c==1)=solutionFBA.f;
end

model.lb(model.lb==-1000)=-inf;
model.ub(model.ub==1000)=inf;

QP.F = [ F1 sparse(nRxns,nEnzs) ;
         sparse(nEnzs,nRxns) F2 ];
QP.c = [ c1 ; zeros(nEnzs,1) ];
QP.A = [ model.S sparse(nMets,nEnzs) ;
         sparse(1:nRxnsGPR,GPR_ID,-1,nRxnsGPR,nRxns) Rxn2Enz ;
         sparse(1:nRxnsGPRrev,GPRrev_ID,1,nRxnsGPRrev,nRxns) Rxn2Enz(GPRrev_ID2,:) ];
QP.b = zeros(nMets+nRxnsGPR+nRxnsGPRrev,1);
QP.lb = [ model.lb ; zeros(nEnzs,1) ];
QP.ub = [ model.ub ; inf*ones(nEnzs,1) ];
QP.csense(1:nMets) = 'E';
QP.csense(nMets+1:nMets+nRxnsGPR+nRxnsGPRrev) = 'G';
QP.osense = 1;

% Solve the problem
fprintf('Solving reference state estimation problem..\n');

if strcmp(solver,'CobraQP')
    % Solve the reference state estimation problem using CobraQP solver
    QPsol = solveCobraQP(QP, 'printLevel',1);
    if (QPsol.stat > 0)
        solutionRef.w = QPsol.full(1:nRxns);
        solutionRef.w(abs(solutionRef.w) < 1e-6) = 0;
        solutionRef.W = QPsol.full(nRxns+1:nRxns+nEnzs);
        solutionRef.W(solutionRef.W < 1e-6) = 0;
        solutionRef.stat = QPsol.stat;
        fprintf('Done\n\n');
    else
        solutionRef.stat = QPsol.stat;
        warning('Referense state estimation problem is infeasible or unbounded');
    end
elseif strcmp(solver,'cplex_direct')
    % Solve the reference state estimation problem using ILOG CPLEX
    [x,obj,exitflag,output] = cplexqp(QP.F,QP.c,-QP.A(nMets+1:end,:),-QP.b(nMets+1:end),QP.A(1:nMets,:),QP.b(1:nMets),QP.lb,QP.ub);
    if (exitflag > 0)
        solutionRef.w = x(1:nRxns);
        solutionRef.w(abs(solutionRef.w) < 1e-6) = 0;
        solutionRef.W = x(nRxns+1:nRxns+nEnzs);
        solutionRef.W(solutionRef.W < 1e-6) = 0;
        solutionRef.stat = exitflag;
        solutionRef.output = output;
        fprintf('Done\n\n');
    else
        solutionRef.stat = exitflag;
        solutionRef.output = output;
        warning('Referense state estimation problem is infeasible or unbounded');
    end
else error('Unknown QP solver');
end

end


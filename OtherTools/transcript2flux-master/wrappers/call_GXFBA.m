function fluxes = call_GXFBA(model, gene_names, gene_exp, gene_exp_ref, wt_ref)
% Call GX-FBA using the code provided in [Navid and Almaas, BMC Sys Bio, 2012].
%
% Note: The code was modified to re-use the reference in multiple calls.
%
% INPUTS
%       model - cobra model
%       gene_names - gene ids
%       gene_exp - gene expression
%       gene_exp_ref - gene expression of reference condition
%       wt_ref - reference flux distribution and variability
%
% OUTPUTS
%       fluxes - flux distribution
%
% Author: Daniel Machado, 2013 

    levels = log2(gene_exp ./ gene_exp_ref);
    
    sol = gxFBA(wt_ref, model, gene_names, levels);
    fluxes = sol.x;
end

function [gxFBAsolution,genestat,Z,gxFVAmin,gxFVAmax,Zgx,Zgmin,Zgmax] = gxFBA(wt_ref,model2,geneList,exprVect,exprType,threshold,iter,verbFlag)
%
%USAGE
% [gxFBAsolution,genestat,Z,gxFVAmin,gxFVAmax,Zgx,Zgmin,Zgmax] = 
%              gxFBA(model1,model2,geneList,exprVect,exprType,threshold,iter,verbFlag)
%
%DESCRIPTION
% gxFBA performs Gene-Expression FBA analysis following Navid & Almaas, BMC Syst. Biol. (2012).
%   This is gxFBA version 1.0
%
%INPUT
% model1            COBRA model structure, with 1 reaction assigned as 
%                    objective to be maximized, intended for the reference state
%
% model2           COBRA model structure, intended for the perturbed or stressed state.
%
% geneList         A (cell array) list of genes for which there is 
%                    expression data. 
%                    Example: {'g1','g2','g3'}
%
% exprVect         A vector of gene expression data for the genes in the
%                    variable geneList.
%                    Example: [-3.22,2.86,1.27] for log2(ratio), thus
%                    up-regulation corresponds to positive entries.
%
%OPTIONAL INPUTS
% exprType         gene expression data type:
%                      1  - if gene expression ratio, only positive values
%                      2  - if log2(gene expression ratio)  (default)
%
% threshold        minumum fold change to be considered up or down regulated 
%                    (default=0.5, corresponding to 50% increase)
%
% iter             Number of samples of degenerate FBA reference flux 
%                    states to estimate effect of reference state flux 
%                    degeneracy on GX-FBA result (default = 0)
%
% verbFlag         Level of feedback (default true)
% 
%OUTPUT
% gxFBAsolution    structured array object from interior point solve with fields
%    f             Objective value
%    x             Primal
%    y             Dual
%    w             Reduced costs
%    stat          Solver status in standardized form
%                    1   Optimal solution
%                    2   Unbounded solution
%                    0   Infeasible
%                   -1   No solution reported
%
%OPTIONAL OUTPUT
% genestat         structured array object with gene use status fields
%    match         list of genes with match in the model
%    nomatch       list of genes with no match in the model
%
% Z                Vector of GX-FBA objective function coefficients
%
% gxFVAmin         Min Flux Variability values in GX-FBA optimal state
%
% gxFVAmax         Max Flux Variability values in GX-FBA optimal state
% 
% Zgx              Variation in GX-FBA optimal value from sampling flux
%                     degeneracy in the FBA reference state.
% Zgmin            Variation in Growth_min rate in GX-FBA optimal state
%                     due to degeneracy in the FBA ref. state
% Zgmax            Variation in Growth_max rate in GX-FBA optimal state
%                     due to degeneracy in the FBA ref. state
%
%
% Please cite: Navid & Almaas, BMC Syst Biol (2012)
%
% Eivind Almaas  25/09/2012
%
%DEPENDENCIES
%
% optimizeCbModel()
% fluxVariability()
% changeRxnBounds()
%
%%


% FIXES: Daniel M. - don't recompute FVA
% FIXES: Daniel M. - use cobra's optimizeCbModel


fprintf('\nBeginning GX-FBA calculation.\n\n');
if (nargin < 4)
    fprintf('  --> Not enough input arguments. \n  --> Please check and try again.\n\n');
    gxFBAsolution=0;Z=0;gxFVAmin=0;gxFVAmax=0;Zgx=0;Zgmin=0;Zgmax=0;genestat=0;
    return
end
if (nargin < 5)
    exprType = 2;
end
if (nargin < 6)
    threshold = 0.5;
end
if (nargin < 7)
    iter = 0;
end
if (nargin < 8)
    verbFlag = true;
end
if (nargin > 9)
    fprintf('  --> Too many input arguments.\n  --> Please check and try again.\n\n');
    gxFBAsolution=0;Z=0;gxFVAmin=0;gxFVAmax=0;Zgx=0;Zgmin=0;Zgmax=0;genestat=0;
    return
end
threshold = log2(1+threshold);
if (length(geneList) ~= length(exprVect))
    fprintf('  --> The number of entries in gene list and expression list does not agree.\n  --> Please check and try again.\n\n');
    gxFBAsolution=0;Z=0;gxFVAmin=0;gxFVAmax=0;Zgx=0;Zgmin=0;Zgmax=0;genestat=0;
    return
end

%Local variables
eps        = 0.01;    %the smallest allowable gene expression ratio
eps2       = 0.00001; %the smallest flux value used in the GX-FBA modifications
eps3       = 0.1;     %the absolute smallest allowable diff between FVA max and min
                      %  in the GX-FBA robustness calculation
modelgxfba = model2;  %the stressed or perturbed state model
t0         = cputime;

%Find biomass reaction
BM_react = find(model2.c > 0.5);

%Calculate reference state
if (verbFlag)
    fprintf('Calculating FBA reference state ... \n');
end

%Must choose InteriorPoint method, not simplex
%   This function is a fix until cobra toolbox allows more parameters to 
%   be passed directly to the solvers.
%wt  = optimizeCbModel_int(model1,'max');

% FIXED: Daniel, use reference
wt  = wt_ref.wt_sol;
acc = find(abs(wt.x) > eps2); 
lac = length(acc);
if lac > 0
    wtc = sum( wt.x(acc) )/lac;
else
    fprintf('  --> FBA reference state carries no flux. Check problem.\n\n');
    gxFBAsolution=0;Z=0;gxFVAmin=0;gxFVAmax=0;Zgx=0;Zgmin=0;Zgmax=0;genestat=0;
    return
end 
    
%Conduct FVA analysis of FBA reference state
if (verbFlag)
    fprintf('Calculating FVA for FBA optimal state ...');
end

%FVA analysis at 0% of FBA optimum to find characteristic flux scale for
%  given reference environemnt.
% FIX: Daniel - fixed to receive as arguments
wt_minf = wt_ref.wt_minf;
wt_maxf = wt_ref.wt_maxf;
wt_avg             = (wt_maxf+wt_minf)/2;

t1 = cputime;
if (verbFlag)
    fprintf('%10.3f s\n',t1-t0);
end

%Identify T, the set of irreversible reactions for which a consistent
%  gene expression signal is available
if (verbFlag)
    fprintf('Building GX-FBA problem description ...\n');
end

%1. extract genes with match in model for stressed state
[inModel,geneInd] = ismember(geneList,model2.genes);
geneID1           = geneInd(inModel);
expVal1           = exprVect(inModel);
genestat.match    = geneList(inModel);
genestat.nomatch  = geneList( not(inModel) );

%2. extract genes with expression above threshold, convert expression
% log2() data to expression ratios
if (exprType == 1)
    expVal1 = log2(expVal1);
end
pos     = find( (expVal1 < -threshold) | (expVal1 > threshold) );
geneID2 = geneID1(pos);
expVal2 = pow2(expVal1(pos));
 
%3. determine reaction set T: find affected irreversible reactions and 
%   check for consistent signal
rxExMap  = containers.Map('KeyType','int32','ValueType','double');
remRxMap = containers.Map('KeyType','int32','ValueType','int32');     
remGX    = containers.Map('KeyType','int32','ValueType','int32'); 
wtZ      = containers.Map('KeyType','int32','ValueType','int32'); 

for i=1:length(geneID2)
    d = find(model2.rxnGeneMat(:,geneID2(i)) == 1);

    for j=1:length(d)  %each gene may lead to multiple reactions

        %only continue with irreversible, or effectively irreversible, reactions
        if ( model2.rev(d(j)) == 0 || abs(wt_maxf(d(j))-wt_minf(d(j))) < abs(wt_maxf(d(j)))+eps2 )
       
            if isKey(remRxMap,d(j))  %must avoid conflicted reactions
               continue
            end
            
            if isKey(rxExMap,d(j))    %we must check existing entries for 
                                      %  gene expression consistency               
                if (rxExMap(d(j)) < 1) && (expVal2(i) < 1)
                    if rxExMap(d(j)) > expVal2(i)
                        rxExMap(d(j)) = expVal2(i);
                    end
                elseif (rxExMap(d(j)) > 1) && (expVal2(i) > 1)
                    if rxExMap(d(j)) < expVal2(i)
                        rxExMap(d(j)) = expVal2(i);
                    end
                else
                    %conflicting gene expression data, must be excluded
                    remRxMap(d(j)) = 0;
                    remove(rxExMap,d(j));
                end                
            else
                %reaction entry doesn't previously exists
                rxExMap(d(j)) = expVal2(i);
            end
        end
    end
end

%4. Set new upper and lower bounds for reactions in set T and construct
%   new objective function
%   First, deal with irreversible reactions with zero WT flux.
cand = find( (wt.x <= eps2) & (model2.rev == 0) );
rem  = find( wt_maxf(cand) < eps2);
for i=1:length(rem)
    remGX(cand(rem(i))) = 0;
end
keep = find( wt_maxf(cand) > eps2);
for i=1:length(keep)
    wtZ(cand(keep(i))) = 0;
end

%remove irreversible reactions that are unconstrained in the given
%  environemnt, so that they will not enter the gxFBA goal function.
rem  = find( (wt_maxf > 500) & (model2.rev == 0) );  
for i=1:length(rem)
    remGX(rem(i)) = 0;
end

modelgxfba.c = zeros(size(modelgxfba.c));
keyset = keys(rxExMap);
for i=keyset
    rx = i{1};
    if isKey(remGX,rx)  %Skip these reactions
        continue
    end
    if isKey(wtZ,rx)  %zero wt-flux reactions
        if rxExMap(rx) > 1
            modelgxfba.ub(rx) = rxExMap(rx)*wtc;
            %not adding lower bound for down regulation since wt-zero-flux reaction
        end
    else  %normal, irreversible reaction with GX data
        if rxExMap(rx) < 1
            modelgxfba.lb(rx) = rxExMap(rx)*wt.x(rx);
        else        
            modelgxfba.ub(rx) = rxExMap(rx)*wt.x(rx);
        end
    end
    modelgxfba.c(rx) = log2( rxExMap(rx) ) / wt_avg(rx);
end
Z = modelgxfba.c;

%Run GX-FBA simulation
if (verbFlag)
    fprintf('Solving GX-FBA problem ...\n');
end

%gxFBAsolution = optimizeCbModel_int(modelgxfba,'max');
gxFBAsolution = optimizeCbModel(modelgxfba,'max');

end

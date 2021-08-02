function [gxFBAsolution,wt] = gxFBA(model1,model2,geneList,exprVect,exprType,threshold,iter,verbFlag)
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
BM_react = find(model1.c > 0.5);

%Calculate reference state
if (verbFlag)
    fprintf('Calculating FBA reference state ... \n');
end

%Must choose InteriorPoint method, not simplex
%   This function is a fix until cobra toolbox allows more parameters to 
%   be passed directly to the solvers.
% wt  = optimizeCbModel_int(model1,'max');
wt=solveFBACplex(model1);
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
% [wt_minf, wt_maxf] = fluxVariability(model1,0,'max');  
% save('wild_flux_vari.mat','wt_minf', 'wt_maxf')
load('/Users/vikashpandey/Documents/MATLAB/EcoliRelative/GX-FBA2/GXFBA/wild_flux_vari.mat')

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

%
rId=cell2mat(rxExMap.keys);

cand = find( (wt.x <= eps2) & (model2.rev == 0) );
rem  = find( wt_maxf(cand) < eps2);
for i=1:length(rem)
    remGX(cand(rem(i))) = 0;
end
keep = find( wt_maxf(cand) > eps2);
for i=1:length(keep)
    wtZ(cand(keep(i))) = 0;
end
keep1 = find( wt.x(rId)==0);
for i=1:length(keep1)
    wtZ(rId(keep1(i))) = 0;
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
        modelgxfba.c(rx) = log2( rxExMap(rx) ) / wt_avg(rx);
    end
    
end
Z = modelgxfba.c;

%Run GX-FBA simulation
if (verbFlag)
    fprintf('Solving GX-FBA problem ...\n');
end
% gxFBAsolution = optimizeCbModel_int(modelgxfba,'max');
gxFBAsolution=solveFBACplex(modelgxfba)
saveSol=gxFBAsolution;
%% add additional constraint
modelgxfba.mets=[modelgxfba.mets;'newMet'];
modelgxfba.metNames=[modelgxfba.metNames;'newMet'];
modelgxfba.b=[modelgxfba.b;gxFBAsolution.objval];
cons=zeros(1,numel(modelgxfba.c));
cons(find(modelgxfba.c))=-1;
modelgxfba.S=[modelgxfba.S;cons];
modelgxfba.c=ones(numel(modelgxfba.c),1)*-1;


gxFBAsolution=solveFBACplex(modelgxfba);
if isequal(gxFBAsolution.status,3)
    gxFBAsolution=saveSol;
end
%%
% % if (verbFlag)
% %     fprintf('Calculating FVA in GX-FBA optimal state ...');
% % end
% % [gxFVAmin,gxFVAmax] = fluxVariability(modelgxfba,100,'max');
% % t0 = t1;
% % t1 = cputime;
% % if (verbFlag)
% %     fprintf('%10.3f s\n',t1-t0);
% % end
% 
% %Calculate sensitivity / robustness of GX-FBA result to flux degeneracy in 
% %  reference FBA state    
% if iter > 0
%     if (verbFlag)
%         fprintf('\nStarting GX-FBA robustness analysis ...\n')
%     end
%     
%     %identify reactions with degeneracy on FBA optimal surface
%     if (verbFlag)
%         fprintf('Calculating FVA in FBA optimal reference state to asses variability ... ');
%     end
%     [wt_minf, wt_maxf] = fluxVariability(model1,100,'max');  
%     if (verbFlag)
%         fprintf('%10.3f s\n',t1-t0);
%     end
%     wt_diff = find(abs(wt_minf - wt_maxf) > eps3);
%     ent     = length(wt_diff);
%     if (verbFlag)
%         fprintf('   %i reactions are degenerate on FBA optimal surface.\n',ent);
%     end
%     if (ent > 0)
%         Zgx   = zeros(iter,1);
%         Zgmin = zeros(iter,1);
%         Zgmax = zeros(iter,1);
%             
%         h = waitbar(0,'GX-FBA robustness analysis in progress ...');
%         div = floor(0.1*iter);
%         for i=1:iter
%         
%             if (mod(i,div) == 0)
%                 waitbar(i/iter,h);
%             end       
%             gxfba_it   = model2;
%             gxfba_it.c = modelgxfba.c;
% 
%             %Pick random reaction among ones degenerate on FBA optimal surface
%             r        = wt_diff( randi(ent) );
%             flux_r   = wt_minf(r) + ( wt_maxf(r)- wt_minf(r) )*rand;
%             model_it = changeRxnBounds(model1,model1.rxns(r),flux_r,'b');
%             alt_f    = optimizeCbModel_int(model_it,'max');
% 
%             %Identify new set of zero wt-flux reactions
%             wtZ_it   = containers.Map('KeyType','int32','ValueType','int32');
%             remGX_it = containers.Map('KeyType','int32','ValueType','int32'); 
%             cand     = find( (alt_f.x <= eps2) & (model2.rev == 0) );
%             rem      = find( wt_maxf(cand) < eps2);
%             for j=1:length(rem)
%                 remGX_it(cand(rem(j))) = 0;
%             end
%             keep = find( wt_maxf(cand) > eps2);
%             for j=1:length(keep)
%                 wtZ_it(cand(keep(j))) = 0;
%             end
%             rem  = find( (wt_maxf > 500) & (model2.rev == 0) );
%             for j=1:length(rem)
%                 remGX_it(rem(j)) = 0;
%             end
%             
%             %calculate new upper and lower bounds for reactions
%             for j=keyset 
%                 rx = j{1};
%                 if isKey(remGX_it,rx)  %Skip these reactions
%                     continue
%                 end
%                 if isKey(wtZ_it,rx)  %zero wt-flux reactions
%                     if rxExMap(rx) > 1
%                         gxfba_it.ub(rx) = rxExMap(rx)*wtc;
%                         %not adding lower bound for down regulation since wt-zero-flux reaction
%                     end
%                 else  %normal, irreversible reaction with GX data
%                     if rxExMap(rx) < 1
%                         gxfba_it.lb(rx) = rxExMap(rx)*alt_f.x(rx);
%                     else
%                         gxfba_it.ub(rx) = rxExMap(rx)*alt_f.x(rx);
%                     end
%                 end
%             end
% 
%             %Run GX-FBA simulation
%             itsol               = optimizeCbModel_int(gxfba_it,'max');
%             Zgx(i)              = itsol.f;
%             [Zgmin(i),Zgmax(i)] = fluxVariability(gxfba_it,100,'max',model1.rxns(BM_react(1)));
%         end
%         close(h);   
%         t0 = t1;
%         t1 = cputime;
%         if (verbFlag)
%             fprintf('\n ... completed in %10.3f s\n',t1-t0);
%         end
%     else
%         fprintf('No degeneracy in FBA optimal state.\n');
%         Zgx = 0; Zgmin = 0; Zgmax = 0;
%     end
% else
%    Zgx = 0; Zgmin = 0; Zgmax = 0;
% end
fprintf('\nGX-FBA calculation completed.\n\n');
end

% End of GX-FBA. The rest is minor code fix on existing cobra functions to 
% use interior point solver with either glpk or gurobi    


%% optimizeCbModel_int
function [FBAsolution] = optimizeCbModel_int(model,osenseStr,minNorm,allowLoops)

%optimizeCbModel Solve a flux balance analysis problem
%
%MODIFICATION OF ORIGINAL CODE SO THAT INTERIOR POINT SOLVER CAN BE USED
% E.A.
%
% Solves LP problems of the form: max/min c'*v
%                                 subject to S*v = b
%                                            lb <= v <= ub
% FBAsolution = optimizeCbModel(model,osenseStr,minNormFlag)
%
%INPUT
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%
%OPTIONAL INPUTS
% osenseStr      Maximize ('max')/minimize ('min') (opt, default = 'max')
%
% minNorm        {(0), 'one', > 0 , n x 1 vector}, where [m,n]=size(S);
%                0      Default, normal LP
%                'one'  Minimise the Taxicab Norm using LP.
%                                 min |v|
%                                   s.t. S*v = b
%                                        c'v = f
%                                        lb <= v <= ub
%                -----
%                The remaining options work only with a valid QP solver:
%                -----
%                > 0    Minimises the Euclidean Norm of internal fluxes.
%                       Typically 1e-6 works well.
%                                 min ||v||
%                                   s.t. S*v = b
%                                        c'v = f
%                                        lb <= v <= ub
%               n x 1   Forms the diagonal of positive definiate
%                       matrix F in the quadratic program
%                               min 0.5*v'*F*v
%                               st. S*v = b
%                                   c'*v = f
%                                   lb <= v <= ub
%
% allowLoops    {0,(1)} If true, then instead of a conventional FBA,
%               the solver will run an MILP version which does not allow
%               loops in the final solution.  Default is true.
%               Runs much slower when set to false.
%               See addLoopLawConstraints.m to for more info.
%
%OUTPUT
% FBAsolution
%   f         Objective value
%   x         Primal
%   y         Dual
%   w         Reduced costs
%   s         Slacks
%   stat      Solver status in standardized form
%              1   Optimal solution
%              2   Unbounded solution
%              0   Infeasible
%             -1  No solution reported (timelimit, numerical problem etc)
%
%
% Markus Herrgard       9/16/03
% Ronan Fleming         4/25/09  Option to minimises the Euclidean Norm of internal
%                                fluxes using 'cplex_direct' solver
% Ronan Fleming         7/27/09  Return an error if any imputs are NaN
% Ronan Fleming         10/24/09 Fixed 'E' for all equality constraints
% Jan Schellenberger             MILP option to remove flux around loops
% Ronan Fleming         12/07/09 Reworked minNorm parameter option to allow
%                                the full range of approaches for getting
%                                rid of net flux around loops.
% Jan Schellenberger    2/3/09   fixed bug with .f being set incorrectly
%                                when minNorm was set.
% Nathan Lewis          12/2/10  Modified code to allow for inequality
%                                constraints.
% Ronan Fleming         12/03/10 Minor changes to the internal handling of global parameters.
%% Process arguments and set up problem

allowLoops = 1;
minNorm    = 0;

if exist('osenseStr', 'var')
    if isempty(osenseStr)
        osenseStr = 'max';
    end
else
    osenseStr = 'max';
end
% Figure out objective sense
if strcmpi(osenseStr,'max')
    LPproblem.osense = -1;
else
    LPproblem.osense = +1;
end

if exist('minNorm', 'var')
    if isempty(minNorm)
        %use global solver parameter for minNorm
        minNorm = getCobraSolverParams('LP','minNorm');
    end
else
    %use global solver parameter for minNorm
    minNorm = getCobraSolverParams('LP','minNorm');
end
if exist('allowLoops', 'var')
    if isempty(allowLoops)
        allowLoops = true;
    end
else
    allowLoops = true;
end

%use global solver parameter for printLevel
[printLevel,primalOnlyFlag] = getCobraSolverParams('LP',{'printLevel','primalOnly'});

[nMets,nRxns] = size(model.S);

% add csense
%Doing this makes csense a double array.  Totally smart design move.
%LPproblem.csense = [];
if ~isfield(model,'csense')
    % If csense is not declared in the model, assume that all
    % constraints are equalities.
    LPproblem.csense(1:nMets,1) = 'E';
else % if csense is in the model, move it to the lp problem structure
    if length(model.csense)~=nMets,
        warning('Length of csense is invalid! Defaulting to equality constraints.')
        LPproblem.csense(1:nMets,1) = 'E';
    else
        model.csense = columnVector(model.csense);
        LPproblem.csense = model.csense;
    end
end

% Fill in the RHS vector if not provided
if (~isfield(model,'b'))
    LPproblem.b = zeros(size(model.S,1),1);
else
    LPproblem.b = model.b;
end

% Rest of the LP problem
LPproblem.A = model.S;
LPproblem.c = model.c;
LPproblem.lb = model.lb;
LPproblem.ub = model.ub;

%Double check that all inputs are valid:
if ~(verifyCobraProblem(LPproblem, [], [], false) == 1)
    warning('invalid problem');
    return;
end

%%
t1 = clock;
% Solve initial LP
if allowLoops
    solution = solveCobraLP_int(LPproblem);
end

if (solution.stat ~= 1) % check if initial solution was successful.
    if printLevel>0
        warning('Optimal solution was not found');
    end
    FBAsolution.f = 0;
    FBAsolution.x = [];
    FBAsolution.stat = solution.stat;
    FBAsolution.origStat = solution.origStat;
    FBAsolution.solver = solution.solver;
    FBAsolution.time = etime(clock, t1);
    return;
end

objective = solution.obj; % save for later use.

% Store results
if (solution.stat == 1)
    %solution found.
    FBAsolution.x = solution.full(1:nRxns);
    
    %this line IS necessary.
    FBAsolution.f = model.c'*solution.full(1:nRxns); %objective from original optimization problem.
    if abs(FBAsolution.f - objective) > .01
        if strcmp(minNorm,'one')
            display('optimizeCbModel.m warning:  objective appears to have changed while minimizing taxicab norm');
        else
            error('optimizeCbModel.m: minimizing Euclidean norm did not work')
        end
    end
    
    %if (~primalOnlyFlag && allowLoops && any(~minNorm)) % LP rcost/dual only correct if not doing minNorm
    % LP rcost/dual are still meaninful if doing, one simply has to be aware that there is a
    % perturbation to them the magnitude of which depends on norm(minNorm) - Ronan   
    if (~primalOnlyFlag && allowLoops)
        FBAsolution.y = solution.dual;
        FBAsolution.w = solution.rcost;
    end
else
    %some sort of error occured.
    if printLevel>0
        warning('Optimal solution was not found');
    end
    FBAsolution.f = 0;
    FBAsolution.x = [];
end

FBAsolution.stat = solution.stat;
FBAsolution.origStat = solution.origStat;
FBAsolution.solver = solution.solver;
FBAsolution.time = etime(clock, t1);
end
%%


function solution = solveCobraLP_int(LPproblem, varargin)
%solveCobraLP Solve constraint-based LP problems
%
% solution = solveCobraLP(LPproblem, parameters)
%
%INPUT
% LPproblem Structure containing the following fields describing the LP
% problem to be solved
%  A      LHS matrix
%  b      RHS vector
%  c      Objective coeff vector
%  lb     Lower bound vector
%  ub     Upper bound vector
%  osense Objective sense (-1 max, +1 min)
%  csense Constraint senses, a string containting the constraint sense for
%         each row in A ('E', equality, 'G' greater than, 'L' less than).
%
%OPTIONAL INPUTS
% Optional parameters can be entered using parameters structure or as
% parameter followed by parameter value: i.e. ,'printLevel',3)
%
% parameters    Structure containing optional parameters as fields.
%               Setting parameters = 'default' uses default setting set in
%               getCobraSolverParameters.
% printLevel    Printing level
%               = 0    Silent (Default)
%               = 1    Warnings and Errors
%               = 2    Summary information 
%               = 3    More detailed information
%               > 10   Pause statements, and maximal printing (debug mode)
% saveInput     Saves LPproblem to filename specified in field. 
%               i.e. parameters.saveInput = 'LPproblem.mat';
% minNorm       {(0), scalar , n x 1 vector}, where [m,n]=size(S); 
%               If not zero then, minimise the Euclidean length 
%               of the solution to the LP problem. minNorm ~1e-6 should be
%               high enough for regularisation yet maintain the same value for 
%               the linear part of the objective. However, this should be
%               checked on a case by case basis, by optimization with and
%               without regularisation.
% primalOnly    {(0),1} 1=only return the primal vector (lindo solvers)
%               
% optional parameters can also be set through the
% solver can be set through changeCobraSolver('LP', value);
% changeCobraSolverParames('LP', 'parameter', value) function.  This
% includes the minNorm and the printLevel flags
%
%OUTPUT
% solution Structure containing the following fields describing a LP
% solution
%  full     Full LP solution vector
%  obj      Objective value
%  rcost    Reduced costs
%  dual     Dual solution
%  solver   Solver used to solve LP problem
%
%  stat     Solver status in standardized form
%            1   Optimal solution
%            2   Unbounded solution
%            0   Infeasible
%           -1   No solution reported (timelimit, numerical problem etc)
%
%  origStat Original status returned by the specific solver
%  time     Solve time in seconds
%
%
% Markus Herrgard    08/29/06
% Ronan Fleming      11/12/08 'cplex_direct' allows for more refined control
%                             of cplex than tomlab tomrun
% Ronan Fleming      04/25/09 Option to minimise the Euclidean Norm of internal
%                             fluxes using either 'cplex_direct' solver or 'pdco'
% Jan Schellenberger 09/28/09 Changed header to be much simpler.  All parameters
%                             now accessed through 
%                             changeCobraSolverParams(LP, parameter,value)
% Richard Que        11/30/09 Changed handling of optional parameters to use
%                             getCobraSolverParams().
% Ronan Fleming      12/07/09 Commenting of input/output
% Ronan Fleming      21/01/10 Not having second input, means use the parameters as specified in the
%                             global paramerer variable, rather than 'default' parameters
% Steinn Gudmundsson 03/03/10 Added support for the Gurobi solver

%% Process arguments etc

global CBTLPSOLVER
if (~isempty(CBTLPSOLVER))
    solver = CBTLPSOLVER;
else
    error('No solver found.  call changeCobraSolver(solverName)');
end
optParamNames = {'minNorm','printLevel','primalOnly','saveInput', ...
    'feasTol','optTol','EleNames','EqtNames','VarNames','EleNameFun', ...
    'EqtNameFun','VarNameFun','PbName','MPSfilename'};
parameters = '';
if nargin ~=1
    if mod(length(varargin),2)==0
        for i=1:2:length(varargin)-1
            if ismember(varargin{i},optParamNames)
                parameters.(varargin{i}) = varargin{i+1};
            else
                error([varargin{i} ' is not a valid optional parameter']);
            end
        end
    elseif strcmp(varargin{1},'default')
        parameters = 'default';
    elseif isstruct(varargin{1})
        parameters = varargin{1};
    else
        display('Warning: Invalid number of parameters/values')
        solution=[];
        return;
    end
end
[minNorm, printLevel, primalOnlyFlag, saveInput, feasTol, optTol] = ...
    getCobraSolverParams('LP',optParamNames(1:6),parameters);


%Save Input if selected
if ~isempty(saveInput)
    fileName = parameters.saveInput;
    if ~find(regexp(fileName,'.mat'))
        fileName = [fileName '.mat'];
    end
    display(['Saving LPproblem in ' fileName]);
    save(fileName,'LPproblem')
end


[A,b,c,lb,ub,csense,osense] = deal(LPproblem.A,LPproblem.b,LPproblem.c,LPproblem.lb,LPproblem.ub,LPproblem.csense,LPproblem.osense);

% if any(any(~isfinite(A)))
%     error('Cannot perform LP on a stoichiometric matrix with NaN of Inf coefficents.')
% end

% Defaults in case the solver does not return anything
f = [];
x = [];
y = [];
w = [];
origStat = -99;
stat = -99;

t_start = clock;
switch solver
    %% GLPK
    case 'glpk'
        params.msglev = printLevel; % level of verbosity
        params.tolbnd = feasTol; %tolerance
        params.toldj = optTol; %tolerance
        params.lpsolver = 2; %Interior point solver
        if (isempty(csense))
            clear csense
            csense(1:length(b),1) = 'S';
        else
            csense(csense == 'L') = 'U';
            csense(csense == 'G') = 'L';
            csense(csense == 'E') = 'S';
            csense = columnVector(csense);
        end
        %glpk needs b to be full, not sparse -Ronan
        b=full(b);
        [x,f,y,w,stat,origStat] = solveGlpk_int(c,A,b,lb,ub,csense,osense,params);

    case 'gurobi'
        %% gurobi
        % Free academic licenses for the Gurobi solver can be obtained from
        % http://www.gurobi.com/html/academic.html
        %
        % The code below uses Gurobi Mex to interface with Gurobi. It can be downloaded from
        % http://www.convexoptimization.com/wikimization/index.php/Gurobi_Mex:_A_MATLAB_interface_for_Gurobi

        clear opts            % Use the default parameter settings
        if printLevel == 0
           % Version v1.10 of Gurobi Mex has a minor bug. For complete silence
           % Remove Line 736 of gurobi_mex.c: mexPrintf("\n"); 
           opts.Display = 0;
           opts.DisplayInterval = 0;
        else
           opts.Display = 1;
        end

        opts.FeasibilityTol = feasTol;
        opts.OptimalityTol = optTol;
        opts.Method        = 2;  %Interior point solver
        
        if (isempty(csense))
            clear csense
            csense(1:length(b),1) = '=';
        else
            csense(csense == 'L') = '<';
            csense(csense == 'G') = '>';
            csense(csense == 'E') = '=';
            csense = csense(:);
        end
	%gurobi_mex doesn't cast logicals to doubles automatically
	c = double(c);
        [x,f,origStat,output,y] = gurobi_mex(c,osense,sparse(A),b, ...
                                             csense,lb,ub,[],opts);
        if origStat==2
           stat = 1; % Optimal solutuion found
        elseif origStat==3
           stat = 0; % Infeasible
        elseif origStat==5
           stat = 2; % Unbounded
        elseif origStat==4
           stat = 0; % Gurobi reports infeasible *or* unbounded
        else
           stat = -1; % Solution not optimal or solver problem
        end
        
    case 'pdco'
        %-----------------------------------------------------------------------
        % pdco.m: Primal-Dual Barrier Method for Convex Objectives (16 Dec 2008)
        %-----------------------------------------------------------------------
        % AUTHOR:
        %    Michael Saunders, Systems Optimization Laboratory (SOL),
        %    Stanford University, Stanford, California, USA.
        %Interfaced with Cobra toolbox by Ronan Fleming, 27 June 2009
        [nMet,nRxn]=size(LPproblem.A);
        x0 = ones(nRxn,1);
        y0 = ones(nMet,1);
        z0 = ones(nRxn,1);
 
        %setting d1 to zero is dangerous numerically, but is necessary to avoid 
        %minimising the Euclidean norm of the optimal flux. A more
        %numerically stable way is to use pdco via solveCobraQP, which has
        %a more reasonable d1 and should be more numerically robust. -Ronan
        d1=0; 
        d2=1e-6;
        options = pdcoSet;
        options.FeaTol    = 1e-12;
        options.OptTol    = 1e-12;
        %pdco is a general purpose convex optization solver, not only a
        %linear optimization solver. As such, much control over the optimal
        %solution and the method for solution is available. However, this
        %also means you may have to tune the various parameters here,
        %especially xsize and zsize (see pdco.m) to get the real optimal
        %objective value
        xsize = 1000;
        zsize = 10000;
        
        options.Method=2; %QR
        options.MaxIter=100;
        options.Print=printLevel;
        [x,y,w,inform,PDitns,CGitns,time] = ...
            pdco(osense*c*10000,A,b,lb,ub,d1,d2,options,x0,y0,z0,xsize,zsize);
        f= c'*x;
        % inform = 0 if a solution is found;
%        = 1 if too many iterations were required;
%        = 2 if the linesearch failed too often;
%        = 3 if the step lengths became too small;
%        = 4 if Cholesky said ADDA was not positive definite.
        if (inform == 0)
            stat = 1;
        elseif (inform == 1 || inform == 2 || inform == 3)
            stat = 0;
        else
            stat = -1;
        end
        origStat=inform;
    case 'mps'
        %% BuildMPS
        % This calls buildMPS and generates a MPS format description of the
        % problem as the result
        % Build MPS Author: Bruno Luong
        % Interfaced with CobraToolbox by Richard Que (12/18/09)
        display('Solver set to MPS. This function will output an MPS matrix string for the LP problem');
        %Get optional parameters
        [EleNames,EqtNames,VarNames,EleNameFun,EqtNameFun,VarNameFun,PbName,MPSfilename] = ...
            getCobraSolverParams('LP',{'EleNames','EqtNames','VarNames','EleNameFun','EqtNameFun','VarNameFun','PbName','MPSfilename'},parameters);
        %split A matrix for L and E csense
        Ale = A(csense=='L',:);
        ble = b(csense=='L');
        Aeq = A(csense=='E',:);
        beq = b(csense=='E');
        
        %%%%Adapted from BuildMPS%%%%%
        [neq nvar]=size(Aeq);
        nle=size(Ale,1);
        if isempty(EleNames)
            EleNames=arrayfun(EleNameFun,(1:nle),'UniformOutput', false);
        end
        if isempty(EqtNames)
            EqtNames=arrayfun(EqtNameFun,(1:neq),'UniformOutput', false);
        end
        if isempty(VarNames)
            VarNames=arrayfun(VarNameFun,(1:nvar),'UniformOutput', false);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [solution] = BuildMPS(Ale, ble, Aeq, beq, c, lb, ub, PbName,'MPSfilename',MPSfilename,'EleNames',EleNames,'EqtNames',EqtNames,'VarNames',VarNames);
        
        
    otherwise
        error(['For GX-FBA, please use either glpk or gurobi instead of: ' solver]);
        
end
if ~strcmp(solver,'cplex_direct') && ~strcmp(solver,'mps')
    %% Assign solution
    t = etime(clock, t_start);
    if ~exist('basis','var'), basis=[]; end
    [solution.full,solution.obj,solution.rcost,solution.dual,solution.solver,solution.stat,solution.origStat,solution.time,solution.basis] = ...
        deal(x,f,w,y,solver,stat,origStat,t,basis);
end
end

%% solveGlpk Solve actual LP problem using glpk and return relevant results
function [x,f,y,w,stat,origStat] = solveGlpk_int(c,A,b,lb,ub,csense,osense,params)

% Old way of calling glpk
%[x,f,stat,extra] = glpkmex(osense,c,A,b,csense,lb,ub,[],params);
[x,f,origStat,extra] = glpk(c,A,b,lb,ub,csense,[],osense,params);
y = extra.lambda;
w = extra.redcosts;
% Note that status handling may change (see glplpx.h)
if (origStat == 180 || origStat == 5)
    stat = 1; % Optimal solution found
elseif (origStat == 182 || origStat == 183 || origStat == 3 || origStat == 110)
    stat = 0; % Infeasible
elseif (origStat == 184 || origStat == 6)
    stat = 2; % Unbounded
else
    stat = -1; % Solution not optimal or solver problem
    fprintf('\nYour version of the mex file: glpk does not support interior point.\n ');
    fprintf('Suggestion: 1. recompile the glpk mex yourself\n');
    fprintf('            2. use gurobi: its standard mex works with interior point.\n\n\n')
end
end






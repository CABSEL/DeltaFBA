function fluxes = call_iMAT(model, gene_names, gene_exp, lower_threshold, upper_threshold, eps_param)
% Implements iMAT as defined in [Shlomi et al, Nat Biotech, 2008].
% Adapted from the implementation provided in the cobra toolbox.
%
% INPUTS
%       model - cobra model
%       gene_names - genes ids
%       gene_exp - genes expression
%       lower_threshold - lower expression threshold
%       upper_threshold - upper expression threshold
%       eps_param - flux activation threshold
%
% OUTPUTS
%       fluxes - flux distribution
%
% Author: Daniel Machado, 2013 

    discrete_levels = zeros(size(gene_exp));
    discrete_levels(gene_exp > upper_threshold) = 1;
    discrete_levels(gene_exp < lower_threshold) = -1;
    reaction_levels = gene_to_reaction_levels(model, gene_names, discrete_levels, @min, @max);
    RHindex = find(reaction_levels > 0);
    RLindex = find(reaction_levels < 0);
    
    fluxes = shlomi(model, RHindex, RLindex, eps_param);

end

function v_sol = shlomi(model, RHindex, RLindex, eps_param)
% Implementation from the cobra toolbox (createTissueSpecificModel.m)

    S = model.S;
    lb = model.lb;
    ub = model.ub;

    % Creating A matrix
    A = sparse(size(S,1)+2*length(RHindex)+2*length(RLindex),size(S,2)+2*length(RHindex)+length(RLindex));
    [m,n,s] = find(S);
    for i = 1:length(m)
        A(m(i),n(i)) = s(i); %#ok<SPRIX>
    end

    for i = 1:length(RHindex)
        A(i+size(S,1),RHindex(i)) = 1; %#ok<SPRIX>
        A(i+size(S,1),i+size(S,2)) = lb(RHindex(i)) - eps_param; %#ok<SPRIX>
        A(i+size(S,1)+length(RHindex),RHindex(i)) = 1; %#ok<SPRIX>
        A(i+size(S,1)+length(RHindex),i+size(S,2)+length(RHindex)+length(RLindex)) = ub(RHindex(i)) + eps_param; %#ok<SPRIX>
    end

    for i = 1:length(RLindex)
        A(i+size(S,1)+2*length(RHindex),RLindex(i)) = 1; %#ok<SPRIX>
        A(i+size(S,1)+2*length(RHindex),i+size(S,2)+length(RHindex)) = lb(RLindex(i)); %#ok<SPRIX>
        A(i+size(S,1)+2*length(RHindex)+length(RLindex),RLindex(i)) = 1; %#ok<SPRIX>
        A(i+size(S,1)+2*length(RHindex)+length(RLindex),i+size(S,2)+length(RHindex)) = ub(RLindex(i)); %#ok<SPRIX>
    end

    % Creating csense
    csense1(1:size(S,1)) = 'E';
    csense2(1:length(RHindex)) = 'G';
    csense3(1:length(RHindex)) = 'L';
    csense4(1:length(RLindex)) = 'G';
    csense5(1:length(RLindex)) = 'L';
    csense = [csense1 csense2 csense3 csense4 csense5];

    % Creating lb and ub
    lb_y = zeros(2*length(RHindex)+length(RLindex),1);
    ub_y = ones(2*length(RHindex)+length(RLindex),1);
    lb = [lb;lb_y];
    ub = [ub;ub_y];

    % Creating c
    c_v = zeros(size(S,2),1);
    c_y = ones(2*length(RHindex)+length(RLindex),1);
    c = [c_v;c_y];

    % Creating b
    b_s = zeros(size(S,1),1);
    lb_rh = lb(RHindex);
    ub_rh = ub(RHindex);
    lb_rl = lb(RLindex);
    ub_rl = ub(RLindex);
    b = [b_s;lb_rh;ub_rh;lb_rl;ub_rl];

    % Creating vartype
    vartype1(1:size(S,2),1) = 'C';
    vartype2(1:2*length(RHindex)+length(RLindex),1) = 'B';
    vartype = [vartype1;vartype2];

    MILPproblem.A = A;
    MILPproblem.b = b;
    MILPproblem.c = c;
    MILPproblem.lb = lb;
    MILPproblem.ub = ub;
    MILPproblem.csense = csense;
    MILPproblem.vartype = vartype;
    MILPproblem.osense = -1;
    MILPproblem.x0 = [];

    params.timeLimit = 100;
    solution = solveCobraMILP(MILPproblem, params);

    x = solution.cont;
    for i = 1:length(x)
        if abs(x(i)) < 1e-6
            x(i,1) = 0;
        end
    end

    v_sol = x;

end
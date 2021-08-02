function fluxes = call_Lee12(model, gene_names, gene_exp, gene_exp_sd, gene_scale_rxn, flux_scale_rxn, flux)
% Calls the implementation of Lee-12 provided in [Lee et al, BMC Sys Bio, 2012].
% 
% Note: fixes to the original code are tagged with FIX
%
% INPUTS
%       model - cobra model
%       gene_names - genes ids
%       gene_exp - gene expression
%       gene_exp_sd - gene expression std
%       gene_scale_rxn - reaction to scale gene expression
%       flux_scale_rxn - reaction to scale flux distribution
%       flux - flux of scaling reaction
%
% OUTPUTS
%       fluxes - flux distribution
%
% Author: Daniel Machado, 2013 

    scale_gene_idx = find(strcmp(gene_scale_rxn, model.rxns));
    scale_rxn_idx = find(strcmp(flux_scale_rxn, model.rxns));

    
    % fill in zero entries
    unmeasured = setdiff(model.genes, gene_names);
    gene_names = [gene_names; unmeasured];
    gene_exp = [gene_exp; zeros(length(unmeasured),1)];
    
    if isempty(gene_exp_sd)
        gene_exp_sd = zeros(size(gene_exp));
        disp('warning: no gene expression std given')
    else
        gene_exp_sd = [gene_exp_sd; zeros(length(unmeasured),1)];
    end
    
    [rxn_exp,rxn_exp_sd] = geneToReaction(model,gene_names,gene_exp,gene_exp_sd);
    if rxn_exp(scale_gene_idx) == 0
        error('Error: gene expression for scaling enzyme is zero.');
    end
    rxn_exp = rxn_exp/rxn_exp(scale_gene_idx);
    if any(rxn_exp_sd>0)
        rxn_exp_sd(rxn_exp_sd == 0) = min(rxn_exp_sd(rxn_exp_sd>0))/2;
        rxn_exp_sd = rxn_exp_sd/rxn_exp(scale_gene_idx);
    end
    
    % FIX: scale all bounds accordingly, not just glc uptake
    model.lb = model.lb / abs(flux);
    model.ub = model.ub / abs(flux);

    % Gene expression constraint FBA
    fluxes = dataToFlux(model, rxn_exp, rxn_exp_sd);

    % scale back by uptk rate
    if fluxes(scale_rxn_idx) == 0
        error('Error: flux for scaling reaction is zero.');
    end
    fluxes = fluxes * abs(flux / fluxes(scale_rxn_idx));

end


function [r,r_sd] = geneToReaction(m,g,t,t_sd)

    % kieran: 16 sep 11

    r       = zeros(size(m.rxns));
    r_sd    = zeros(size(m.rxns));

    for k = 1:length(g)
        g{k} = strrep(g{k},'-','_');
    end

    for k = 1:length(m.rxns)
        ga = m.grRules{k};
        ga = strrep(ga,'-','_');
        w = regexp(ga,'\<\w*\>','match'); 
        w = setdiff(w,{'and','or'});
        for kk = 1:length(w)
            j = find(strcmp(w{kk},g));
            n = t(j);
            n_sd = t_sd(j);        
            ga = regexprep(ga,['\<',w{kk},'\>'],[num2str(n),'pm',num2str(n_sd)]);
        end
        [n,n_sd] = addGeneData(ga);
        r(k) = n;
        r_sd(k) = n_sd;
    end

end

function str = AandB(str1,str2) %#ok<DEFNU>

    %ApmB = '([0-9\.])+pm([0-9\.]+)';

    %FIX Daniel
    ApmB = '([0-9\.\-e])+pm([0-9\.\-e]+)';

    match_expr      = ApmB;
    m1              = eval(regexprep(str1,match_expr,'$1'));
    s1              = eval(regexprep(str1,match_expr,'$2'));
    m2              = eval(regexprep(str2,match_expr,'$1'));
    s2              = eval(regexprep(str2,match_expr,'$2'));

    [m,j] = min([m1,m2]);

    if j == 1
        s = s1;
    else
        s = s2;
    end

    str = [num2str(m),'pm',num2str(s)];

end

function str = AorB(str1,str2) %#ok<DEFNU>

    %ApmB = '([0-9\.])+pm([0-9\.]+)';

    %FIX Daniel
    ApmB = '([0-9\.\-e])+pm([0-9\.\-e]+)';

    match_expr      = ApmB;
    m1              = eval(regexprep(str1,match_expr,'$1'));
    s1              = eval(regexprep(str1,match_expr,'$2'));
    m2              = eval(regexprep(str2,match_expr,'$1'));
    s2              = eval(regexprep(str2,match_expr,'$2'));

    m = m1 + m2;

    s = sqrt(s1^2 + s2^2);

    str = [num2str(m),'pm',num2str(s)];
end

function [n,n_sd] = addGeneData(g)

% kieran: 22 july 11

    n = nan;
    n_sd = nan;

    %ApmB = '[0-9\.]+pm[0-9\.]+';
    %FIX Daniel
    ApmB = '([0-9\.\-e])+pm([0-9\.\-e]+)';
    
    tries = 0;

    f_and = @(a,b)AandB(a,b);
    f_or = @(a,b)AorB(a,b);
    
    if ~isempty(g)
        while isnan(n)
            
            tries = tries+1;
            if tries > 1000
                fprintf(1, 'warning: stuck at loop evaluating %s\n', g);
                break
            end

            try 
                match_expr      = ApmB; %'([0-9\.])+pm([0-9\.]+)';
                g_av            = regexprep(g,match_expr,'$1');
                g_sd            = regexprep(g,match_expr,'$2');
                n               = eval(g_av);
                n_sd            = eval(g_sd);

            catch %#ok<CTCH>

                % replace brackets 
                match_expr      = ['\(\s*(',ApmB,')\s*\)']; % FIX 
                %match_expr      = ['\((',ApmB,')\)'];

                replace_expr    = '$1';
                g = regexprep(g,match_expr,replace_expr);

                % replace and
                match_expr      = ['(',ApmB,')\s+and\s+(',ApmB,')']; % FIX 
                %match_expr      = ['(',ApmB,') and (',ApmB,')'];
                replace_expr    = '${f_and($1,$2)}';
                %replace_expr    = '${AandB($1,$2)}';
                g = regexprep(g,match_expr,replace_expr,'once');

                % replace or
                match_expr      = ['(',ApmB,')\s+or\s+(',ApmB,')']; % FIX 
                %match_expr      = ['(',ApmB,') or (',ApmB,')'];
                replace_expr    = '${f_or($1,$2)}';
                %replace_expr    = '${AorB($1,$2)}';
                g = regexprep(g,match_expr,replace_expr,'once');

            end
        end
    end

end

function v_sol = dataToFlux(m,r,r_sd)

% kieran: 21 sep 11

    rev = false(size(m.rxns));

    nR_old = 0;

    % m.lb(~ismember(m.lb,[-inf,0,inf,-1000,1000])) = -inf;
    % m.ub(~ismember(m.ub,[-inf,0,inf,-1000,1000])) = inf;

    v_sol = zeros(size(m.rxns));

    tries = 0;
    
    while sum(~m.rev) > nR_old
        
        %fprintf(1, 'irrev old %d, irrev new %d \n', nR_old, sum(~m.rev));
        
        tries = tries+1;
        if tries > 100
            fprintf(1, 'warning: stuck at loop. irrev old %d, irrev new %d \n', nR_old, sum(~m.rev));
            break
        end

        nR_old = sum(~m.rev);

        % 1. fit to data

        N = m.S;    
        L = m.lb;
        U = m.ub;
        f = zeros(size(m.rxns))';
        b = zeros(size(m.mets));

        for k = 1:length(m.rxns)
            d = r(k);
            s = r_sd(k);

            if ~m.rev(k) && ~isnan(d) && s>0
                [s1,s2] = size(N);
                N(s1+1,k) = 1; N(s1+1,s2+1) = -1; N(s1+1,s2+2) = 1;
                L(s2+1) = 0; L(s2+2) = 0;
                U(s2+1) = inf; U(s2+2) = inf;
                b(s1+1) = d;
                f(s2+1) = - 1/s;
                f(s2+2) = - 1/s;
            end
        end

        [v,fOpt,conv] = easyLP(f,N,b,L,U);

        if conv
            v_sol = v(1:length(m.rxns));
            for k = 1:length(m.rxns)
                if rev(k), v_sol(k) = - v_sol(k); end
            end

            % 2. run FVA

            N = [N; f]; %#ok<AGROW>
            b = [b(:); fOpt];

            for k = 1:length(m.rxns)

                if m.rev(k)

                    f = zeros(size(L));
                    f(k) = -1;

                    [~,fOpt,conv] = easyLP(f,N,b,L,U);

                    if conv && (-fOpt >= 0) % irreversibly forward

                        m.lb(k) = max(m.lb(k),0);
                        m.rev(k) = 0;

                    else
                        f(k) = 1;
                        [~,fOpt,conv] = easyLP(f,N,b,L,U);

                        if conv && abs(fOpt)<=0 % irreversibly backward

                            m.S(:,k) = - m.S(:,k);

                            m.ub(k) = - m.ub(k);
                            m.lb(k) = - m.lb(k);

                            ub = m.ub(k);
                            m.ub(k) = m.lb(k);
                            m.lb(k) = ub;

                            m.lb(k) = max(m.lb(k),0);
                            m.rev(k) = 0;

                            rev(k) = ~rev(k);

                        end
                    end
                end
            end
        end
    end

end

function [v,fOpt,conv] = easyLP(f,a,b,vlb,vub)
%
%easyLP
%
% solves the linear programming problem: 
%   max f'x subject to 
%   a x = b
%   vlb <= x <= vub. 
%
% Usage: [v,fOpt,conv] = easyLP(f,a,b,vlb,vub)
%
%   f           objective coefficient vector
%   a           LHS matrix
%   b           RHS vector
%   vlb         lower bound
%   vub         upper bound
%
%   v           solution vector
%   fOpt        objective value
%   conv        convergence of algorithm [0/1]
%
% the function is a wrapper for on the "solveCobraLP" script provided with
% the COBRA (COnstraint-Based Reconstruction and Analysis) toolbox 
% http://opencobra.s.net/
%
%kieran, 20 april 2010

    % matlab can crash if inputs nan
    if any(isnan(f))||any(any(isnan(a)))||any(isnan(b))...
            ||any(isnan(vlb))||any(isnan(vub))  
        error('nan inputs not allowed');
    end

    % initialize
    v = zeros(size(vlb));
    v = v(:);
    f = full(f(:));
    vlb = vlb(:);
    vub = vub(:);

    % remove any tight contstraints as some solvers require volume > 0
    j1 = (vlb ~= vub); 
    j2 = (vlb == vub);
    v(j2) = vlb(j2);
    b = b(:) - a*v;
    a(:,j2) = []; 
    vlb(j2) = []; 
    vub(j2) = [];
    f0 = f;
    f(j2) = [];
    fOpt = nan;

    % solve
    csense(1:length(b)) = 'E';

    solution = solveCobraLP(...
        struct('A',a,'b',b,'c',f,'lb',vlb,'ub',vub,'osense',-1,'csense',csense));

    % define outputs
    conv = solution.stat == 1;

    if conv
        v0 = solution.full;
        v(j1) = v0;
        fOpt = f0'*v;
    end

end
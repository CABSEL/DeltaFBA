function model_binary = createBinaryUseVariable(model_delNet, epsilon, M_prime, regRxns,regRxnsRatio)
% Function for taking a model that has been created for delta and adding
% use variables for each reaction. We use a big M strategy to set binary
% variables Z1,Z2 and Z3 for every variable (reaction fluxes). We also add
% constraints that will make binary varibale Z1 (active or =1) when the
% reaction delta flux > epsilon. Z2 becomes active when delta flux < epsilon. Z3 
% activity will indicate a no change in delta. We also
% add constraints that prevent reversible reactions to have fluxes in the
% forward and reverse direction.
%
%
% Please refer to example_ishii.m for example on what is needed
%
% Inputs:
% 1. model_delNet - GEM transformed to irreversible model, modified bounds
% and with additional netflux variables added in.
% 2. Epsilon - Threshold for defining up and downregulated reactions
% 3. M_prime - A very large number. Typically 3 order of magnitude larger
% than biggest flux.
% 4. regRxns - a cell array of differentially expressed (can be both up and
% downregulated reactions in any order)
% 5. regRxnsRatio - Corresponding thresholding modifications that could be
% applied specifically to define when a reaction is upregulated of
% downregulated. Eg.- a reaction is changing 2 fold, multiplying the
% epsilon by 2 would make sure that the reaction is only upregulated when
% its higher by a factor of its gene expression ratio. If empty then use a
% numeric vector of 1's, the same length as regRxns.
%
%
% Outputs:
% 1. model_binary - Model with binary descriptors for every reaction. 





model = model_delNet;

eps = epsilon;
M = M_prime;

[num_eqs, num_vars] = size(model.A);
[~,id] = ismember(model.rxns, regRxns);

% for Z1, we use the following equations
% V_i - M*Z1_i > eps - M
% V_i - M*Z1_i < eps

% if Z1_i = 0       if Z1_i = 1
% V_i > eps-M       V_i > eps
% V_i < eps         V_i < eps+M

% for Z2, we use the following equations
% V_i + M*Z2_i < -eps + M
% V_i + M*Z2_i > -eps

% if Z2_i = 0        if Z2_i = 1
% V_i < -eps+M       V_i < -eps
% V_i > -eps         V_i > -eps-M

% To account for Z3, we use the following eq
% V_i + M*Z3_i < M
% V_i - M*Z3_i > -M

% if Z3_i = 0       if Z3_i = 1
% V_i < M           V_i < 0
% V_i > -M          V_i > 0
% This will make V_i = 0 when Z3_i = 1;

% for each variable, we add 3 binary variable and 6 equations 

for i=1:size(model.S,2)
    model.A(num_eqs+1, num_vars+1) = -M;
    model.A(num_eqs+1,i) = 1;
    model.eqtype(num_eqs+1) = '>';
    if id(i)~=0
        model.rhs(num_eqs+1) = (eps*regRxnsRatio(id(i)))-M;
    else
        model.rhs(num_eqs+1) = eps-M;
    end
    model.eqNames{(num_eqs+1),1} = strcat('1_Z1_', model.varNames{i,1});
    model.varNames{(num_vars+1),1} = strcat('Z1_', model.varNames{i,1});
    model.varlb(num_vars+1) = 0;
    model.varub(num_vars+1) = 1;
    model.vtype(num_vars+1) = 'B';
    model.obj(num_vars+1) = 0;
    
    model.A(num_eqs+2, num_vars+1) = -M;
    model.A(num_eqs+2,i) = 1;
    model.eqtype(num_eqs+2) = '<';
    if id(i)~=0
        model.rhs(num_eqs+2) = eps*regRxnsRatio(id(i));
    else
        model.rhs(num_eqs+2) = eps;
    end
    model.eqNames{(num_eqs+2),1} = strcat('2_Z1_', model.varNames{i,1});
    
    model.A(num_eqs+3, num_vars+2) = M;
    model.A(num_eqs+3,i) = 1;
    model.eqtype(num_eqs+3) = '<';
    if id(i)~=0
        model.rhs(num_eqs+3) = (-eps*(1/regRxnsRatio(id(i))))+M;
    else
        model.rhs(num_eqs+3) = -eps+M;
    end
    
    model.eqNames{(num_eqs+3),1} = strcat('1_Z2_', model.varNames{i,1});
    model.varNames{(num_vars+2),1} = strcat('Z2_', model.varNames{i,1});
    model.varlb(num_vars+2) = 0;
    model.varub(num_vars+2) = 1;
    model.vtype(num_vars+2) = 'B';
    model.obj(num_vars+2) = 0;
    
    model.A(num_eqs+4, num_vars+2) = M;
    model.A(num_eqs+4,i) = 1;
    model.eqtype(num_eqs+4) = '>';
    if id(i)~=0
        model.rhs(num_eqs+4) = -eps*(1/regRxnsRatio(id(i)));
    else
        model.rhs(num_eqs+4) = -eps;
    end
    model.eqNames{(num_eqs+4),1} = strcat('2_Z2_', model.varNames{i,1});
          
    model.A(num_eqs+5, num_vars+3) = M;
    model.A(num_eqs+5,i) = 1;
    model.eqtype(num_eqs+5) = '<';
    model.rhs(num_eqs+5) = M;
    model.eqNames{(num_eqs+5),1} = strcat('1_Z3_', model.varNames{i,1});
    model.varNames{(num_vars+3),1} = strcat('Z3_', model.varNames{i,1});
    model.varlb(num_vars+3) = 0;
    model.varub(num_vars+3) = 1;
    model.vtype(num_vars+3) = 'B';
    model.obj(num_vars+3) = 0;
    
    model.A(num_eqs+6, num_vars+3) = -M;
    model.A(num_eqs+6,i) = 1;
    model.eqtype(num_eqs+6) = '>';
    model.rhs(num_eqs+6) = -M;
    model.eqNames{(num_eqs+6),1} = strcat('2_Z3_', model.varNames{i,1});
    
    
    [num_eqs, num_vars] = size(model.A);

end

num_rev = numel(find(model.rev))/2;
num_irrev = size(model.S,2)-num_rev;
rev_idx = model.match(size(model.S,2)-num_rev+1:end);
ori_rxns = num_rev+num_irrev;

% For every reversible reaction, we add 6 equation to restrict the activity
% of Z1 or Z2 in both forward and reverse

% for Z1                  for Z2
%          
% Z1_f + Z1_b < 1         Z2_f + Z2_b < 1

% And also to prevent that if a forward reaction is active in the positive
% direction, then the backward flux must either be 0 or in the negative
% direction

% Z1_f < Z2_b + Z3_b
% Z2_f < Z1_b + Z3_b
% Z1_b < Z2_f + Z3_f
% Z2_b < Z1_f + Z3_f

for i = 1:num_rev
    
    model.A(num_eqs+1,(ori_rxns+num_irrev+((rev_idx(i)-1)*3)+1)) = 1; %forward Z1
    model.A(num_eqs+1,(ori_rxns+num_irrev+((num_irrev+i-1)*3)+1)) = 1; % Reverse Z1
    model.eqtype(num_eqs+1) = '<';
    model.rhs(num_eqs+1) = 1;
    model.eqNames{(num_eqs+1),1} = strcat('Activity1_', model.varNames{rev_idx(i),1}(1:end-2));
    
    model.A(num_eqs+2,(ori_rxns+num_irrev+((rev_idx(i)-1)*3)+2)) = 1; %forward Z2
    model.A(num_eqs+2,(ori_rxns+num_irrev+((num_irrev+i-1)*3)+2)) = 1; % Reverse Z2
    model.eqtype(num_eqs+2) = '<';
    model.rhs(num_eqs+2) = 1;
    model.eqNames{(num_eqs+2),1} = strcat('Activity2_', model.varNames{rev_idx(i),1}(1:end-2));
    
    model.A(num_eqs+3,(ori_rxns+num_irrev+((rev_idx(i)-1)*3)+1)) = 1; %Z1_f
    model.A(num_eqs+3,(ori_rxns+num_irrev+((num_irrev+i-1)*3)+2)) = -1; %Z2_b
    model.A(num_eqs+3,(ori_rxns+num_irrev+((num_irrev+i-1)*3)+3)) = -1; %Z3_b
    model.eqtype(num_eqs+3) = '<';
    model.rhs(num_eqs+3) = 0;
    model.eqNames{(num_eqs+3),1} = strcat('Activity3_', model.varNames{rev_idx(i),1}(1:end-2));
    
    model.A(num_eqs+4,(ori_rxns+num_irrev+((rev_idx(i)-1)*3)+2)) = 1; %Z2_f
    model.A(num_eqs+4,(ori_rxns+num_irrev+((num_irrev+i-1)*3)+1)) = -1; %Z1_b
    model.A(num_eqs+4,(ori_rxns+num_irrev+((num_irrev+i-1)*3)+3)) = -1; %Z3_b
    model.eqtype(num_eqs+4) = '<';
    model.rhs(num_eqs+4) = 0;
    model.eqNames{(num_eqs+4),1} = strcat('Activity4_', model.varNames{rev_idx(i),1}(1:end-2));
    
    model.A(num_eqs+5,(ori_rxns+num_irrev+((num_irrev+i-1)*3)+1)) = 1; %Z1_b
    model.A(num_eqs+5,(ori_rxns+num_irrev+((rev_idx(i)-1)*3)+2)) = -1; %Z2_f
    model.A(num_eqs+5,(ori_rxns+num_irrev+((rev_idx(i)-1)*3)+3)) = -1; %Z3_f
    model.eqtype(num_eqs+5) = '<';
    model.rhs(num_eqs+5) = 0;
    model.eqNames{(num_eqs+5),1} = strcat('Activity5_', model.varNames{rev_idx(i),1}(1:end-2));
    
    model.A(num_eqs+6,(ori_rxns+num_irrev+((num_irrev+i-1)*3)+2)) = 1; %Z2_b
    model.A(num_eqs+6,(ori_rxns+num_irrev+((rev_idx(i)-1)*3)+1)) = -1; %Z1_f
    model.A(num_eqs+6,(ori_rxns+num_irrev+((rev_idx(i)-1)*3)+3)) = -1; %Z3_f
    model.eqtype(num_eqs+6) = '<';
    model.rhs(num_eqs+6) = 0;
    model.eqNames{(num_eqs+6),1} = strcat('Activity6_', model.varNames{rev_idx(i),1}(1:end-2));
    
    [num_eqs, num_vars] = size(model.A);
end

model_binary = model;
end

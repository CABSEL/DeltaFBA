function model_binary_gene = createBinaryUseVariable_enzyme(m_GPRT, epsilon, M_prime, regGenes, regGeneRatio)
% Function for taking a model that has been created for delta and adding
% use variables for each enzyme usage reaction.We use a big M strategy to set binary
% variables Z1,Z2 and Z3 for every variable (reaction fluxes). We also add
% constraints that will make binary varibale Z1 (active or =1) when the
% reaction delta flux > epsilon. Z2 becomes active when delta flux < epsilon. Z3 
% activity will indicate a no change in delta. We also
% add constraints that prevent reversible reactions to have fluxes in the
% forward and reverse direction.

model = m_GPRT;

eps = epsilon;
M = M_prime;

[num_eqs, num_vars] = size(model.A);
[~,enzyme_usage] = ismember(strcat('u_',model.genes),model.varNames);
[~,sel] = ismember(strcat('u_',regGenes), model.varNames(enzyme_usage));
Generatio = ones(numel(enzyme_usage),1);
Generatio(sel) = regGeneRatio;

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

for i=1:numel(enzyme_usage)
    
    model.A(num_eqs+1, num_vars+1) = -M;
    model.A(num_eqs+1,enzyme_usage(i)) = 1;
    model.eqtype(num_eqs+1) = '>';
    model.rhs(num_eqs+1) = (eps*Generatio(i))-M;
    model.eqNames{(num_eqs+1),1} = strcat('1_Z1_', model.varNames{enzyme_usage(i),1});
    model.varNames{(num_vars+1),1} = strcat('Z1_', model.varNames{enzyme_usage(i),1});
    model.varlb(num_vars+1) = 0;
    model.varub(num_vars+1) = 1;
    model.vtype(num_vars+1) = 'B';
    model.obj(num_vars+1) = 0;
    
    model.A(num_eqs+2, num_vars+1) = -M;
    model.A(num_eqs+2,enzyme_usage(i)) = 1;
    model.eqtype(num_eqs+2) = '<';
    model.rhs(num_eqs+2) = eps*Generatio(i);
    model.eqNames{(num_eqs+2),1} = strcat('2_Z1_', model.varNames{enzyme_usage(i),1});
    
    model.A(num_eqs+3, num_vars+2) = M;
    model.A(num_eqs+3,enzyme_usage(i)) = 1;
    model.eqtype(num_eqs+3) = '<';
    model.rhs(num_eqs+3) = (-eps*(1/Generatio(i)))+M;
    model.eqNames{(num_eqs+3),1} = strcat('1_Z2_', model.varNames{enzyme_usage(i),1});
    model.varNames{(num_vars+2),1} = strcat('Z2_', model.varNames{enzyme_usage(i),1});
    model.varlb(num_vars+2) = 0;
    model.varub(num_vars+2) = 1;
    model.vtype(num_vars+2) = 'B';
    model.obj(num_vars+2) = 0;
    
    model.A(num_eqs+4, num_vars+2) = M;
    model.A(num_eqs+4,enzyme_usage(i)) = 1;
    model.eqtype(num_eqs+4) = '>';
    model.rhs(num_eqs+4) = -eps*(1/Generatio(i));
    model.eqNames{(num_eqs+4),1} = strcat('2_Z2_', model.varNames{enzyme_usage(i),1});
    
    model.A(num_eqs+5, num_vars+3) = M;
    model.A(num_eqs+5,enzyme_usage(i)) = 1;
    model.eqtype(num_eqs+5) = '<';
    model.rhs(num_eqs+5) = M;
    model.eqNames{(num_eqs+5),1} = strcat('1_Z3_', model.varNames{enzyme_usage(i),1});
    model.varNames{(num_vars+3),1} = strcat('Z3_', model.varNames{enzyme_usage(i),1});
    model.varlb(num_vars+3) = 0;
    model.varub(num_vars+3) = 1;
    model.vtype(num_vars+3) = 'B';
    model.obj(num_vars+3) = 0;
    
    model.A(num_eqs+6, num_vars+3) = -M;
    model.A(num_eqs+6,enzyme_usage(i)) = 1;
    model.eqtype(num_eqs+6) = '>';
    model.rhs(num_eqs+6) = -M;
    model.eqNames{(num_eqs+6),1} = strcat('2_Z3_', model.varNames{enzyme_usage(i),1});
    
    
    [num_eqs, num_vars] = size(model.A);
    
end

model_binary_gene = model;
end

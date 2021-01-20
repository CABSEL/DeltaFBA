function model = addNetFluxVariablesNEW(model)
% adds NetFluxes variables to the model and their associated constraints

[num_mets num_rxns] = size(model.S);
[num_constraints,org_num_vars] = size(model.A);

% F_vi=getAllVar(model,{'F'});
% R_vi=getAllVar(model,{'R'});

model.A=[model.A zeros(num_constraints,num_rxns)];
[num_constraints,num_vars] = size(model.A);

tmp_matrix=zeros(num_rxns,num_vars);

for i=1:num_rxns
    % Find the indices of F_x and R_x variables in the varNames for each
    % reaction x
    F_vi=find(ismember(model.varNames,strcat('F_',model.rxns{i})));
    R_vi=find(ismember(model.varNames,strcat('R_',model.rxns{i})));
    
    % If there exist such F_x variable
    if ~isempty(F_vi)
        % Keep track of the this variable
        tmp_matrix(i,F_vi)=1;
        % Increase the number of variables by 1, and set its "ub" equal to
        % the "ub" of the F_x variable
        model.var_ub(org_num_vars+i)=model.var_ub(F_vi);
        % Make the bounds of the old variable F_x wide enough, so that they
        % do not affect the model when we update the NF_x variable.
        model.var_lb(F_vi)=0;
        model.var_ub(F_vi)=100;
    else
        % Set the upper bound of one extra variable to zero
        model.var_ub(org_num_vars+i)=0;
    end
    
    % If there exist such R_x variable
    if ~isempty(R_vi)
        % Keep track of the this variable
        tmp_matrix(i,R_vi)=-1;
        % Increase the number of variables by 1, and set its "lb" equal to
        % the "-ub" of the R_x variable
        model.var_lb(org_num_vars+i)=-1*model.var_ub(R_vi);
        % Make the bounds of the old variable R_x wide enough, so that they
        % do not affect the model when we update the NF_x variable.
        model.var_lb(R_vi)=0;
        model.var_ub(R_vi)=100;
    else
        % Set the lower bound of one extra variable to zero
        model.var_lb(org_num_vars+i)=0;
    end
    
    tmp_matrix(i,org_num_vars+i)=-1;
end

model.A=[model.A;tmp_matrix];
NF_varNames=strcat('NF_',model.rxns);

model.varNames=[model.varNames;NF_varNames];
model.constraintNames=[model.constraintNames;NF_varNames];
model.rhs=[model.rhs;zeros(num_rxns,1)];
model.constraintType=[model.constraintType;cellstr(repmat('=',num_rxns,1))];

model.vartypes=[model.vartypes;cellstr(repmat('C',num_rxns,1))];
model.f(length(model.f)+1:length(model.var_ub),1) = 0;

NF_vi=getAllVar(model,{'NF'});
model.var_lb(NF_vi)=model.lb;
model.var_ub(NF_vi)=model.ub;


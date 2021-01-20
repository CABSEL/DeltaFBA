function [sol]= solveFBACplex(model,sense)
    % model cobra format
    if nargin<2
        sense='maximize';
    end
    prob=Cplex();
    prob.Model.A=model.S;
    prob.Model.lb=model.lb;
    prob.Model.ub=model.ub;
    prob.Model.lhs=model.b;
    prob.Model.rhs=model.b;
    prob.Model.obj=model.c;
    prob.Model.sense=sense;
    sol=prob.solve();
    
end
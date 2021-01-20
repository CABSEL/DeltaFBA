function [MILP] = create1normSolutionMILP(MILPstructure, NFind, Max_Consistency)
% Function for taking an MILPstructure for maximixing consistency and
% create L1 Norms for all Net fluxes or delta V's.


MILP = MILPstructure;
MILP.lb(find(MILP.obj)) = Max_Consistency;
MILP.obj(find(MILP.obj)) = 0;
for i = 1:numel(NFind)
    MILP.varNames{(end+1),1} = strcat('Abs_',MILP.varNames{(NFind(i)),1});
    MILP.lb(end+1) = 0;
    MILP.ub(end+1) = 1000;
    MILP.vtype(end+1) = 'C';
    MILP.obj(end+1) = 0;
    MILP.A(:,end+1) = 0;
    MILP.genconabs(i).resvar = size(MILP.A,2);
    MILP.genconabs(i).argvar = NFind(i);
end

MILP.A(end+1,end+1) = 1;
MILP.A(end,(size(MILP.A,2)-1-numel(NFind)+1):(size(MILP.A,2)-1)) = -1;
MILP.sense(end+1) = '=';
MILP.rhs(end+1) = 0;
MILP.eqNames{end+1,1} = 'Total_Abs_Delta';
MILP.varNames{end+1,1} = 'Total_Abs_Delta';
MILP.lb(end+1) = 0;
MILP.ub(end+1) = 1e7;
MILP.vtype(end+1) = 'C';
MILP.obj(end+1) = 1;

MILP.modelsense = 'min';

end

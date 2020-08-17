function model_GPRT = createGPRTmodel(model_delNet, IrGPRs)
% Function for taking a model that has been created for delta and adding
% enzyme usage variable. Single enzyme reactions will have only one row and
% column added (-1 for reaction; 1 for continuour enzyme usage variable);
% promiscuous enzymes (one enzyme, multiple reactions) will have 1 row and
% column added (-1 for all reactions and 1 for enzyme usage variable);
% isozymes (multiple enzymes to 1 reaction) will have multiple columns and
% multiple rows added for each enzyme to reaction. Isozymes get special
% attention as 1 reaction is now split into multiple. For each split, we
% will add another row to sum all these reaction to its parents. Enzyme
% complexes will have multiple rows but only 1 column


gpr = IrGPRs;
model = model_delNet;

[num_eqs, num_vars] = size(model.A);
lb = min(model.varlb(1:size(gpr,2)));
ub = max(model.varub(1:size(gpr,2)));

for i=1:size(gpr,2) %every gpr for ir model
    if length(gpr{1,i})==1 %complexes or single enzymes
        for j=1:length(gpr{1,i}{1,1})
            k = find(~cellfun(@isempty,regexp(model.varNames,strcat('u_',gpr{1,i}{1,1}{j},'_\d'))));
            model.A(num_eqs+1, num_vars+1) = 1;
            model.A(num_eqs+1,i) = -1;
            model.eqtype(num_eqs+1) = '=';
            model.rhs(num_eqs+1) = 0;
            model.eqNames{(num_eqs+1),1} = strcat('e_', gpr{1,i}{1,1}{j},'_',model.rxns{i,1});
            model.varNames{(num_vars+1),1} = strcat('u_', gpr{1,i}{1,1}{j},'_',num2str(numel(k)+1));
            model.varlb(num_vars+1) = lb;
            model.varub(num_vars+1) = ub;
            model.vtype(num_vars+1) = 'C';
            model.obj(num_vars+1) = 0;
            [num_eqs, num_vars] = size(model.A);
            
        end
    elseif length(gpr{1,i})>1 %isozymes
        for z=1:numel(gpr{1,i})
            model.A(num_eqs, num_vars+1) = 0;
            model.varNames{(num_vars+1),1} = strcat('Split_',num2str(z),'_',model.rxns{i,1});
            model.varlb(num_vars+1) = lb;
            model.varub(num_vars+1) = ub;
            model.vtype(num_vars+1) = 'C';
            model.obj(num_vars+1) = 0;
            [num_eqs, num_vars] = size(model.A);
            if iscell(gpr{1,i}{1,z})
                for j=1:length(gpr{1,i}{1,z})
                    k = find(~cellfun(@isempty,regexp(model.varNames,strcat('u_',gpr{1,i}{1,z}{j},'_\d'))));
                    model.A(num_eqs+1, num_vars+1) = 1;
                    model.A(num_eqs+1, find(ismember(model.varNames, strcat('Split_',num2str(z),'_',model.rxns{i,1})))) = -1;
                    model.eqtype(num_eqs+1) = '=';
                    model.rhs(num_eqs+1) = 0;
                    model.eqNames{(num_eqs+1),1} = strcat('e_', gpr{1,i}{1,z}{j},'_',model.rxns{i,1});
                    model.varNames{(num_vars+1),1} = strcat('u_', gpr{1,i}{1,z}{j},'_',num2str(numel(k)+1));
                    model.varlb(num_vars+1) = lb;
                    model.varub(num_vars+1) = ub;
                    model.vtype(num_vars+1) = 'C';
                    model.obj(num_vars+1) = 0;
                    [num_eqs, num_vars] = size(model.A);
                    
                end
            else
                k = find(~cellfun(@isempty,regexp(model.varNames,strcat('u_', gpr{1,i}{z},'_\d'))));
                
                model.A(num_eqs+1, num_vars+1) = 1;
                model.A(num_eqs+1, find(ismember(model.varNames, strcat('Split_',num2str(z),'_',model.rxns{i,1})))) = -1;
                model.eqtype(num_eqs+1) = '=';
                model.rhs(num_eqs+1) = 0;
                model.eqNames{(num_eqs+1),1} = strcat('e_', gpr{1,i}{z},'_',model.rxns{i,1});
                model.varNames{(num_vars+1),1} = strcat('u_', gpr{1,i}{z},'_',num2str(numel(k)+1));
                model.varlb(num_vars+1) = lb;
                model.varub(num_vars+1) = ub;
                model.vtype(num_vars+1) = 'C';
                model.obj(num_vars+1) = 0;
                [num_eqs, num_vars] = size(model.A);
                
            end
        end
        model.A(num_eqs+1, find(~cellfun(@isempty,regexp(model.varNames,strcat('Split_\d','_',model.rxns{i,1}))))) = 1;
        model.A(num_eqs+1, i) = -1;
        model.eqtype(num_eqs+1) = '=';
        model.rhs(num_eqs+1) = 0;
        model.eqNames{(num_eqs+1),1} = strcat('Total_Split_', model.rxns{i,1});
        [num_eqs, num_vars] = size(model.A);
    else
        [num_eqs, num_vars] = size(model.A);
    end
end




for i = 1:numel(model.genes)
        ids = find(~cellfun(@isempty,regexp(model.varNames,strcat('u_', model.genes{i,1},'_\d'))));
    
        model.A(num_eqs+1, num_vars+1) = 1;
        model.A(num_eqs+1,ids) = -1;
        model.eqtype(num_eqs+1) = '=';
        model.rhs(num_eqs+1) = 0;
        model.eqNames{(num_eqs+1),1} = strcat('e_',model.genes{i,1});
        model.varNames{(num_vars+1),1} = strcat('u_',model.genes{i,1});
        model.varlb(num_vars+1) = lb;
        model.varub(num_vars+1) = ub;
        model.vtype(num_vars+1) = 'C';
        model.obj(num_vars+1) = 0;
        [num_eqs, num_vars] = size(model.A);
    
%     ids = find(~cellfun(@isempty,regexp(model.varNames,strcat('u_', model.genes{i,1},'_\d'))));
%     model.A(num_eqs+1, num_vars+1) = 1;
%     model.eqtype(num_eqs+1) = '=';
%     model.rhs(num_eqs+1) = 0;
%     model.eqNames{(num_eqs+1),1} = strcat('e_',model.genes{i,1});
%     model.varNames{(num_vars+1),1} = strcat('u_',model.genes{i,1});
%     model.varlb(num_vars+1) = lb;
%     model.varub(num_vars+1) = ub;
%     model.vtype(num_vars+1) = 'C';
%     model.obj(num_vars+1) = 0;
%     clear ids2
%     for b = 1:numel(ids)
%         ids2(b) = find(model.A(find(model.A(:,ids(b))),:)<0);
%     end
%     ids2 = ids2';
%     strs = erase(erase(model.varNames(ids2),'_f'),'_b');
%     count = 1;
%     for b= 1:numel(ids2)
%         str = strs(b);
%         ind = find(ismember(strs,str));
%         if (numel(ind)>1 && ind(1)>=b)
%             model.A(num_eqs+1,ids(ind(1))) = -1;
%             
%             model.A(end+1,ids(ind(1))) = 1;
%             model.A(end,ids(ind(2))) = -1;
%             model.eqtype(end+1) = '=';
%             model.rhs(end+1) = 0;
%             model.eqNames{(end+1),1} = strcat('e_',model.genes{i,1},'_Rev_', num2str(count));
%             count = count+1;
%         else
%             model.A(num_eqs+1,ids(ind(1))) = -1;
%         end
%     end
%     
%     [num_eqs, num_vars] = size(model.A);
    
end

model_GPRT = model;
end

function GPRs = GPR2cell(model)
% This function is used for taking the model's GPR rules in grRules
% structure. The rules are structured as 'or' / 'and'. 342 344 507 572 750
% 1215
GPRs = {}; 

for i=1:numel(model.rxns)
    if length(model.grRules{i})>0
        checkOR = isempty(strfind(model.grRules{i},'or'));
        %checkOK will be 0 when there is or and 1 when not
        
        if checkOR==0
            split = strtrim(regexp(model.grRules{i},'or','split'));
            store = {};
            for j=1:length(split)
                checkBracket = isempty(strfind(split{j},'('));
                    if checkBracket == 0
                        splitAND = strtrim(regexp(split{j},'and','split'));
                        splitAND = strrep(splitAND,'(','');
                        splitAND = strrep(splitAND,')','');
                        store(end+1) = {splitAND};
                        GPRs{i} = store;
                    else
                        store(end+1) = {split{j}};
                        GPRs{i} = store;
                    end
            end
        elseif checkOR==1
            checkAND = isempty(strfind(model.grRules{i},'and'));
            if checkAND==0
                splitAND2 = strtrim(regexp(model.grRules{i},'and','split'));
                splitAND2 = strrep(splitAND2,'(','');
                splitAND2 = strrep(splitAND2,')','');
                GPRs{i} = {splitAND2};
            else
                single = strtrim(model.grRules{i});
                single=strrep(single,'','');
                single=strrep(single,')','');
                test = {single};
                GPRs{i} = {test};
            end
        end
    end
end
end

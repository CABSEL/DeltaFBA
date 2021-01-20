function [parsedGPR,store] = parsingGPRs(model)
    
    %     warning off all
    store=[];
    parsedGPR = containers.Map();
    
    for i = 1:length(model.rxns)
        if length(model.grRules{i}) > 1
            % this script is for checking whether grRule does not contian
            % any or  rule
            flagOR=isempty(strfind(model.grRules{i},'or'));
            if flagOR==0
                [parsing] = strtrim(regexp(model.grRules{i},'or','split'));
                store_AND={};
                for p=1:numel(parsing)
                    % see bracket is available
                    flag=isempty(strfind(parsing{p},'('));
                    if flag==0
                        [parsing1] = strtrim(regexp(parsing{p},'or','split'));
                        
                        if ~isequal(parsing1{1},parsing{p})
                            disp('or coming after or')
                            store=[store;i]; % its store indexes of reactions in which or comming after or.
                        else
                            parseAND = strtrim(regexp(parsing{p},'and','split'));
                            parseAND=strrep(parseAND,'(','');
                            parseAND=strrep(parseAND,')','');
                            store_AND{end+1}=parseAND;
                            
                            parsedGPR(model.rxns{i}) = store_AND;
                        end
                    else
                        % means no bracket
                        store_AND{end+1}={parsing{p}};
                        parsedGPR(model.rxns{i}) = store_AND;
                    end
                    
                end
            elseif flagOR==1;
                flagAND=isempty(strfind(model.grRules{i},'and'));
                if flagAND==0% when rule has only and no or
                    parseAND = strtrim(regexp(model.grRules{i},'and','split'));
                    parseAND=strrep(parseAND,'(','');
                    parseAND=strrep(parseAND,')','');
                    
                    parsedGPR(model.rxns{i}) = {parseAND};
                else
                    singleParse=strtrim(model.grRules{i});
                    singleParse=strrep(singleParse,'(','');
                    singleParse=strrep(singleParse,')','');
                    parsedGPR(model.rxns{i})={{singleParse}};
                end
            end
        end
        
    end
end

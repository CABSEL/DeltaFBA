function [min,max] = parfor_runTMMBig(cplex,varInd,TimeLimitSecs,file,d_num)

% cplex model
% varIdx
% time for cplex
% save file 
% d_num devide in sub groups

cplex.Param.timelimit.Cur = TimeLimitSecs;
min=nan(numel(varInd),1);
max=nan(numel(varInd),1);
xInt=floor(numel(varInd)/d_num);
f_zero=zeros(numel(cplex.Model.obj),1);
for count=1:xInt
%     min=nan(d_num,1);
%     max=nan(d_num,1);
    count
    parfor i=(count-1)*d_num+1:count*d_num
        m=cplex;
        m.Model.obj=f_zero;
        m.Model.obj(varInd(i))=1;
        m.Model.sense='maximize';
        m.Param.timelimit.Cur = TimeLimitSecs;
        sol=m.solve();
        if isfield(sol,'objval')
            max(i) = sol.objval;
        else
            max(i) = nan;
        end
        m.Model.sense='minimize';
        m.Param.timelimit.Cur = TimeLimitSecs;
        sol=m.solve();
        if isfield(sol,'objval')
            min(i) = sol.objval;
        else
            min(i) = nan;
        end
    end
    
    mm=[min max];
    save(file,'mm')
end

%% last part

    
    parfor i=count*d_num+1:numel(varInd)
        m=cplex;
        m.Model.obj=f_zero;
        m.Model.obj(varInd(i))=1;
        m.Model.sense='maximize';
        m.Param.timelimit.Cur = TimeLimitSecs;
        sol=m.solve();
        if isfield(sol,'objval')
            max(i) = sol.objval;
        else
            max(i) = nan;
        end
        m.Model.sense='minimize';
        m.Param.timelimit.Cur = TimeLimitSecs;
        sol=m.solve();
        if isfield(sol,'objval')
            min(i) = sol.objval;
        else
            min(i) = nan;
        end
    end
mm=[min max];
save(file,'mm')

end
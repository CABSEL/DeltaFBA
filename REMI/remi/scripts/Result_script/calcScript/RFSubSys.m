function [result]=RFSubSys(mm1,mm2,tol)
% mm1:  minmax for reference
% mm2:  minmax for second condition
mm1((mm1>-tol & mm1<tol))=0;
mm2((mm2>-tol & mm2<tol))=0;
perRedFR=nan(size(mm1,1),1);
for i=1:size(mm1,1)
    
    
    
    if (abs(mm1(i,1))>tol || abs(mm1(i,2))>tol)  && ((mm1(i,2)-mm1(i,1))>tol) % to see whether ther are not fixed
        % check fixed range
        if (mm2(i,2)-mm2(i,1))>(mm1(i,2)-mm1(i,1))
            mm2(i,:)=mm1(i,:); % removing some TFA error
        end
%         perRedFR(i,1)=(1-((mm2(i,2)-mm2(i,1))/(mm1(i,2)-mm1(i,1))))*100;
        perRedFR(i,1)=((mm2(i,2)-mm2(i,1))/(mm1(i,2)-mm1(i,1)));
        
    end
end

perRedFR(perRedFR<tol & perRedFR>-tol)=0; % numerical error
inds=find(~(isnan(perRedFR) | perRedFR<0));
result.perRedFR=perRedFR;
result.notNanIdx=inds;
result.notNanVals=perRedFR(inds);
if numel(result.notNanVals)>0
    result.avg=mean(result.notNanVals);
else
    result.avg=nan;
end
    
    
    
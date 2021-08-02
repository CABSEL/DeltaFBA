function [bdr]=IdentifyBDR(mm,tol)
bi=[];
uni=[];
block=[];
mm((mm>-tol & mm<tol))=0;
for i=1:size(mm,1)
    if (mm(i,1)<-tol & mm(i,2)>tol)
        bi=[bi;i];
    elseif (mm(i,1)<-tol & mm(i,2)<=0)||(mm(i,1)>=0 & mm(i,2)>tol)
        uni=[uni;i];
    else
        block=[block;i];
    end
end
bdr.bi=bi;
bdr.uni=uni;
bdr.block=block;
bdr.mm=mm;
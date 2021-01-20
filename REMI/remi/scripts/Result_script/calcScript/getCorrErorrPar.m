function [storeRes1,storeRes2,storeRes3,storeRes4,storeXX]=getCorrErorrPar(model,rxns,z_matrix,FLUX_LIST1,MEASURED_FLUX1,FLUX_LIST2,MEASURED_FLUX2)
    % this function we used to compute fist solution for each alternative 
    % then we computed Pearson correlation and then percentage error
    [~,NFind]=ismember(strcat('NF_',rxns),model.varNames); % get index of wild type
    [~,PNFind]=ismember(strcat('PERTURB_NF_',rxns),model.varNames); % get index of perturbed state
    
    % this is for storing results 
    storeRes1=nan(1,size(z_matrix,2)); % correlation for wild type
    storeRes2=nan(1,size(z_matrix,2)); % correlation for mutant type
    storeRes3=nan(1,size(z_matrix,2)); % correlation diffrence
    storeRes4=nan(1,size(z_matrix,2)); % this is for the percentage error
    %%
    fid = fopen(FLUX_LIST1); % open the file which contains fluxomics data
    match_list = textscan(fid, '%s');
    fclose(fid);
    ll=numel(match_list{1});
    
    %%
  
    storeXX=nan(ll,4,size(z_matrix,2));
    parfor eId=1:size(z_matrix,2)
        model_tmp=model;
        %objIndex=model_tmp.objIndex1B; % in case of expression 
        %activeIdx=model_tmp.objIndex1B(find(z_matrix(:,eId)>0.99));
        %objIndex=[model_tmp.objIndex1B;model_tmp.metB] % in case of expression and metabolite
        objIndex=[model_tmp.metB] % in case of expression and metabolite
        activeIdx=objIndex(find(z_matrix(:,eId)>0.99));
        model_tmp.var_lb(activeIdx)=1;
        
        f=zeros(numel(model_tmp.f),1);
        f(NFind)=-1;
        f(PNFind)=-1;
        model_tmp.f=f;
        sol=solveTFBAmodel(model_tmp);
        if isempty(sol.x)
            storeRes1(eId)=nan;
            storeRes2(eId)=nan;
            storeRes3(eId)=nan;
            storeRes4(eId)=nan;
        else
            flux = sol.x(NFind);
            [corrCond1,m1,p1]=compareExpSim(rxns,flux,FLUX_LIST1,MEASURED_FLUX1); % this function will compute pearson correlation
            flux = sol.x(PNFind);
            [corrCond2,m2,p2]=compareExpSim(rxns,flux,FLUX_LIST2,MEASURED_FLUX2);
            xx=[m1 m2 p1 p2];
            storeXX(:,:,eId)=xx;
            %% skip nan zero rows
            xx(any(isnan(xx),2),:) = [];
            flux_cut=1e-05;
            xx(xx<flux_cut & xx>-flux_cut)=0;
            exp=xx(:,[1 2]);
            sim=xx(:,[3 4]);
            d_exp=(exp(:,2)-exp(:,1))./(abs(exp(:,2))+abs(exp(:,1)));
            % find zero rows
            d_exp(all(exp==0,2))=0;
            d_sim=(sim(:,2)-sim(:,1))./(abs(sim(:,2))+abs(sim(:,1)));
            % find zero rows
            d_sim(all(sim==0,2))=0;
            % find correlation
            corrDiff = dot(d_sim,d_exp)/(norm(d_sim)*norm(d_exp));
            percentErorr=mean(abs(d_sim-d_exp)); % here we compute percentage error 
            storeRes1(eId)=corrCond1;
            storeRes2(eId)=corrCond2;
            storeRes3(eId)=corrDiff;
            storeRes4(eId)=percentErorr;
        end
    end
end

function findAltCombi(numsol,model,index_updown,path_save)
    
   % numsol: number of solutions
   % model: in TFA format
   % index_updown: indexes of varaiables for which we need alternate
   % combinations
   % path save='loaction of file in which u want to save results'
   %output save: it saves solution as well as model and Z_matrix
   %combinations
   z_matrix=[];
   sol_matrix=[];
   sol=solveTFBAmodel(model);
   z_matrix=[z_matrix sol.x(index_updown)];
   sol_matrix=[sol_matrix sol.x];
   store_obj=[];
   % find alternate solutions
   obj=sol.val;
   store_obj(end+1)=obj;
%    numsol=200;
   [num_cons,num_vars] = size(model.A);
   % add size constraint
%    mat=[zeros(1,num_vars)];
%    mat(index_updown)=1;
%    model.A=[model.A;mat];
%    model.constraintNames{end+1}=['Size_cut'];
%    model.constraintType((end+1)) = {'>'};
%    model.rhs(end+1)=obj-2;
   %
   
   for i=1:numsol
       mat=[zeros(1,num_vars)];
       nonzero_ind=find(sol.x(index_updown)>0.98);
       mat(index_updown(nonzero_ind))=1;
       model.A=[model.A;mat];
       model.constraintNames{end+1}=['integercut' num2str(i)];
       model.constraintType((end+1)) = {'<'};
       model.rhs(end+1)=numel(nonzero_ind)-1;
       sol=solveTFBAmodel(model);
       save(path_save,'z_matrix','store_obj','sol_matrix','model');
       if isempty(sol.x)
           break
       else
           z_matrix=[z_matrix sol.x(index_updown)];
           store_obj(end+1)=sol.val;
           sol_matrix=[sol_matrix sol.x];
       end
   end
   save(path_save,'z_matrix','store_obj','sol_matrix','model');
   
   %
   
end
function [correl_m,m,p]=compareExpSim(rxns,eflux2_flux,FLUX_LIST,MEASURED_FLUX)
    %% this model is beacuse data are prepared for
    
    fid = fopen(FLUX_LIST);
    match_list = textscan(fid, '%s');
    fclose(fid);
    
    % process the predicted fluxes according to the meausred flux - the predicted flux association
    modified_correl = zeros(size(match_list{1}));
    
    for i = 1:length(match_list{1})
        silico_vitro = match_list{1}{i};
     
        lv = 1;
        j = 1;
        op = [];
        val = [];
        
        while j <= length(silico_vitro)
            
            if silico_vitro(j) == '('
                lv = lv + 1;
                j = j + 1;
                
            elseif silico_vitro(j) == ')'
                newval = val(lv);
                lv = lv - 1;
                
                if  length(op) >= lv && op(lv) > 0
                    
                    if  op(lv) == 1
                        val(lv) = val(lv) + newval;
                        
                    elseif op(lv) == 2
                        val(lv) = min([val(lv); newval]);
                        
                    end
                    op(lv) = 0;
                else
                    val(lv) = newval;
                    
                end
                j = j + 1;
                
            elseif strcmp(silico_vitro(j:(j + 3)), '_OR_')
                
                op(lv) = 1;
                j = j + 4;
                
            elseif strcmp(silico_vitro(j:(j + 4)), '_AND_')
                op(lv) = 2;
                j = j + 5;
                
            else
                
                k = strfind(silico_vitro(j:end),')');
                if ~isempty(k)
                    flux_string = silico_vitro(j:(j + k(1) - 2));
                else
                    flux_string = silico_vitro(j:end);
                end
                
                % find and assign eflux value of corresponding flux name
                if strncmp('-', flux_string, 1) == 1 % flux name has a negative value in front of itself
                    str_to_compare = [ flux_string(2:end)];
                    iflux = find(strcmp(str_to_compare, rxns));
                    if ~isempty(iflux) % matched with the flux name
                        newval = -eflux2_flux(iflux); % assign 'negative' value of eflux
                        
                    else
                        newval = inf;
                        
                    end
                else
                    iflux = find(strcmp([flux_string], rxns));
                    if ~isempty(iflux)
                        newval = eflux2_flux(iflux);
                        
                    else
                        newval = inf;
                        
                    end
                end
                
                if  length(op) >= lv && op(lv) > 0
                    if  op(lv) == 1
                        val(lv) = val(lv) + newval; % OR relationship among fluxes => sum
                        
                    elseif op(lv) == 2
                        val(lv) = min([val(lv); newval]); % AND relationhsip among fluxes => minimum
                        
                    end
                    op(lv) = 0;
                    
                else
                    val(lv) = newval;
                end
                
                j = j + length(flux_string);
                
            end
            
        end
        
        modified_correl(i) = val(1);
        
    end
    
    
    %%
    insilico_flux = modified_correl;
    % import measured fluxes
    
    fid = fopen(MEASURED_FLUX);
    measured_flux = fscanf(fid, '%f');
    fclose(fid);
    
    % calculate correl b/w measured flux & predicted flux - except for NaN elements in the measured fluxes
    m = measured_flux(~isnan(measured_flux))/100;
    p = insilico_flux(~isnan(measured_flux));
    correl_m = dot(p,m)/(norm(p)*norm(m));
    m=measured_flux/100;
    p=insilico_flux;
end
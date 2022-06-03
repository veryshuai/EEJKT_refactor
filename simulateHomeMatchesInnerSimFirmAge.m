%%                         simulateHomeMatchesInnerSimFirmAge    

% calculate each firm's time (in periods) in domestic market 

    flr = max(iterH_in.flrlag,t*iterH_in.new_firm(:,t));  
    % floor resets to current year for new firms and count begins from there. 
    age = t*ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1) - flr;    
    iterH_in.flrlag = flr ;
    iterH_in.cumage = cat(2,iterH_in.cumage,age);
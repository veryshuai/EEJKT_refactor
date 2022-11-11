%%                         simulateHomeMatchesInnerSimFirmAge    

% calculate each firm's time (in periods) in domestic market 

function iterH_in = simulateHomeMatchesInnerSimFirmAge(iterH_in,mm)

    flr = max(iterH_in.flrlag, iterH_in.t*iterH_in.new_firm(:,iterH_in.t));  
    % floor resets to current year for new firms and count begins from there. 
    age = iterH_in.t*ones(mm.sim_firm_num_by_prod_succ_type(iterH_in.pt_ndx),1) - flr;    
    % Age in periods. age=0 all year for firms with no shipment previous two years
    iterH_in.flrlag = flr ;
    iterH_in.cumage = cat(2,iterH_in.cumage,age);
    iterH_in.age    = age;
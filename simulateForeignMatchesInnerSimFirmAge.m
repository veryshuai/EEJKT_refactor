%%                         simulateForeignMatchesInnerSimFirmAge    

% calculate each firm's time (in periods) in foreign market 

function [iter_in, age] = simulateForeignMatchesInnerSimFirmAge(iter_in, mm)
    % Calculate time (in periods) in export market
    flr = max(iter_in.flrlag,iter_in.t*iter_in.new_firm(:,iter_in.t)); 
    % floor resets to current period for new exporters.
    age = iter_in.t*ones(mm.sim_firm_num_by_prod_succ_type(iter_in.pt_ndx),1) - flr;     
    % age in periods. age=0 all year for firms with no shipment previous year
    
    iter_in.flrlag    = flr ; % carry floor forward for continuing matches
    iter_in.cumage    = cat(2,iter_in.cumage,age);
end
function [nship_obs,ln_ships,match_count] = match_shpt_cntr(matches,max_match)

% match = mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age w/in yr, firm age] 

ff2_obs = find(matches(:,3) >0); % obs. with positive current year shipments
    nship_obs = size(ff2_obs,1);
    ln_ships = sum(log(matches(ff2_obs,3)));
    
%% construct match count vector
    
    if sum(matches(:,3))>0       
        % add up matches with positive shipments, firm by firm:
        num_match = sum(dummyvar(matches(:,1)).*(matches(:,3)>0)); 
        % count firms with 0 matches, 1 match, etc.
%  to drop firms with more than max_match matches:
%         temp = sum(dummyvar([num_match'+1;max_match+1]))';
%         temp = min(temp,max_match+1); 
%         match_count = temp(2:max_match+1);
%         match_count(max_match) = match_count(max_match) - 1;
        
%  alternatively, to topcode firms with more than max_match matches:
        num_match = min(num_match,max_match);
        temp = sum(dummyvar([num_match'+1;max_match+1]))';
        match_count = temp(2:max_match+1);
        match_count(max_match) = match_count(max_match) - 1;
    else
        match_count = zeros(max_match,1);
    end
end
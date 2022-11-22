function [nship_obs,ln_ships,match_count] = match_shpt_cntr(matches,max_match)

% This function is called from simulateForeignMatchesInnerMoments, which
% passes the arguments iter_in.mat_yr_sales and mm.max_match

% match = mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age, firm age] 
% match_count is a vector that counts the number of matches for each active firm

%% find number of shipments for each match, take logs and sum

ff1_obs = find((floor(matches(:,1))-matches(:,1)==0).*(matches(:,3) >0)); % firms with integer IDs
ff2_obs = find((floor(matches(:,1))-matches(:,1)~=0).*(matches(:,3) >0)); % new firms with +0.5 IDs 

nship_obs = length(ff1_obs)+length(ff2_obs); % number of matches w/ shipments > 0
ln_ships  = sum(log(matches(ff1_obs,3))) + sum(log(matches(ff2_obs,3)));
    
%% construct match count vector for each firm_ID
match_count = double.empty(0,1);

      if length(unique(ff1_obs))>0
         incumb_match = sum(dummyvar(matches(ff1_obs,1)))';     
         match_count = sortrows(incumb_match(incumb_match>0));
     end
     
     if length(unique(ff2_obs))>0
        newfirm_match = sum(dummyvar(matches(ff2_obs,1)))';   
        match_count = [match_count; sortrows(newfirm_match(newfirm_match>0))];
        match_count = sortrows(match_count);
     end
     
end
function [nship_obs,ln_ships,match_count] = match_shpt_cntr(matches,max_match)

% match = mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age w/in yr, firm age] 

% if size(matches,1) > 10
%     'pause here'
% end

ff1_obs = find((floor(matches(:,1))-matches(:,1)==0).*(matches(:,3) >0)); % firms with integer IDs
ff2_obs = find((floor(matches(:,1))-matches(:,1)~=0).*(matches(:,3) >0)); % new firms with +0.5 IDs 

nship_obs = length(ff1_obs)+length(ff2_obs);
ln_ships  = sum(log(matches(ff1_obs,3))) + sum(log(matches(ff2_obs,3)));
    
%% construct match count vector
match_count = double.empty(0,2);

     if length(unique(ff1_obs))>1
        incumb_match = [1:max(matches(ff1_obs,1));sum(dummyvar(matches(ff1_obs,1)))]';     
        match_count = [match_count; incumb_match];
     end
     if length(unique(ff1_obs))==1
        incumb_match = [matches(ff1_obs,1);sum(dummyvar(matches(ff1_obs,1)))]'; 
        match_count = [match_count; incumb_match];
     end
     
     if length(unique(ff2_obs))>1
        newfirm_match = [1:2*max(matches(ff2_obs,1));sum(dummyvar(2*matches(ff2_obs,1)))]';   
        match_count = [match_count; [newfirm_match(:,1)./2, newfirm_match(:,2)]];
     end
     if length(unique(ff2_obs))==1
        newfirm_match = [matches(ff2_obs,1);sum(dummyvar(matches(ff2_obs,1)))]'; 
        match_count = [match_count; [newfirm_match(:,1)./2, newfirm_match(:,2)]];
     end
    
   % match count: [firm ID, number of matches in year]
        match_count = sortrows(match_count(match_count(:,2)>0,:),1);


%  to drop firms with more than max_match matches:
%         temp = sum(dummyvar([num_match'+1;max_match+1]))';
%         temp = min(temp,max_match+1); 
%         match_count = temp(2:max_match+1);
%         match_count(max_match) = match_count(max_match) - 1;
        
%  alternatively, to topcode firms with more than max_match matches:
%         num_match = min(num_match,max_match);
%         temp = sum(dummyvar([num_match'+1;max_match+1]))';
%         match_count = temp(2:max_match+1);
%         match_count(max_match) = match_count(max_match) - 1;

end
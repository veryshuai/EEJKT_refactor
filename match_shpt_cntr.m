function [nship_obs,ln_ships,match_count,match_countD] = match_shpt_cntr(iter_in,mm)

% This function is called from simulateForeignMatchesInnerMoments, which
% passes the arguments iter_in.mat_yr_sales and mm.max_match

% match = mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age, firm age] 
% match_count is a vector that counts the number of matches for each active firm

yr_tlag   = iter_in.t-mm.pd_per_yr;
matches  = iter_in.mat_yr_sales;
duds     = iter_in.cur_duds(:,yr_tlag+1:iter_in.t);

%% find number of shipments for each match, take logs and sum
 % firms with integer IDs & shipments>0
ff1_obs = find((floor(matches(:,1))-matches(:,1)==0).*(matches(:,3) >0));
 % new firms with +0.5 IDs & shipments>0
ff2_obs = find((floor(matches(:,1))-matches(:,1)~=0).*(matches(:,3) >0));
% number of matches w/ shipments > 0
nship_obs = length(ff1_obs)+length(ff2_obs); 
ln_ships  = sum(log(matches(ff1_obs,3))) + sum(log(matches(ff2_obs,3)));
    
%% match counts by frequency across firm_ID, excluding duds
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
     
%% stack dud matches with others to get alternative match count

% Note: this dud count will miss duds at new firms during their first year

dud_count   = [(1:size(duds,1))',sum(duds,2)];
ff_dudfinder = find(dud_count(:,2)>0);
dud_count(ff_dudfinder,:);

dud_matches = double.empty(0,7);
for ii = 1:size(duds,1)
    Ndud = dud_count(ii,2);
  if Ndud >0  
  dud_matches = [dud_matches; [dud_count(ii,1)*ones(Ndud,1), ones(Ndud,6)]];
  end
end
assert(sum(dud_count(:,2)) == size(dud_matches,1))

matchesD = [matches;dud_matches];
match_countD = double.empty(0,1);

ff1D_obs = find((floor(matchesD(:,1))-matchesD(:,1)==0).*(matchesD(:,3) >0));
 % new firms with +0.5 IDs & shipments>0
ff2D_obs = find((floor(matchesD(:,1))-matchesD(:,1)~=0).*(matchesD(:,3) >0));
% number of matches w/ shipments > 0

      if length(unique(ff1D_obs))>0        
         incumb_matchD = sum(dummyvar(matchesD(ff1D_obs,1)))';     
         match_countD = sortrows(incumb_matchD(incumb_match>0));
     end
     
     if length(unique(ff2D_obs))>0
        newfirm_matchD = sum(dummyvar(matchesD(ff2D_obs,1)))';   
        match_countD = [match_countD; sortrows(newfirm_matchD(newfirm_match>0))];
        match_countD = sortrows(match_countD);
     end
     
end
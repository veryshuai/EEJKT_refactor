function [nship_obs,ln_ships,match_count,match_countD,dud_matches] = match_shpt_cntr(iter_in,mm)

% This function is called from simulateForeignMatchesInnerMoments, which
% passes the arguments iter_in.mat_yr_sales and mm.max_match

% match = mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age, firm age] 
% match_count is a vector that counts the number of matches for each active firm

yr_tlag      = iter_in.t-mm.pd_per_yr;
matches      = iter_in.mat_yr_sales;
all_duds     = [(1:size(iter_in.cur_duds,1))',iter_in.cur_duds(:,yr_tlag+1:iter_in.t)];
all_firm_age = iter_in.cumage(:,end);

%% split firmID into multiple dud rows if hotel room turns over during the year

% A firm is new until it makes a match. 
[aa,bb] = find(iter_in.new_firm(:,yr_tlag+1:iter_in.t) - iter_in.new_firm(:,yr_tlag:iter_in.t-1)==-1);
if ~isempty(aa)
    newfirm_duds = zeros(length(aa),mm.pd_per_yr+1);
    newfirm_age  = zeros(length(aa),1);
    incumb_duds = all_duds;
    for ii=1:length(aa)
        newfirm_duds(ii,:) = [aa(ii)+0.5,zeros(1,bb(ii)-1),all_duds(aa(ii),bb(ii):mm.pd_per_yr)];
        incumb_duds(aa(ii),:) = [incumb_duds(aa(ii),1:bb(ii)),zeros(1,mm.pd_per_yr - bb(ii)+1)]; 
        newfirm_age(ii) = mm.pd_per_yr - bb(ii); 
    end
    all_duds = [incumb_duds;newfirm_duds];
    all_firm_age = [all_firm_age;newfirm_age];
end

%% find number of shipments for each match, take logs and sum
 % firms with integer IDs & shipments>0:
ff1_obs = find((floor(matches(:,1))-matches(:,1)==0).*(matches(:,3) >0));
 % new firms with +0.5 IDs & shipments>0:
ff2_obs = find((floor(matches(:,1))-matches(:,1)~=0).*(matches(:,3) >0));
 % number of matches w/ shipments > 0
nship_obs = length(ff1_obs)+length(ff2_obs); 
 % sum of log shipment counts across matches
ln_ships  = sum(log(matches(ff1_obs,3))) + sum(log(matches(ff2_obs,3)));
    
%% match counts by frequency across firm_ID, excluding duds
match_count = double.empty(0,1);

      if ~isempty(ff1_obs) 
        incumb_match = sum(dummyvar(matches(ff1_obs,1)))';     
        match_count = sortrows(incumb_match(incumb_match>0));
      end
      
     if ~isempty(ff2_obs)
        newfirm_match = sum(dummyvar(2*matches(ff2_obs,1)))';   
        match_count = [match_count; sortrows(newfirm_match(newfirm_match>0))];
        match_count = sortrows(match_count);
     end

     
%% stack dud matches with others to get alternative match count

dud_count   = [all_duds(:,1),sum(all_duds(:,2:end),2),all_firm_age];

dud_matches = double.empty(0,7);

% impute dud revenue using profit function evaluated at Zcut, macro state, and 1 shipment
scale     = mm.scale_f;
macro_shk = mm.X_f(iter_in.macro_state_f(iter_in.t));
dud_rev   = exp(scale + (mm.eta-1)*mm.Phi(mm.pt_type(iter_in.pt_ndx,1)) + macro_shk).*exp(iter_in.Zcut_eoy); 

% create records of dud matches and stack them
for ii = 1:size(all_duds,1)
  Ndud = dud_count(ii,2);
  if Ndud >0  
  % dud matches: firm_ID, shipment=1, sale=dud_rev, bop Z = eop Z = match_age = 0, and actual firm age     
  dud_matches = [dud_matches; [dud_count(ii,1)*ones(Ndud,1), ones(Ndud,1).*[dud_rev 1 0 0 0 dud_count(ii,3)] ] ];
  end
end
assert(sum(dud_count(:,2)) == size(dud_matches,1))

% stack dud matches with others to crate alternative set of match records
matchesD = [matches;dud_matches];
% matchesD: [firmID, revenue, shipments, boy Z, eoy Z, match age, firm age]

% initialize vector that will count firms in each match count category
match_countD = double.empty(0,1);

% incumbent firms with shipments>0
ff1D_obs = find((floor(matchesD(:,1))-matchesD(:,1)==0).*(matchesD(:,3) >0));
% new firms with +0.5 IDs & shipments>0
ff2D_obs = find((floor(matchesD(:,1))-matchesD(:,1)~=0).*(matchesD(:,3) >0));

% number of matches w/ shipments > 0
      
    if ~isempty(ff1D_obs) % incumbent firm match counts       
        % Note: had to double match ID bec. dummyvar only works with integers
        % All incumbents become even firmIDs, all new firms become odd IDs
        incumb_matchD = sum(dummyvar(2*matchesD(ff1D_obs,1)))';     
        match_countD = sortrows(incumb_matchD(incumb_matchD>0));
     end
     
     if~isempty(ff2D_obs) % new firm counts stacked with incumbent counts
        newfirm_matchD = sum(dummyvar(2*matchesD(ff2D_obs,1)))';   
        match_countD = [match_countD; sortrows(newfirm_matchD(newfirm_matchD>0))];
        match_countD = sortrows(match_countD);
     end
     
end
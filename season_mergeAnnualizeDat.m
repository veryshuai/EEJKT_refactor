function [mat_yr_sales,firm_yr_sales,iterX_in] = season_mergeAnnualizeDat(all_seas, som_seas, mm, mat_cols,iterX_in)
 
 % all_seas and som_seas: 
 %  (1) t, (2) season, (3) year, (4) initial state, (5) exporter id, (6) ending state,
 %  (7) match revenue,(8) #shipments,(9) exporter age (#periods), (10) match age w/in year

 t = iterX_in.t;
  
 a = find(sum(all_seas,2)>0); % non-zero rows of all_seas
 s = find(sum(som_seas,2)>0); % non-zero rows of som_seas
 try
 fndx_a = zeros(size(a,1),mm.pd_per_yr);
 fndx_s = zeros(size(s,1),mm.pd_per_yr);  
  for ss=1:mm.pd_per_yr 
   fndx_a(:,ss) = all_seas(a,ss*mat_cols-5); % each season, firm ID is 6th col. from end, if active
   fndx_s(:,ss) = som_seas(s,ss*mat_cols-5);
  end
 firm_a = max(fndx_a,[],2); % get firm IDs for each row of all_seas
 firm_s = max(fndx_s,[],2); % get firm IDs for each row of som_seas
 catch
     'problem in season_mergeAnnualizeDat lines 12-19'
 end
%% aggregate match-specific sales and shipments across seasons
 pick_sales     = kron(ones(mm.pd_per_yr,1),[zeros(6,2);eye(2);zeros(2,2)]); % to pick off sales and shipments from all seasons
 mat_yr_sales_a = all_seas(a,:)*pick_sales; % add up within-yr. sales for matches that are active at year's end
 mat_yr_sales_s = som_seas(s,:)*pick_sales; % add up within-yr. sales for matches that die before year's end
 

%% create indicator for season in which firm flips
  firm_age = iterX_in.cumage(:,t-mm.pd_per_yr+1:t); % age in periods, current year  
 if t > mm.pd_per_yr
  firm_age_lag = iterX_in.cumage(:,t-mm.pd_per_yr:t-1); % age in periods, last year 
 else
  firm_age_lag = zeros(size(firm_age,1),1);
 end
 flip_seas = (firm_age-firm_age_lag)<0; % marks season in which firm flips
 flip_firm = sum(flip_seas,2)>0; % equals 1 when age drops during year
 
 % iterX_in.flip_ndx gives the season (1-12), if any, after a firm exit.
 % (It is 0 for non-flippers.)
 iterX_in.flip_ndx = zeros(length(flip_seas),1);
 [rr,cc] = find(flip_seas);
 for jj=1:length(rr)
 iterX_in.flip_ndx(rr(jj)) = cc(jj);
 end

 %% Find firm IDs, match age, and firm age for populated rows, by season
 for ss=1:mm.pd_per_yr 

% temp0 = sortrows([[firm_a,age_aa];[firm_s,age_ss]],[1 2]); % stack all match obs. on firm, age and sort
     
 % assign max of firm age for each firm_ID to all matches for continuing firms
 
 %% find within-year entry month for each match 
 z_a = zeros(length(a),mm.pd_per_yr);
 z_s = zeros(length(s),mm.pd_per_yr);
 for ss = 1:mm.pd_per_yr
     z_a(:,ss) = (all_seas(a,10*ss-6)>0)*ss;
     z_s(:,ss) = (som_seas(s,10*ss-6)>0)*ss;
 end
 

%% stack all matches in mat_yr_sales 
% find max age (in periods) within year for each match
matage_pick = zeros(mm.pd_per_yr*mat_cols,mat_cols);
firmage_pick = zeros(mm.pd_per_yr*mat_cols,mat_cols);
  for kk=1:mm.pd_per_yr
      matage_pick(kk*mat_cols,kk) = 1;
      firmage_pick(kk*mat_cols-1,kk) = 1; 
  end
  
if size(z_a,1) + size(z_s,1) > 0
    mat_age_a = max(all_seas(a,:)*matage_pick,[],2);
    mat_age_s = max(som_seas(s,:)*matage_pick,[],2);
    firm_age_a = all_seas(a,mat_cols*mm.pd_per_yr-1);  % firm age at end of year (ongoing matches)
    firm_age_s = max(som_seas(s,:)*firmage_pick,[],2); % firm age in match's final period (dying matches)
     
 mat_yr_sales =...
     sortrows([[firm_a,mat_yr_sales_a,all_seas(a,4),all_seas(a,mat_cols*mm.pd_per_yr-4),...
                mat_age_a ,firm_age_a];...
               [firm_s,mat_yr_sales_s,som_seas(s,4),som_seas(s,mat_cols*mm.pd_per_yr-4),...
                mat_age_s,firm_age_s]]);
          
% mat_yr_sales: [firm ID, w/in yr sales, w/in yr shpmts., boy Z, eoy Z, w/in yr match age, firm age]


%% Distinguish match of entering firms from matches of exiting firms  
%  (Add 0.5 to the firm IDs associated with the former matches.)
  F_active = unique(mat_yr_sales(:,1)); 
  entry_month = iterX_in.flip_ndx;
  % OLD: entry_month = sortrows([[firm_a,entry_seas_a];[firm_s,entry_seas_s]],1);

try
  if sum(flip_firm)>0
  for j=1:length(F_active)
     j_type = mat_yr_sales(:,1)==F_active(j); 
     % adjust firm ID to distinguish post-flip matches
     mat_yr_sales(j_type,1) =...
              mat_yr_sales(j_type,1) + 0.5.*flip_firm(F_active(j))...
           .*(entry_month(F_active(j))>iterX_in.flip_ndx(F_active(j)));

  end
  end
catch
    'problem in season_mergeAnnualizeDat lines 118-122'
end
%% Sum over matches for each active firm

 % tenure  = mat_yr_sales(:,7)>mm.pd_per_yr; % matches with firm age > 1 year
 
  firm_yr_sales = double.empty(0,4);    
  F_active2 = unique(mat_yr_sales(:,1)); % Now may include IDs with 0.5 added
  for ff=1:length(F_active2)
      % incumbents
    select = find(mat_yr_sales(:,1)==F_active2(ff));
    if size(select,1)>0
     firm_yr_sales = [firm_yr_sales;...
       [F_active2(ff),sum(mat_yr_sales(select,2:3),1),max(mat_yr_sales(select,7))] ];
    end
  end
   

else
  firm_yr_sales = double.empty(0,4);
  mat_yr_sales  = double.empty(0,7);
end

% fprintf('\ryear and pt_ndx =%4.0f %4.0f\n', [iterX_in.year, iterX_in.pt_ndx] )

% if t==540
%      debug_scraps
%     'pause at end of season_mergeAnnualizeDat'  
% end

 end
 
end
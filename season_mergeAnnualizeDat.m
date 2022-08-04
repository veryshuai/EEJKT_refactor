function [mat_yr_sales, firm_yr_sales,iterX_in] = season_mergeAnnualizeDat(all_seas, som_seas, mm, mat_cols,iterX_in)
 
 % all_seas and som_seas: 
 %  (1) t, (2) season, (3) year, (4) initial state, (5) exporter id, (6) ending state,
 %  (7) match revenue,(8) #shipments,(9) exporter age (#periods), (10) match age w/in year

 t = iterX_in.t;
 
%  if t == 72
%      'pause here'
%      pickID = [ 0 0 0 0 1 0 0 0 0 0];
%      pickIDyr = kron(ones(1,12),pickID);
%      test2a = all_seas.*(ones(size(all_seas,1),1)*pickIDyr);
%      firmIDa = max(test2a,[],2);
%      test2s = som_seas.*(ones(size(all_seas,1),1)*pickIDyr);
%      firmIDs = max(test2s,[],2);
%  end
 
 a = find(sum(all_seas,2)>0); % non-zero rows of all_seas
 s = find(sum(som_seas,2)>0); % non-zero rows of som_seas


 %% create indicator for season in which firm flips
  firm_age = iterX_in.cumage(:,t-mm.pd_per_yr+1:t); % age in periods, current year  
 if t > mm.pd_per_yr
  firm_age_lag = iterX_in.cumage(:,t-mm.pd_per_yr:t-1); % age in periods, last year 
 else
  firm_age_lag = zeros(size(firm_age,1),1);
 end
 flip_seas = (firm_age-firm_age_lag)<0; 
 flip_firm = sum(flip_seas,2)>0; % equals 1 when age drops during year
 
 % diagnostics only:
 % find(sum(flip_seas,2))
 
 % iterX_in.flip_ndx gives the season (if any) in which a firm slot turns over
 % (It identifies the first period after a firm exit. It is 0 for non-flippers.
 iterX_in.flip_ndx = zeros(length(flip_seas),1);
 [rr,cc] = find(flip_seas);
 for jj=1:length(rr)
 iterX_in.flip_ndx(rr(jj)) = cc(jj);
 end
           
 %% aggregate match-specific sales and shipments across seasons
 pick_sales     = kron(ones(mm.pd_per_yr,1),[zeros(6,2);eye(2);zeros(2,2)]); % to pick off sales and shipments from all seasons
 mat_yr_sales_a = all_seas(a,:)*pick_sales; % add up within-yr. sales for matches that are active at year's end
 mat_yr_sales_s = som_seas(s,:)*pick_sales; % add up within-yr. sales for matches that die before year's end
 
 % create firm tags for each match
 fndx_a = zeros(size(a,1),mm.pd_per_yr);
 fndx_s = zeros(size(s,1),mm.pd_per_yr);  
 age_a  = zeros(size(a,1),mm.pd_per_yr);
 age_s  = zeros(size(s,1),mm.pd_per_yr);
 
 match_age_s = zeros(size(s,1),mm.pd_per_yr);
 
 age_old = max(firm_age.*(1-flip_seas),[],2);
 age_new = max(firm_age.*flip_seas,[],2);

 %% Find firm IDs, match age, and firm age for populated rows by season
 for ss=1:mm.pd_per_yr 
   fndx_a(:,ss) = all_seas(a,ss*mat_cols-5); % each season, firm ID is 6th col. from end, if active
   fndx_s(:,ss) = som_seas(s,ss*mat_cols-5);

   age_a(:,ss) = zeros(length(a),1);
   age_s(:,ss) = zeros(length(s),1);
   
     for jj=1:size(firm_age,1)
       firm_age2 = age_old(jj)*(1-flip_seas(jj,ss)) + ...
                   age_new(jj)*flip_seas(jj,ss);
       age_a(:,ss)  = age_a(:,ss) + (fndx_a(:,ss)==jj)*firm_age2;
       age_s(:,ss)  = age_s(:,ss) + (fndx_s(:,ss)==jj)*firm_age2;
     end       
   
   % w/in yr match age is last column, ea. season:  
   match_age_s(:,ss)  = som_seas(s,ss*mat_cols);  
%  match_age_a(:,ss)  = all_seas(a,ss*mat_cols);
   
 end

 firm_a = max(fndx_a,[],2); % get firm IDs for each row of all_seas
 firm_s = max(fndx_s,[],2); % get firm IDs for each row of som_seas
 age_aa = max(age_a,[],2);  % find max firm age during yr., all_seas
 age_ss = max(age_s,[],2);  % find max firm age during yr., som_seas
 
 match_age_s = max(match_age_s,[],2); % find within-period match age for matches that fail by eoy
%  match_age_a = max(match_age_a,[],2); % find within-period match age surviving at eoy
 
 temp0 = sortrows([[firm_a,age_aa];[firm_s,age_ss]],[1 2]); % stack all match obs. on firm, age and sort
     
 % assign max of firm age for each firm_ID to all matches for continuing firms
 
 %% find within-year entry month for each match 
 z_a = zeros(length(a),mm.pd_per_yr);
 z_s = zeros(length(s),mm.pd_per_yr);
 for ss = 1:mm.pd_per_yr
     z_a(:,ss) = (all_seas(a,10*ss-6)>0)*ss;
     z_s(:,ss) = (som_seas(s,10*ss-6)>0)*ss;
 end
 
 entry_seas_a = zeros(size(z_a,1),1); % will contain matches' season of entry 
 for ii=1:size(z_a,1)
     if sum(z_a(ii,:),2) >0
     entry_seas_a(ii) = min(z_a(ii,z_a(ii,:)>0),[],2) - 1;
     end
 end
  entry_seas_s = zeros(size(z_s,1),1); % will contain season of entry 
 for ii=1:size(z_s,1)
     if sum(z_s(ii,:),2) >0
     entry_seas_s(ii) = min(z_s(ii,z_s(ii,:)>0),[],2) - 1;
     end
 end

%% stack all matches in mat_yr_sales     
if size(temp0,1) > 0

 mat_yr_sales =...
     sortrows([[firm_a,mat_yr_sales_a,all_seas(a,4),all_seas(a,mat_cols*mm.pd_per_yr-4),all_seas(a,mat_cols*mm.pd_per_yr),age_a(:,mm.pd_per_yr)];...
               [firm_s,mat_yr_sales_s,som_seas(s,4),som_seas(s,mat_cols*mm.pd_per_yr-4),match_age_s, max(age_s,[],2)]],1);

entry_month = sortrows([[firm_a,entry_seas_a];[firm_s,entry_seas_s]],1);
             
% mat_yr_sales: [firm ID, w/in yr sales, w/in yr shpmts., boy Z, eoy Z, w/in yr match age, firm age]


%% Distinguish match of entering firms from matches of exiting firms  
%  (Add 0.5 to the firm IDs associated with the former matches.)
  F_active = unique(mat_yr_sales(:,1)); 

try
  if sum(flip_firm)>0
  for j=1:length(F_active)
     j_type = mat_yr_sales(:,1)==F_active(j); 
     % adjust firm ID to distinguish post-flip matches
     mat_yr_sales(j_type,1) =...
              mat_yr_sales(j_type,1) + 0.5.*flip_firm(F_active(j))...
           .*(entry_month(j_type,2)>=iterX_in.flip_ndx(F_active(j)));
     % adjust firm age   
     mat_yr_sales(j_type,7) = mat_yr_sales(j_type,7)...
             .*(entry_month(j_type,2)<iterX_in.flip_ndx(F_active(j)))...
             + (entry_month(j_type,2)>=iterX_in.flip_ndx(F_active(j)))...
              .*(mm.pd_per_yr - iterX_in.flip_ndx(F_active(j)));
  end
  end
catch
    'problem in season_mergeAnnualizeDat lines 118-122'
end
%% Sum over matches for each active firm

 % tenure  = mat_yr_sales(:,7)>mm.pd_per_yr; % matches with firm age > 1 year
 
  firm_yr_sales = double.empty(0,4);    
  F_active = unique(mat_yr_sales(:,1)); 
  for ff=1:length(F_active)
      % incumbents
      select = find(mat_yr_sales(:,1)==F_active(ff));
      if size(select,1)>0
      firm_yr_sales = [firm_yr_sales;...
        [F_active(ff),sum(mat_yr_sales(select,2:3),1),max(mat_yr_sales(select,7))] ];
      end
  end
   

else
  firm_yr_sales = double.empty(0,4);
  mat_yr_sales  = double.empty(0,7);
end

fprintf('\ryear and pt_ndx =%4.0f %4.0f\n', [iterX_in.year, iterX_in.pt_ndx] )

if t==540
    'pause at end of season_mergeAnnualizeDat'
end

end
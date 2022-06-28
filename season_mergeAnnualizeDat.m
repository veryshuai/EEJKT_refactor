function [mat_yr_sales, firm_yr_sales] = season_mergeAnnualizeDat(all_seas, som_seas, mm, mat_cols)
 a = find(sum(all_seas,2)>0); % non-zero rows of all_seas
 s = find(sum(som_seas,2)>0); % non-zero rows of som_seas

 % aggregate match-specific sales and shipments across seasons
 pick_sales     = kron(ones(mm.pd_per_yr,1),[zeros(6,2);eye(2);zeros(2,2)]); % to pick off sales and shipments from all seasons
 mat_yr_sales_a = all_seas(a,:)*pick_sales; % add up within-yr. sales for matches that are active at year's end
 mat_yr_sales_s = som_seas(s,:)*pick_sales; % add up within-yr. sales for matches that die before year's end
 
 % create firm tags for each match
 fndx_a = zeros(size(a,1),mm.pd_per_yr);
 fndx_s = zeros(size(s,1),mm.pd_per_yr);  
 age_a  = zeros(size(a,1),mm.pd_per_yr);
 age_s  = zeros(size(s,1),mm.pd_per_yr);
 match_age_s = zeros(size(s,1),mm.pd_per_yr);
 for ss=1:mm.pd_per_yr % find firm IDs, match age, and firm age for populated rows by season
   fndx_a(:,ss)       = all_seas(a,ss*mat_cols-5); % each season, firm ID is 6th col. from end, if active
   fndx_s(:,ss)       = som_seas(s,ss*mat_cols-5);
   age_a(:,ss)        = all_seas(a,ss*mat_cols-1); % each season, firm age is next to last column, if active
   age_s(:,ss)        = som_seas(s,ss*mat_cols-1);
   match_age_s(:,ss)  = som_seas(s,ss*mat_cols);   % each season, match age is last column
 end
 firm_a = max(fndx_a,[],2); % get firm IDs for each row of all_seas
 firm_s = max(fndx_s,[],2); % get firm IDs for each row of som_seas
 age_aa = max(age_a,[],2);  % find max age for each firm during yr., all_seas
 age_ss = max(age_s,[],2);  % find max age for each firm during yr., som_seas
 match_age_s = max(match_age_s,[],2); % find match age for matches that die in current year
 % takes care of odd case where matches generate no shipments, thus firm age = 0
age_ss = max(age_ss,match_age_s);
temp0 = sortrows([[firm_a,age_aa];[firm_s,age_ss]],[1 2]); % stack all match obs. on firm, age and sort
 
%  if size(temp0,1) > 20
%      'pause here'
%  end
     
 % assign max of firm age for each firm_ID to all matches for that firm
 if size(temp0,1) > 0
  ff_age = unique(temp0(:,1));
     age = temp0(:,2);
  for f_ID = 1:size(ff_age,1)
      mask = temp0(:,1) == ff_age(f_ID);
      mx_age = max(age(mask));
      age(mask) = mx_age;
  end
      
 % age = age./mm.pd_per_yr
 % added following line to get integer initial years:
 % age = floor(age./mm.pd_per_yr) + 1;
 
% stack match observations in all_seas and in som_seas, then sort by firm
% and tag firm age on the end (which is already sorted)
 mat_yr_sales = [sortrows([[firm_a,mat_yr_sales_a,all_seas(a,4),all_seas(a,mat_cols*mm.pd_per_yr-4),all_seas(a,mat_cols*mm.pd_per_yr)];...
                  [firm_s,mat_yr_sales_s,som_seas(s,4),som_seas(s,mat_cols*mm.pd_per_yr-4), match_age_s]],1),age];
% mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age in periods (w/in year), firm age in periods] 
 
% NOTE: sometimes, for firms in their first year (age<=12), match age can
% exceed firm age by one period. This has no effect on the results, since
% match age is only used once it is converted to years.


% finally sum over matches for each active firm
  temp    = dummyvar(mat_yr_sales(:,1));  % create firm_ID dummies for each match
  temp2   = sum(temp,1);                  % count match obs. on each firm
  temp3   = find(temp2>0);                % get rid of firm_IDs that don't occur in mat_yr_sales
  firmdum = temp(:,temp3);                % dummy matrix for firms with matches only
  firm_yr_sales = [mat_yr_sales(:,1)'*firmdum./temp2(:,temp3);...
      mat_yr_sales(:,2:3)'*firmdum; mat_yr_sales(:,7)'*(firmdum./sum(firmdum))]';
% firm_yr_sales: [firmID,sales,#shipments,firm age]
else
  firm_yr_sales = double.empty(0,4);
  mat_yr_sales  = double.empty(0,7);
end
end
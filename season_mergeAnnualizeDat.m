function [mat_yr_sales, firm_yr_sales,iterX_in] = season_mergeAnnualizeDat(all_seas, som_seas, mm, mat_cols,iterX_in)
 
 % all_seas and som_seas: 
 %  (1) t, (2) season, (3) year, (4) initial state, (5) exporter id, (6) ending state,
 %  (7) match revenue,(8) #shipments,(9) exporter age (#periods), (10) match age w/in year

 t = iterX_in.t;
 
 a = find(sum(all_seas,2)>0); % non-zero rows of all_seas
 s = find(sum(som_seas,2)>0); % non-zero rows of som_seas


 %% create indicator for season in which firm flips
  firm_age = iterX_in.cumage(:,t-mm.pd_per_yr+1:t); % age in periods, current year  
 if t > 12
     firm_age_lag = iterX_in.cumage(:,t-mm.pd_per_yr:t-1); % age in periods, current year 
 else
     firm_age_lag = zeros(size(firm_age,1),1);
 end
 flip_seas = (firm_age-firm_age_lag)<0; 
 flip_firm = sum(flip_seas,2)>0; % equals 1 when age drops during year
 
 % iterX_in.flip_ndx gives the season (if any) in which a firm slot turns over
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
     
   match_age_s(:,ss)  = som_seas(s,ss*mat_cols);   % w/in yr match age is last column, ea. season
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
     entry_seas_a = min(z_a,[],2);
     entry_seas_s = min(z_s,[],2);

%% stack all matches in mat_yr_sales     
if size(temp0,1) > 0

 mat_yr_sales =...
     sortrows([[firm_a,mat_yr_sales_a,all_seas(a,4),all_seas(a,mat_cols*mm.pd_per_yr-4),all_seas(a,mat_cols*mm.pd_per_yr),age_a(:,mm.pd_per_yr)];...
               [firm_s,mat_yr_sales_s,som_seas(s,4),som_seas(s,mat_cols*mm.pd_per_yr-4),match_age_s, max(age_s,[],2)]],1);

entry_month = sortrows([[firm_a,entry_seas_a];[firm_s,entry_seas_s]],1);
             
% mat_yr_sales: [firm ID, w/in yr sales, w/in yr shpmts., boy Z, eoy Z, w/in yr match age, firm age]


%% Distinguish match of entering firms from matches of exiting firms by 
%  adding 0.5 to the firm IDs associated with the former matches
try
    
if sum(flip_firm)>0
  for j=1:length(flip_firm)
%   mat_yr_sales(:,1) = mat_yr_sales(:,1)...
%       + 0.5*(mat_yr_sales(:,1)==nn(j)).*(mat_yr_sales(:,7)==age_new(nn(j)));

   mat_yr_sales(:,1) = mat_yr_sales(:,1)...
      + 0.5.*flip_firm(j).*(mat_yr_sales(:,1)==j).*(entry_month(:,2)>iterX_in.flip_ndx(j)) ...
      .*(mat_yr_sales(:,7)==age_new(j));
  end
end

catch
    'problem in season_mergeAnnualizeDat lines 96-106'
end

%% Sum over matches for each active firm

  tenure  = mat_yr_sales(:,7)>mm.pd_per_yr; % matches with firm age > 1 year
 
  firm_yr_sales_incb = double.empty(0,4); % for incumbent firms
  firm_yr_sales_new  = double.empty(0,4); % for new firms
  
 try
  F_active = unique(mat_yr_sales(:,1)); 
  for ff=1:length(F_active)
      % incumbents
      select = find(mat_yr_sales(:,1).* tenure==F_active(ff));
      if size(select,1)>0
      firm_yr_sales_incb = ...
        [firm_yr_sales_incb;...
        [F_active(ff),sum(mat_yr_sales(select,2:3),1),max(mat_yr_sales(select,7))] ];
      end
      % new entrants 
      select = find(mat_yr_sales(:,1).*(1-tenure)==F_active(ff));
      if size(select,1)>0
      firm_yr_sales_new = ...
        [firm_yr_sales_new;...
        [F_active(ff),sum(mat_yr_sales(select,2:3),1),max(mat_yr_sales(select,7))] ];
      end
  end
  firm_yr_sales = [ firm_yr_sales_incb; firm_yr_sales_new];
  
  catch
      'problem in season_mergeAnnualizeDat line 96-112'
  end    

else
  firm_yr_sales = double.empty(0,4);
  mat_yr_sales  = double.empty(0,7);
end

fprintf('\ryear and pt_ndx =%4.0f %4.0f\n', [iterX_in.year, iterX_in.pt_ndx] )

end
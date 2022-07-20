function [mat_yr_sales, firm_yr_sales,flip_firm] = season_mergeAnnualizeDat(all_seas, som_seas, mm, mat_cols,iterX_in)
 
 % all_seas and som_seas: 
 %  (1) t, (2) season, (3) year, (4) initial state, (5) exporter id, (6) ending state,
 %  (7) match revenue,(8) #shipments,(9) exporter age (#periods), (10) match age w/in year

 a = find(sum(all_seas,2)>0); % non-zero rows of all_seas
 s = find(sum(som_seas,2)>0); % non-zero rows of som_seas
 
 t        = iterX_in.t;
 firm_age = iterX_in.cumage(:,t-mm.pd_per_yr+1:t);     
 flip_seas = [zeros(size(firm_age,1),1),...
              firm_age(:,2:mm.pd_per_yr)-firm_age(:,1)<0]; 
 flip_firm = sum(flip_seas,2)>0; % equals 1 when age drops during year
 
 % Allows each firm slot to turnover up to once per year
 % (Multiple turnovers have same effect on age in years as a single turnover.)
          
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
 
 age_old = max(firm_age.*(1-flip_seas),[],2);
 age_new = max(firm_age.*flip_seas,[],2);

 for ss=1:mm.pd_per_yr % find firm IDs, match age, and firm age for populated rows by season
   fndx_a(:,ss) = all_seas(a,ss*mat_cols-5); % each season, firm ID is 6th col. from end, if active
   fndx_s(:,ss) = som_seas(s,ss*mat_cols-5);
  
%    age_a(:,ss)        = all_seas(a,ss*mat_cols-1); % each season, firm age is next to last column, if active
%    age_s(:,ss)        =  firm_age(:,ss);
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
if size(temp0,1) > 0

 mat_yr_sales =...
     sortrows([[firm_a,mat_yr_sales_a,all_seas(a,4),all_seas(a,mat_cols*mm.pd_per_yr-4),all_seas(a,mat_cols*mm.pd_per_yr),age_a(:,mm.pd_per_yr)];...
               [firm_s,mat_yr_sales_s,som_seas(s,4),som_seas(s,mat_cols*mm.pd_per_yr-4),match_age_s, max(age_s,[],2)]],1);

% mat_yr_sales: [firm ID, w/in yr sales, w/in yr shpmts., boy Z, eoy Z, w/in yr match age, firm age]

% to distinguish matches pertaining to flipped firm slots, add 0.5 to their firm IDs
 
F_active = unique(mat_yr_sales(:,1)); 
nn = F_active(flip_firm); 
for j=1:length(nn)
  mat_yr_sales(:,1) = mat_yr_sales(:,1)...
      + 0.5*(mat_yr_sales(:,1)==nn(j)).*(mat_yr_sales(:,7)==age_new(nn(j)));
end

% finally sum over matches for each active firm

  tenure  = mat_yr_sales(:,7)>mm.pd_per_yr; % matches with firm age > 1 year
 
  firm_yr_sales_incb = double.empty(0,4); % for incumbent firms
  firm_yr_sales_new  = double.empty(0,4); % for new firms
  
 try
  F_active = unique(mat_yr_sales(:,1)); 
  for ff=1:length(F_active)
      
      select = find(mat_yr_sales(:,1).* tenure==F_active(ff));
      if size(select,1)>0
      firm_yr_sales_incb = ...
        [firm_yr_sales_incb;...
        [F_active(ff),sum(mat_yr_sales(select,2:3),1),max(mat_yr_sales(select,7))] ];
      end
      
      select = find(mat_yr_sales(:,1).*(1-tenure)==F_active(ff));
      if size(select,1)>0
      firm_yr_sales_new = ...
        [firm_yr_sales_new;...
        [F_active(ff),sum(mat_yr_sales(select,2:3),1),max(mat_yr_sales(select,7))] ];
      end
  end
  firm_yr_sales = [ firm_yr_sales_incb; firm_yr_sales_new];
  
  catch
      'problem in season_mergeAnnualizeDat line 71-78'
  end    
      

else
  firm_yr_sales = double.empty(0,4);
  mat_yr_sales  = double.empty(0,7);
end
end
function [mat_yr_sales, firm_yr_sales] = season_mergeAnnualizeDat(all_seas, som_seas, mm, mat_cols,iterX_in)
 
 % all_seas and som_seas: 
 %  (1) t, (2) season, (3) year, (4) initial state, (5) exporter id, (6) ending state,
 %  (7) match revenue,(8) #shipments,(9) exporter age (#periods), (10) match age w/in year

 a = find(sum(all_seas,2)>0); % non-zero rows of all_seas
 s = find(sum(som_seas,2)>0); % non-zero rows of som_seas
 
 t        = iterX_in.t;
 firm_age = iterX_in.cumage(:,t-mm.pd_per_yr+1:t);
% exit_yr  = sum(iterX_in.exit_firm(:,t-mm.pd_per_yr+1:t),2);

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
   fndx_a(:,ss) = all_seas(a,ss*mat_cols-5); % each season, firm ID is 6th col. from end, if active
   fndx_s(:,ss) = som_seas(s,ss*mat_cols-5);
  
%    age_a(:,ss)        = all_seas(a,ss*mat_cols-1); % each season, firm age is next to last column, if active
%    age_s(:,ss)        =  firm_age(:,ss);
   age_a(:,ss) = zeros(length(a),1);
   age_s(:,ss) = zeros(length(s),1);
   for jj=1:size(firm_age,1)
     age_a(:,ss)  = age_a(:,ss) + (fndx_a(:,ss)==jj)*firm_age(jj,ss);
     age_s(:,ss)  = age_s(:,ss) + (fndx_s(:,ss)==jj)*firm_age(jj,ss);
   end 
   match_age_s(:,ss)  = som_seas(s,ss*mat_cols);   % w/in yr match age is last column, ea. season
   match_age_a(:,ss)  = all_seas(a,ss*mat_cols);
 end

 firm_a = max(fndx_a,[],2); % get firm IDs for each row of all_seas
 firm_s = max(fndx_s,[],2); % get firm IDs for each row of som_seas
 age_aa = max(age_a,[],2);  % find max age for each firm during yr., all_seas
 age_ss = max(age_s,[],2);  % find max age for each firm during yr., som_seas
 % eop_age_a = age_a(:,mm.pd_per_yr); % end of period firm age

 match_age_s = max(match_age_s,[],2); % find within-period match age for matches 
 match_age_a = max(match_age_a,[],2);
 
 % take care of odd case where matches generate no shipments, thus firm age = 0:
%   age_aa = max(age_aa,match_age_a);
%   age_ss = max(age_ss,match_age_s);

 temp0 = sortrows([[firm_a,age_aa];[firm_s,age_ss]],[1 2]); % stack all match obs. on firm, age and sort
     
 % assign max of firm age for each firm_ID to all matches for that firm
if size(temp0,1) > 0

 mat_yr_sales =...
     sortrows([[firm_a,mat_yr_sales_a,all_seas(a,4),all_seas(a,mat_cols*mm.pd_per_yr-4),all_seas(a,mat_cols*mm.pd_per_yr),age_a(:,mm.pd_per_yr)];...
               [firm_s,mat_yr_sales_s,som_seas(s,4),som_seas(s,mat_cols*mm.pd_per_yr-4),match_age_s, max(age_s,[],2)]],1);

% finally sum over matches for each active firm

% temp    = dummyvar(mat_yr_sales(:,1));  % create firm_ID dummies for each match
% temp    = dummyvar(mm.periods*mat_yr_sales(:,1)+ mat_yr_sales(:,7));  
% temp    = dummyvar(mat_yr_sales(:,1));  
  tenure  = mat_yr_sales(:,7)>mm.pd_per_yr; % matches with firm age > 1 year
%   newagg  = tenure(:,1);  % new firm dummy
%   oldagg  = tenure(:,2);  % incumbent firm dummy
  
  F_active = unique(mat_yr_sales(:,1)); % 
  firm_yr_sales_inc = double.empty(0,4); % for incumbent firms
  firm_yr_sales_new = double.empty(0,4); % for new firms
  
 try
  
  for ff=1:length(F_active)
      
      select = find(mat_yr_sales(:,1).* tenure==F_active(ff));
      if size(select,1)>0
      firm_yr_sales_inc = ...
        [firm_yr_sales_inc;...
        [F_active(ff),sum(mat_yr_sales(select,2:3),1),max(mat_yr_sales(select,7))] ];
      end
      
      select = find(mat_yr_sales(:,1).*(1-tenure)==F_active(ff));
      if size(select,1)>0
      firm_yr_sales_new = ...
        [firm_yr_sales_new; [F_active(ff),sum(mat_yr_sales(select,2:3),1),max(mat_yr_sales(select,7))] ];
      end
  end
  firm_yr_sales = [ firm_yr_sales_inc; firm_yr_sales_new];
  
  catch
      'problem in season_mergeAnnualizeDat line 71-78'
  end    
      
%   temp2  = zeros(2,length(Nfirms));
%   temp3  = zeros(2,length(Nfirms));
  

%   for jj = 1:2
%   temp2(jj,:) = sum(temp.*tenure(:,jj),1); % jj=1 => new firm; jj=2 => incumbent
%   temp3(jj,:) = temp2(jj,:)>0;
%    firmdum = temp.*(tenure(:,jj).*temp3(jj,:));
%   end

  
  % create firm_ID dummies for each match, recognizing possible entry/exit
  % and age reset                                                                     
%   temp2   = sum(temp,1);                  % count match obs. on each firm
%   temp3   = find(temp2>0);                % get rid of firm_IDs that don't occur in mat_yr_sales
%   firmdum = temp(:,temp3);                % dummy matrix for firms with matches only
%   
%   firm_yr_sales = [mat_yr_sales(:,1)'*firmdum./temp2(:,temp3);...
%       mat_yr_sales(:,2:3)'*firmdum; mat_yr_sales(:,7)'*(firmdum./sum(firmdum))]';
% firm_yr_sales: [firmID,sales,#shipments,firm age]
else
  firm_yr_sales = double.empty(0,4);
  mat_yr_sales  = double.empty(0,7);
end
end
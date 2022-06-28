function [mat_yr_sales,firm_yr_sales] = season_merge(iterX_in,mm)

%  This function takes a year's worth of season-specific match
%  outcomes for a particular type of firm and organizes the information
%  into annual aggregates at the match and firm level. These are passed
%  back to discrete_sim.m, where they are used to construct moments.

%  iterX_in.seas_tran: [t,season,year,initial state,firm_ID,new state,rev,shipments,firm age (in periods)];

%% build within-yr trajectories for all matches--continuing, new, and dying

  N_firms = mm.sim_firm_num_by_prod_succ_type(iterX_in.pt_ndx);

  mat_cols = size(iterX_in.seas_tran{1},2)+1;            % the +1 makes room for a match age variable
  all_seas = zeros(iterX_in.N_match,mat_cols*mm.pd_per_yr); % to hold data on matches present at end of year
  som_seas = zeros(iterX_in.N_match,mat_cols*mm.pd_per_yr); % to hold data on matches that die before end of year
  mat_age  = ones(size(iterX_in.seas_tran{1},1),1);      % initial match age within year (season 1)
  
  smat_tran  = [sortrows(iterX_in.seas_tran{1},[5 6]),mat_age]; 
  % iterX_in.seas_tran: [t, season, year, mat_tran, #shipments, exporter age(in periods)]; where
  % mat_tran:  [initial state, exporter id, ending state, match revenue]
  
  % smat_tran: 
  %  (1) t, (2) season, (3) year, (4) initial state, (5) exporter id, (6) ending state,
  %  (7) match revenue,(8) #shipments,(9) exporter age (#periods), (10) match age w/in year
  
  nrt        = size(smat_tran,1);
  ff_die     = find(smat_tran(:,6)==0); % death of existing match--endog. & exog. 
  ff_cont    = find(smat_tran(:,6)>0);  % find matches that continue next period, new and existing
  
  smat_die   = smat_tran(ff_die,:);    
  smat_cont  = smat_tran(ff_cont,:);
  
  all_cntr = size(ff_cont,1);   % counts the rows filled in all_seas
  som_cntr = size(ff_die,1);    % counts the rows filled in som_seas

  som_seas(1:som_cntr,1:mat_cols) = smat_die; 
  all_seas(1:all_cntr,1:mat_cols) = smat_cont;
 
  nrt_lag = nrt;

% count the number of b.o.y. matches for each firm (firms are columns):
  match_count_lag = sum(smat_tran(:,5).*ones(1,N_firms) - (1:1:N_firms)==0,1);
           
for ss=2:mm.pd_per_yr
   smat_tran = iterX_in.seas_tran{ss};
   nrt = size(smat_tran,1);
   lcb = mat_cols*(ss-1)+1;  % lower column bound for horizontal additions to all_seas and some_seas
   ucb = mat_cols*ss;        % upper column bound for horizontal additions to all_seas and some_seas       
%%     
   match_count = zeros(1,N_firms);  
   if nrt == 0  % no matches in the current season--move all matches to the 
                % som_seas matrix and empty the all_seas matrix

      ff_move =  find(all_seas(:,lcb-4)+all_seas(:,lcb-5)+all_seas(:,lcb-6)>0); % rows populated last period
      if size(ff_move,1)>0
      som_seas(som_cntr+1:som_cntr+size(ff_move,1),:) = all_seas(ff_move,:); 
      som_cntr = size(som_seas,1);      
      all_seas = zeros(iterX_in.N_match,mat_cols*mm.pd_per_yr); 
      all_cntr = 0;
      end
            
   else  % nrt>0 positive number of matches in the current season 
       
%      Firms that go to zero active matches (exit the market) have no records
%      carried to the next period, even if they have positive lagged z values.
%      Need to move these matches to som_seas matrix.
       
      if nrt_lag > 0 % there were some matches in previous season
       match_count = sum(smat_tran(:,5).*ones(1,N_firms) - (1:1:N_firms)==0,1); % current match count, by firm
       ff_firm_exit = find((match_count==0).*(match_count_lag>0));  % firm exits between previous and current season
       N_firm_exit = sum(ff_firm_exit>0); 

% NOTE: This code treats a firm as exiting if it goes to zero matches at any 
%       point during the year, even if it makes subsequent matches before 
%       year's end. This is not the way exit is defined in the data, where 
%       firms must go 12 months without a shipment to be flagged as exiting.
       
           if N_firm_exit>0    
             for jj = 1:N_firm_exit
                 % identify previous season matches of exiting firms, allowing for multiple matches:
                 ff_matfirm_exit = find(ones(size(all_seas,1),1).*ff_firm_exit(1,jj)==all_seas(:,lcb-6)); 
                 % last period matches with exiting firm_ID. Move these 
                 % matches from all_seas to som_seas and increment som_cntr:
                 som_seas(som_cntr+1:som_cntr+size(ff_matfirm_exit,1),:) = all_seas(ff_matfirm_exit,:); 
                 all_seas(ff_matfirm_exit,:) = zeros(size(ff_matfirm_exit,1),mat_cols*mm.pd_per_yr);
                 som_cntr = som_cntr+size(ff_matfirm_exit,1);
             end
           end
      end
%%         
    ff_new     = find(smat_tran(:,4)==0); % current season new matches   
    ff_incum   = find(smat_tran(:,4)>0);  % current season inherited matches 
    
    ff_all_active  = find(all_seas(:,lcb-5)>0); % previous season matches that should have passed to current season (eop Z>0)
    all_cntr       = size(ff_all_active,1);
     
    temp1 = sortrows(smat_tran(ff_incum,:),[5 4]); % sort current matches, incumbents, by exporter ID, bop Z
    temp2 = sortrows(all_seas(ff_all_active,:),[lcb-6 lcb-5]); % sort continuing matches from previous period by exporter ID, eop Z
    
    try
      assert(size(ff_incum,1)==size(ff_all_active,1)); % do counts match?
      assert(sum((temp1(:,5)-temp2(:,lcb-6)).^2)==0);  % do firms match?
      assert(sum((temp1(:,4)-temp2(:,lcb-5)).^2)==0);  % do Z's match?
    catch
      warning('problem with splicing of firms across seasons')
      temp1(:,4:5) 
      temp2(:,lcb-6:lcb-5)  
    end
    
%    load current season continuing matches into all_seas, incrementing match age by 1
     all_seas(1:all_cntr,1:ucb) = [temp2(:,1:lcb-1),temp1,temp2(:,lcb-1)+ones(size(temp1,1),1)] ;
%    add new rows to all_seas for new matches. Match age is one for all new matches.
     all_seas(all_cntr+1:all_cntr+size(ff_new,1),lcb:ucb) = [smat_tran(ff_new,:),ones(size(ff_new,1),1)];
     all_cntr = size(ff_all_active,1)+size(ff_new,1);      
%    move history of matches in their last season to som_seas    
     ff_die = find(all_seas(1:all_cntr,lcb+5)==0); % eop Z ==0
     som_seas(som_cntr+1:som_cntr+size(ff_die,1),1:ucb) = all_seas(ff_die,1:ucb);
     som_cntr = som_cntr + size(ff_die,1);
%    clear history of last-season matches out of all_seas and put surviving matches in first rows   
     ff_live = find(all_seas(:,lcb+5)>0); % continuing matches (eop Z>0)
     all_cntr = size(ff_live,1);    
     all_seas(1:all_cntr,:) = all_seas(ff_live,:); % moving survivors to first rows    
     empty_mat = zeros(iterX_in.N_match-all_cntr,mat_cols*mm.pd_per_yr);
     all_seas(all_cntr+1:iterX_in.N_match,:) = empty_mat;   % clear remaining rows

   end % end nrt >0 if block 
    match_count_lag = match_count; 
    nrt_lag = nrt; % JT: line wasn't in previous version of code but didn't matter 
end
%%  Package up the match info. for use in regressions

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


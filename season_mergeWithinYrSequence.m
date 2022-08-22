function [mat_cols, all_seas, som_seas] = season_mergeWithinYrSequence(mm, iterX_in)

  N_firms = mm.sim_firm_num_by_prod_succ_type(iterX_in.pt_ndx);

  mat_cols = size(iterX_in.seas_tran{1},2)+1;  % Holds all matches. The +1 makes room for a match age variable
  all_seas = zeros(iterX_in.N_match,mat_cols*mm.pd_per_yr); % To hold data on matches present at end of year
  som_seas = zeros(iterX_in.N_match,mat_cols*mm.pd_per_yr); % To hold data on matches that die before end of year
  mat_age  = ones(size(iterX_in.seas_tran{1},1),1); % Initial match age within year (season 1)
  
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

% NOTE: This function treats a firm as exiting if it goes to zero matches at any 
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
      fprintf('\r\n firm type = %4.0', iterX.in_pt_ndx); 
    
      temp1(:,4:5) 
      temp2(:,lcb-6:lcb-5)  
    end
    try
%    load current season continuing matches into all_seas, incrementing match age by 1
     all_seas(1:all_cntr,1:ucb) = [temp2(:,1:lcb-1),temp1,temp2(:,lcb-1)+ones(size(temp1,1),1)] ;
%    add new rows to all_seas for new matches. Match age is zero for all new matches.
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
    catch
        'problem in mergeWithinYrSequence lines 102-116 of seasonMergeWithinYrSequence'
    end
   end % end nrt >0 if block 
    match_count_lag = match_count; 
    nrt_lag = nrt; 
end
end
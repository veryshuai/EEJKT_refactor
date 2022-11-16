function [mat_cols, all_seas, som_seas, Zcut_eoy] = season_mergeWithinYrSequence(mm, iterX_in)

  N_firms = mm.sim_firm_num_by_prod_succ_type(iterX_in.pt_ndx);

    if isempty(iterX_in.seas_tran{1}) == 1
      smat_tran  = double.empty(0,10);  
      mat_cols   = 10;
    else
      mat_age  = ones(size(iterX_in.seas_tran{1},1),1); % Initial match age within year (season 1)
      smat_tran = [sortrows(iterX_in.seas_tran{1},[5 6]),mat_age]; 
      mat_cols = size(iterX_in.seas_tran{1},2)+1;  % Holds all matches. The +1 makes room for a match age variable
    end
    
  all_seas = zeros(iterX_in.N_match,mat_cols*mm.pd_per_yr); % To hold data on matches present at end of year
  som_seas = zeros(iterX_in.N_match,mat_cols*mm.pd_per_yr); % To hold data on matches that die before end of year
 
  % iterX_in.seas_tran: [t, season, year, mat_tran, #shipments, exporter age(in periods)]; where
  % mat_tran:  [initial state, exporter id, ending state, match revenue]

  % smat_tran: 
  %  (1) t, (2) season, (3) year, (4) initial state, (5) exporter id, (6) ending state,
  %  (7) match revenue,(8) #shipments,(9) exporter age (#periods), (10) match age w/in year

try
  Zcut_eoy_lag = iterX_in.Zcut_eoy_lag;
catch
 fprintf('\rIn SimulateHomeMatchesInnerSimUpdZHotel, t =%4.0f, pt_ndx = %4.0f\n', [ iterX_in.t,  iterX_in.pt_ndx] )
end
    
  nrt       = size(smat_tran,1);
  Zcut      = iterX_in.seas_Zcut(1);
  ff_keep   = find(max((smat_tran(:,4)>Zcut_eoy_lag),(smat_tran(:,4)==0))); % matches that didn't die last period
  smat_tran = smat_tran(ff_keep,:);
  ff_die    = find(smat_tran(:,6)<=Zcut); % death of existing match--endog. & exog. 
  ff_cont   = find(smat_tran(:,6)>Zcut);  % find matches that continue next period, new and existing
  
  smat_die  = smat_tran(ff_die,:);    
  smat_cont = smat_tran(ff_cont,:);
  
  all_cntr  = size(ff_cont,1);   % counts the rows filled in all_seas
  som_cntr  = size(ff_die,1);    % counts the rows filled in som_seas

  som_seas(1:som_cntr,1:mat_cols) = smat_die; 
  all_seas(1:all_cntr,1:mat_cols) = smat_cont;
 
  nrt_lag = nrt;

% count the number of b.o.y. matches for each firm (firms are columns):
  match_count_lag = sum(smat_tran(:,5).*ones(1,N_firms) - (1:1:N_firms)==0,1);
           
for ss=2:mm.pd_per_yr
   if isempty(iterX_in.seas_tran{ss}) == 0
     smat_tran = iterX_in.seas_tran{ss};
   else
     smat_tran  = double.empty(0,9);   
   end
   
   Zcut      = iterX_in.seas_Zcut(ss);
   Zcut_lag  = iterX_in.seas_Zcut(ss-1);
   nrt = size(smat_tran,1);
   lcb = mat_cols*(ss-1)+1;  % lower column bound for horizontal additions to all_seas and some_seas
   ucb = mat_cols*ss;        % upper column bound for horizontal additions to all_seas and some_seas       
%%     
   match_count = zeros(1,N_firms);  

   if isempty(iterX_in.seas_tran{ss})==1
%  if nrt == 0  % no matches in the current season--move all matches to the 
                % som_seas matrix and empty the all_seas matrix
   try
      ff_move =  find(all_seas(:,lcb-4)+all_seas(:,lcb-5)+all_seas(:,lcb-6)>0); % rows populated last period 
      if size(ff_move,1)>0
      som_seas(som_cntr+1:som_cntr+size(ff_move,1),:) = all_seas(ff_move,:); 
      som_cntr = size(som_seas,1);      
      all_seas = zeros(iterX_in.N_match,mat_cols*mm.pd_per_yr); 
      all_cntr = 0;
      end
   catch
     'problem in season_mergeWithinYrSequence lines 62-73' 
   end    
        
   else  % nrt>0 positive number of matches in the current season 
       
%      Firms that go to zero active matches (exit the market) have no records
%      carried to the next period, even if they have positive lagged z values.
%      Need to move these matches to som_seas matrix.
       
      if nrt_lag > 0 % there were some matches in previous season
%      ff_cont   = find(smat_tran(:,6)>Zcut);
%      smat_tran = smat_tran(ff_keep,:);
       match_count = sum((smat_tran(:,5).*ones(1,N_firms)) - (1:1:N_firms)==0,1); % current match count, by firm
       ff_firm_exit = find((match_count==0).*(match_count_lag>0));  % firm exits between previous and current season
       N_firm_exit  = sum(ff_firm_exit>0); 

% NOTE: This block treats a firm as exiting if it goes to zero matches at any 
%       point during the year, even if it makes subsequent matches before 
%       year's end. This is not the way exit is defined in the data, where 
%       firms must go 12 months without a shipment to be flagged as exiting.
       
           if N_firm_exit>0    
             for jj = 1:N_firm_exit
                 % identify previous season matches of exiting firms, allowing for multiple matches:
                 ff_matfirm_exit = find(ones(size(all_seas,1),1).*ff_firm_exit(1,jj)==all_seas(:,lcb-6)); 
                 % last period matches with an exiting firm #. Move these 
                 % matches from all_seas to som_seas and increment som_cntr:
                 som_seas(som_cntr+1:som_cntr+size(ff_matfirm_exit,1),:) = all_seas(ff_matfirm_exit,:); 
                 all_seas(ff_matfirm_exit,:) = zeros(size(ff_matfirm_exit,1),mat_cols*mm.pd_per_yr);
                 som_cntr = som_cntr+size(ff_matfirm_exit,1);
             end
           end
      end
%%         
    ff_new     = find(smat_tran(:,4)==0); % current season new matches   
    ff_incum   = find(smat_tran(:,4)>Zcut_lag);  % current season inherited matches

    ff_all_active  = find(all_seas(:,lcb-5)>Zcut_lag); % previous season matches that should have passed to current season (eop Z>Zcut)
    all_cntr       = size(ff_all_active,1);
     
    temp1 = sortrows(smat_tran(ff_incum,:),[5 4]); % sort current matches, incumbents, by exporter ID, bop Z
    temp2 = sortrows(all_seas(ff_all_active,:),[lcb-6 lcb-5]); % sort continuing matches from previous period by exporter ID, eop Z
    
    try
      assert(size(ff_incum,1)==size(ff_all_active,1)); % do counts match?
      assert(sum((temp1(:,5)-temp2(:,lcb-6)).^2)==0);  % do firms match?
      assert(sum((temp1(:,4)-temp2(:,lcb-5)).^2)==0);  % do Z's match?
    catch
       rows1 = size(ff_incum,1);
       rows2 = size(ff_all_active,1);
      fprintf('\r\n problem with splicing of matches across seasons: season_merge temp1 & temp2 \n')
      fprintf('\r\n rows of temp1 = %.2f, rows of temp2 = %.2f\n', [rows1 rows2]) 
      fprintf('\r\n period = %.2f, firm type = %.2f, market = %.2f\n', [iterX_in.t iterX_in.pt_ndx iterX_in.mkt])

        fileID4 = fopen('results/EEJKT_error_log.txt','a');
        fprintf(fileID4,'\r\n  ');
        fprintf(fileID4,'\r\n problem splicing matches across seasons in seasonMergeWithinYrSequence');
        fprintf(fileID4,'\r\n period = %.2f, firm type = %.2f, market = %.2f', [iterX_in.t iterX_in.pt_ndx iterX_in.mkt]);
        fprintf(fileID4,'\r\n params = ');
        fprintf(fileID4,'\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',mm.param_vec(1:6));
        fprintf(fileID4,'\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',mm.param_vec(7:12));
        fprintf(fileID4, '\r\n  ');   
        fclose(fileID4);
        firm_type = iterX_in.pt_ndx;
        problem_month = iterX_in.t;
        problem_market = iterX_in.mkt;
        params = mm.param_vec;
        save 'mismat_recs.mat' 'temp1' 'temp2' 'Zcut' 'smat_tran' 'params' 'firm_type' 'problem_month' 'problem_market','-append';  

    end

     try  
        
%    load current season continuing matches into all_seas, incrementing match age by 1
     all_seas(1:all_cntr,1:ucb) = [temp2(:,1:lcb-1),temp1,temp2(:,lcb-1)+ones(size(temp1,1),1)] ;
%    add new rows to all_seas for new matches. Match age is zero for all new matches.
     all_seas(all_cntr+1:all_cntr+size(ff_new,1),lcb:ucb) = [smat_tran(ff_new,:),ones(size(ff_new,1),1)];
     all_cntr = size(ff_all_active,1)+size(ff_new,1);  
      
%    move history of matches in their last season to som_seas    
     ff_die = find(all_seas(1:all_cntr,lcb+5)<=Zcut); % eop Z <= Zcut     
%    ff_die = find(all_seas(1:all_cntr,lcb+5)==0); % eop Z == 0

     som_seas(som_cntr+1:som_cntr+size(ff_die,1),1:ucb) = all_seas(ff_die,1:ucb);
     som_cntr = som_cntr + size(ff_die,1);
     
%    clear history of last-season matches out of all_seas and put surviving matches in first rows   
     ff_live = find(all_seas(:,lcb+5)>Zcut); % continuing matches (eop Z>Zcut)
%    ff_live = find(all_seas(:,lcb+5)>0); % continuing matches (eop Z>0)
     
     all_cntr = size(ff_live,1);    
     all_seas(1:all_cntr,:) = all_seas(ff_live,:); % moving survivors to first rows    
     empty_mat = zeros(iterX_in.N_match-all_cntr,mat_cols*mm.pd_per_yr);
     all_seas(all_cntr+1:iterX_in.N_match,:) = empty_mat;   % clear remaining rows
%    move matches that die at end of year out of all_seas and into som_seas     
     if ss==mm.pd_per_yr
         ff_eoy_die =  find((all_seas(1:all_cntr,lcb+5)<=Zcut).*(all_seas(1:all_cntr,lcb+3)>0));
         ff_eoy_live = find(all_seas(1:all_cntr,lcb+5) > Zcut);
%        ff_eoy_die =  find((all_seas(1:all_cntr,lcb+5)==0).*(all_seas(1:all_cntr,lcb+3)>0));
%        ff_eoy_live = find(all_seas(1:all_cntr,lcb+5) > 0);          
         
         som_cntr = som_cntr + size(ff_eoy_die,1);
         all_cntr = size(ff_eoy_live,1);
         som_seas(som_cntr+1:som_cntr+size(ff_eoy_die,1),1:ucb) = all_seas(ff_eoy_die,1:ucb);
         all_seas(1:all_cntr,:) = all_seas(ff_eoy_live,:); 
         all_seas(all_cntr+1:iterX_in.N_match,:) = zeros(iterX_in.N_match-all_cntr,mat_cols*mm.pd_per_yr);  %
     end
    catch
      fprintf('\r\n period = %.2f, firm type = %.2f, market =%.2f\n',[iterX_in.t iterX_in.pt_ndx iterX_in.mkt]) 
      fprintf('\r\n problem with splicing in seasonMergeWithinYrSequence: lines 151-183\n')
      if size(temp1,1) >20
      fileID3 = fopen('results/EEJKT_error_log.txt','a');
      fprintf(fileID3,'\r\n  ');
      fprintf(fileID3,'\r\n problem stacking matches in seasonMergeWithinYrSequence');
      fprintf(fileID3,'\r\n period = %.2f, firm type = %.2f, market = %.2f', [iterX_in.t iterX_in.pt_ndx iterX_in.mkt]);
      fprintf(fileID3,'\r\n params = ');
      fprintf(fileID3,'\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',mm.param_vec(1:6));
      fprintf(fileID3,'\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',mm.param_vec(7:12));
      fclose(fileID3);
      firm_type = iterX_in.pt_ndx;
      problem_month = iterX_in.t;
      problem_market = iterX_in.mkt;
      Xvec = mm.param_vec;
      save 'mismat2_recs.mat' 'Xvec' 'firm_type' 'problem_month' 'problem_market';    
      end
    end
   end % end nrt >0 if block 
    match_count_lag = match_count; 
    nrt_lag = nrt; 
    Zcut_eoy = Zcut;
end

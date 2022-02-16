% This function simulates panel data on domestic sales for all N_firms firms of a
% particular type. Type is determined by productivity and foreign theta.

function [agg_firm_yr_sales,agg_x_fsales_h,agg_y_fsales_h,agg_theta_h_firm,...
          agg_fmoms_h_xx,agg_fmoms_h_xy,agg_fysum_h,agg_fnobs_h,my_flag]...
    = matchdat_gen_h(N_firms,N_Z,N_Phi,N_mic_state,pd_per_yr,frac_of_year,periods,...
           max_ships,typemat,macro_state_h,theta_ndx,prod_ndx,mm,pmat_cum_h,...
           c_val_f_orig,cum_pz,poisCDF,prob_mdeath,prob_fdeath,...
           haz_ship,pmat_cum_z,N_theta1,th1_cdf)  

    %% create home theta draws  
    
    theta_df = (N_theta1+1)*ones(N_firms,1) - ...
    sum(ones(N_firms,1)*th1_cdf > rand(N_firms,1)*ones(1,N_theta1),2);
%   try
%     assert(max(theta_df)<=N_theta1 & min(theta_df>=1))
%   catch 
%      'theta_df out of bounds in discrete_sim_parfor3'
%      [min(theta_df),max(theta_df)]
%   end
  
  % list the non-zero values and their frequencies    
   [uv,~,idx]   = unique(theta_df);
   theta1_cntr  = [uv,accumarray(idx(:),1)];
   theta_h      = sortrows(theta_df);
  try
    assert(max(theta1_cntr(:,1))<=N_theta1 & min(theta1_cntr(:,1)>=1))
  catch 
     'theta1_cntr out of bounds in discrete_sim_parfor4'
     [min(theta_df),max(theta_df)]
  end

       
       %% Initialize matrices

  seas_tran = cell(1,pd_per_yr); % cells will hold one year's worth of season- and match-specific outcomes for all firms w/in type
  seas_Zcut = zeros(1,pd_per_yr);    % elements will hold season-specifics Z cut-offs for endog. drops
% Each firm begins with zero trials zero successes, macro state at median position
  lag_row = ones(N_firms,1);

  cur_cli_cnt  = zeros(N_firms,periods,1); % clients active in the current period
  add_cli_cnt  = zeros(N_firms,periods,1); % gross additions to client count
%  cum_meets    = zeros(N_firms,periods,1); % cumulative number of meetings
  cum_succ     = zeros(N_firms,periods,1); % cumulative number of successes
  new_firm     = zeros(N_firms,periods);   % will mark first period of a new firm
%  tot_ships    = zeros(N_firms,periods,1); % total number of shipments
  exog_deaths  = zeros(N_firms,periods,1); % number of exogenous match deaths
  micro_state  = ones(N_firms,periods,1); % scalar indices for #success/#meetings 
%  cur_cli_zst  = zeros(N_firms,N_Z);  % breaks down current client counts by z state
  lag_cli_zst  = zeros(N_firms,N_Z);  % breaks down lagged clients by z state
  new_cli_zst  = zeros(N_firms,N_Z);  % breaks down new client counts by z state
  die_cli_zst  = zeros(N_firms,N_Z);  % breaks down client death counts by z state
  surviv_zst   = zeros(N_firms,N_Z);  % breaks down surviving client counts by z state
  trans_zst    = zeros(N_firms,N_Z);  % counts survival types by firm after z innovations
%  drop_cnt     = zeros(N_firms,1);    % number of clients endogenously dropped
  flrlag       = ones(N_firms,1);     % initializing vector for age debugging
  cumage       = zeros(N_firms,1);    % initializing vector for age debugging 
%  cont_expr    = ones(N_firms,1);     % continuing exporter indicator
  
% create first observation on firm-year level aggregates (will concatenate below)
% max_match         = 50; % upper bound on number of matches to be counted 
% agg_match_count   = zeros(max_match,1);
agg_mat_yr_sales  = double.empty(0,9);
agg_firm_yr_sales = double.empty(0,6);
agg_theta_h_firm  = double.empty(0,1);
% agg_time_gaps     = double.empty(0,7);
% mkt_exit          = zeros(1,3);

agg_x_fsales_h    = zeros(0,2); % will contain rows of x matrix for regression
agg_y_fsales_h    = zeros(0,1); % will contain rows of y matrix for regression

tic
firm_cntr = 0;  % simulated firm counter, all types combined
chk1 = 0;
chk2 = 0;

% match level moment aggregators
agg_moms_xx   = zeros(4,4);
agg_moms_xy   = zeros(4,1);
agg_ysum      = 0;
agg_nobs      = 0;
agg_ship_obs  = 0;
agg_ln_ships = 0;

% firm level moment aggregators
agg_fmoms_h_xx = zeros(2,2);
agg_fmoms_h_xy = zeros(2,1);
agg_fysum_h    = 0;
agg_fnobs_h    = 0; 

agg_exit_moms_xx = zeros(6,6);
agg_exit_moms_xy = zeros(6,1);
agg_sum_succ_rate = 0;
agg_exit_obs = 0;
agg_sum_exits = 0;

agg_mat_exit_moms_xx  = zeros(5,5);
agg_mat_exit_moms_xy  = zeros(5,1);
agg_mat_obs       = 0;
agg_nmat_exit     = 0;

  trans_count    = zeros(N_Z+1,N_Z+1,N_firms); % counts transitions across buyer types, 
% for each seller type. New buyer types are considered type 0 at beginning of period, hence the +1. 
% Exiting firms are considered to move to type 0 at the the end of the period.
% columns: (1) initial z-state (2) new z-state (3) firm index, given type (4) firm type   

%  trans_tot    = zeros(N_Z+1,N_Z+1,periods,N_types); % transition counts through time, aggregeted across firms

%% Initialize stuff

% keep_cli will be used to select matches that are endogenously dropped.
  keep_cli      = ones(1,N_Z); % applies to clients existing in period 1
  keep_cli(1:5) =  zeros(1,5); % implying worst 5 client types from period 1 are dropped.

% load relevant transition probabilites for current micro type & macro
% state. pmat_cum_t holds cumulative transition probs across #success/#meeting pairs.
% It's loaded from pmat_cum, which is constructed in inten_sim using the relevant Q matrix.       
  year = 1;

%% TIME LOOP BEGINS HERE
  tlag = 1;
  season = 1;
  N_match = 0;
  N_match_lag = 0;
  firm_yr_sales_lag = zeros(N_firms,4);   
  % firm_yr_sales_lag will contain: [firmID,sales,#shipments,firm age]

  tic
  my_flag = 0;
for t = 2:1:periods

% update year, type and pmat_cum_t
    year = floor((t-1)/pd_per_yr);
    
% reset season when previous period completes a year
     if abs(floor((t-1)/pd_per_yr) - (t-1)/pd_per_yr) <= 0.001 
        season = 1;  
     end

%% gross additions to clients, before drops and deaths between t-1 and t

   stay     = ones(N_firms,1);  % will flag firms that continue, t-1 to t
   
%  Identify firms that had at least one shipment last year. For others, cont_expr=0.    
   if sum(firm_yr_sales_lag(:,3))>0 % at least one firm had positive shipments
    temp      = ones(N_firms,1).*(1:1:N_firms);
    temp2     = temp(:,firm_yr_sales_lag(firm_yr_sales_lag(:,3)>0,1)); 
    cont_expr = sum(((1:N_firms)' - temp2 == 0),2);
   else
    cont_expr = zeros(N_firms,1);
   end

   if  N_firms > 0
       trans_rands = zeros(N_firms,N_mic_state);
       cntr = 0;
       
  try
    assert(max(theta1_cntr(:,1))<=7 & min(theta1_cntr(:,1)>=1))
  catch 
     'theta1_cntr out of bounds in matchdat_gen_h1'
     [min(theta_df),max(theta_df)]
  end
       
  for ii = 1:size(theta1_cntr,1) 
           cntr2 = cntr+theta1_cntr(ii,2);
           type = find(typemat(:,2) == macro_state_h(t) & typemat(:,3) == theta1_cntr(ii,1) & typemat(:,4) == prod_ndx,1,'first');
           try  % try and catch can be eliminated once this is clearly working properly
           pmat_cum_ht = pmat_cum_h{type};
           catch
                [macro_state_h(t),theta1_cntr(ii,1),prod_ndx]
%                pmat_cum_h{type}
%                type
           end
               
           trans_rands(cntr+1:cntr2,:) = pmat_cum_ht(micro_state(cntr+1:cntr+theta1_cntr(ii,2),t-1),:)> rand(theta1_cntr(ii,2),1)*ones(1,N_mic_state);
           cntr = cntr2;
  end  % end of home-theta type loop (alternative to end below)
     micro_state(:,t) = int16(N_mic_state + 1 - sum(trans_rands,2)); % drawn new micro states
     cum_succ(:,t)    =  micro_state(:,t) - 1; % cumulative successes, new state, matrix for all firms starting (t+1)
%      cum_meets(:,t)   = Q_index(micro_state(:,t),2) - 1; % trials, new state, matrix for all firms (t+1)
   end   
    
 %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
    add_cli_cnt(:,t) = max(cum_succ(:,t) - cum_succ(:,t-1),0); % 
    % max() resets count to 0 for neg. differences to deal with new exporters 
     
    % identify first period in which a new exporter is active. First
   % new_firm(:,t) = max((cum_meets(:,t) - cum_meets(:,t-1) < 0),(stay == 0));     
    new_firm(:,t) = micro_state(:,t) == 1;  
    incumb        = ones(size(new_firm(:,t),1),1)- new_firm(:,t);     
    % Careful: incumb means no reset after previous period.  It is not the same
    % as cont_expr, which tracks whether an exporter had at least one
    % shipment during the previous year.
   
   
%% Deal with endogenous and exogenous drops (not in pmat_cum_h)

% identify z values at which exporters keep current matches from t-1 to t
  keep_cli = c_val_f_orig(:,prod_ndx,macro_state_h(t-1))' > 0; % = 1 if want to keep type for t
  drop_Zcut = N_Z - sum(keep_cli); % matches dropped at z value <= drop_Zcut
% count endogenous drops (z too low to continue) 
  drop_cnt = sum(lag_cli_zst.*(1-keep_cli),2);
  
% draw the number of exogenous deaths of remaining matches between t-1 and t 
  ddum = find(cur_cli_cnt(:,t-1)-drop_cnt > 0); 
   if sum(ddum)>0
     exog_deaths(ddum,t-1) =...
       random('bino',cur_cli_cnt(ddum,t-1)-drop_cnt(ddum),prob_mdeath);  
   end

%% update current count for new matches, drops, and exogenous deaths  
    cur_cli_cnt(:,t) = add_cli_cnt(:,t) + cur_cli_cnt(:,t-1) ...
                      - drop_cnt - exog_deaths(:,t-1) ; 

%% break down by buyer types (z)
                  
   for i=1:N_firms
    % break down new clients that occur between t-1 and t into z-types 
%      new_cli_zst(i,:) = new_vec(add_cli_cnt(i,t),N_Z,cum_pz); % distribute gross additions  

   new_cli_zst(i,:) = new_vec_C(add_cli_cnt(i,t),N_Z,cum_pz); % distribute gross additions
% The following patch was added because occasionally new_vec_C doesn't work right. Need to clean up  
   if sum(new_cli_zst(i,:)) ~= add_cli_cnt(i,t)
       'C routine for random draws failed in matchdat_gen_h. Trying Matlab version'
        new_cli_zst(i,:) = new_vec(add_cli_cnt(i,t),N_Z,cum_pz);
   end
       
    if exog_deaths(i,t-1) > 0
      % break down exogenous deaths that occur between t-1 and t down by z state: 
      die_cli_zst(i,:) = die_vec(lag_cli_zst(i,:).*keep_cli,exog_deaths(i,t-1),N_Z); 
      % record number of endogenous plus exogenous exits by z state in 1st column of trans_count:           
    end   
    trans_count(2:N_Z+1,1,i) = (lag_cli_zst(i,:).*(1-keep_cli))' + die_cli_zst(i,:)'; 
    % For each firm (i) of a particular type, column 1 of trans_mat now    
    % contains counts of all exiting matches, by buyer type (row).
    
    % Update surviving client counts by z type using transition matrix for
    % z. Do this for those that don't die for endogenous reasons, minus those 
    % that die for exogenous reasons:
     
    surviv_zst(i,:) = lag_cli_zst(i,:).*keep_cli - die_cli_zst(i,:); 
    N_sur = sum(surviv_zst(i,:),2); % number of survivors from t-1, firm i

     if N_sur > 0
        sur_typ = find(surviv_zst(i,:)); % addresses for z-states populated by at least one survivor
        for jj = sur_typ  % loop over initial states of surviving matches, exporter i
         draw = rand(surviv_zst(i,jj),1);  
       % identify destination z states for each surviving client (could be multiple survivors per initial type):
         trans_z = ones(surviv_zst(i,jj),1)*pmat_cum_z(jj,:) > draw;    
       % count # clients in each destination z state for each beginning z state. 
       % Record counts in cols 2:N_Z+1 of trans_mat. Rows are initial states, plus 1:
         trans_count(jj+1,2:N_Z+1,i) = sum(trans_z(:,1:N_Z) - [zeros(size(draw,1),1),trans_z(:,1:N_Z-1)],1); 
       % cumulate over b.o.p. z types to get vector of surviving client e.o.p. types. Rows (i) are exporters:   
         trans_zst(i,:) = trans_zst(i,:) +  trans_count(jj+1,2:N_Z+1,i); 
        end        
     end
    if sum(new_cli_zst(i,:),2)>0
        trans_count(1,2:N_Z+1,i) = new_cli_zst(i,:);
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % The following block is solely for de-bugging, and can be commented out

%    if i == N_firms  % matrices must fully loaded before tests
%     adds = cat(2,sum(new_cli_zst,2),add_cli_cnt(:,t))' ;
%     dies = cat(2,sum(die_cli_zst,2),exog_deaths(:,t-1))'; 
%     surv = cat(2,sum(surviv_zst,2),sum(trans_zst,2))' ;
%     zst_test = cat(2,cur_cli_cnt(:,t-1),reshape(sum(sum(trans_count(2:N_Z+1,:,:),2),1),N_firms,1))';
%    %     [type,prod_ndx,theta_ndx,t,season]
%       try     
%        assert(sum((adds(1,:)-adds(2,:)).^2)==0)
%       catch
%         [t,season,year,prod_ndx,theta_ndx,macro_state_h(t)]
%         warning('new client figures inconsistent: new_cli_zst versus add_cli_cnt')
%         adds
%       end
%       try
%          assert(sum((dies(1,:)-dies(2,:)).^2)==0)
%       catch
%          [t,season,year,prod_ndx,theta_ndx,macro_state_h(t)]
%          warning('dying client figures inconsistent: die_cli_zst versus exog_deaths')
%          dies
%       end
%       try
%          assert(sum((surv(1,:)-surv(2,:)).^2)==0)
%       catch
%          [t,season,year,prod_ndx,theta_ndx,macro_state_h(t)]
%          warning('survivor counts inconsistent: surviv_zst versus trans_zst')
%          surv
%       end
%       try
%         assert(sum((zst_test(1,:)-zst_test(2,:)).^2)==0)
%       catch
%          [t,season,year,prod_ndx,theta_ndx,macro_state_h(t)]
%          warning('transition counts inconsistent: cur_cli_cnt versus trans_count')
%          zst_test
%       end
%     end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   end
 cur_cli_zst = new_cli_zst + trans_zst; 
 
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Some more debugging stuff at this point has been moved to "debug_scraps_discrete_sim.m
%
%   assert(max(abs(cur_cli_cnt(:,t) - sum(cur_cli_zst,2)))==0)
%   assert(min(min(cur_cli_zst))==0)
% %  [cur_cli_cnt(:,t),reshape(sum(sum(trans_mat(:,2:N_Z+1,:,type),2),1),N_firms,1)]'
%   assert(max(abs(cur_cli_cnt(:,t) - reshape(sum(sum(trans_count(:,2:N_Z+1,:),2),1),N_firms,1)))==0) 
% %  [incumb.*add_cli_cnt(:,t),incumb.*reshape(sum(sum(trans_mat(1,2:N_Z+1,:,type),2),1),N_firms,1)]'
%   assert(max(abs(incumb.*add_cli_cnt(:,t) -...
%       incumb.*reshape(sum(sum(trans_count(1,2:N_Z+1,:),2),1),N_firms,1)))==0) ;
%  
% try
%   assert(max(abs(cur_cli_cnt(:,t) - sum(cur_cli_zst,2)))==0)
% catch
%    fprintf('\r Warning: max(abs(cur_cli_cnt(:,t) - sum(cur_cli_zst,2))) ~= 0 %.1f\n',...
%            max(abs(cur_cli_cnt(:,t) - sum(cur_cli_zst,2))));  
% end
% try
%   assert(min(min(cur_cli_zst))==0)
% %  [cur_cli_cnt(:,t),reshape(sum(sum(trans_mat(:,2:N_Z+1,:,type),2),1),N_firms,1)]'
% catch
%   fprintf('\r Warning: minimum cur_cli_zst ~= 0 %.1f\n', min(min(cur_cli_zst))); 
% end
% try
%   assert(max(abs(cur_cli_cnt(:,t) - reshape(sum(sum(trans_count(:,2:N_Z+1,:),2),1),N_firms,1)))==0) 
% %  [incumb.*add_cli_cnt(:,t),incumb.*reshape(sum(sum(trans_mat(1,2:N_Z+1,:,type),2),1),N_firms,1)]'
% catch
%     fprintf('\r Warning: current client counts messed up'); 
% end
% try
%   assert(max(abs(incumb.*add_cli_cnt(:,t) -...
%       incumb.*reshape(sum(sum(trans_count(1,2:N_Z+1,:),2),1),N_firms,1)))==0) ;
% catch
%     fprintf('\r Warning: client addition counts messed up');  
% end
 %% End diagnostics
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 %% Deal with dormant firms: after two years with no clients, swap them out
 if t > 2*pd_per_yr
  dmt = sum(cur_cli_cnt(:,t-2*pd_per_yr+1:t),2)==0; % identify dormant firms
  micro_state(dmt,t)  = 1;  % reset initial micro state to 1 (entrant)
  new_firm(dmt,t)     = 1;  % will mark firms that haven't made a match in 2 yrs
  exit_firm(dmt,t)    = 1;  % will mark last period of exiting firm (same thing)
  cum_meets(dmt,t)    = 0;  % cumulative number of meetings
  cum_succ(dmt,t)     = 0;  % cumulative number of successes
end
 
%% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% calculate time in domestic market (can do this outside t loop--see earlier version of code)
flr       = max(flrlag,t*new_firm(:,t));  % floor resets to current year for new firms. Hence
age       = t*ones(N_firms,1) - flr;    % age=0 all year for firms with no shipment previous year
flrlag    = flr ;
cumage    = cat(2,cumage,age);
yr_age    = floor(age.*frac_of_year);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%% Construct period-specific variables

%  First load season to season transitions into mat_tran, which describes  
%  matches of all N_firms of a particular type for a particular transition (t-1 to t).
mat_tran_all_zeros = ~any(trans_count(:));
if mat_tran_all_zeros
    mat_tran = zeros(0,4);ship_cur = zeros(0,1); age_vec = zeros(0,1);
else
[mat_tran,ship_cur,age_vec] =...
    match_sales(mm.scale_h,mm.eta,trans_count(:,:,:),age,mm.X_h(macro_state_h(t)),t,poisCDF,...
    max_ships,N_Z,mm.Z,mm.Phi(prod_ndx));
end
% mat_tran:  [initial state, exporter id, ending state, match revenue]

if size(mat_tran,1)>0 
    try
    % check that ending Z, if positive, is always greater than Zcut.
    assert(min(  (mat_tran(:,3)>0).*((mat_tran(:,3)>0) - drop_Zcut.*ones(size(mat_tran,1),1) )>0,[],1)>=0)
    catch
     warning('ending Z positive and less than Zcut')
     [mat_tran(:,3),drop_Zcut.*ones(size(mat_tran,1),1)];
    end
    if season > 1
        try
        % check that beginning Z, if positive, is always greater than last period's Zcut.
         assert(min((mat_tran(:,1)>0).*((mat_tran(:,1)>0) - seas_Zcut(season-1).*ones(size(mat_tran,1),1))>0,[],1)>=0)
        catch
         warning('beginning Z positive and less than last period Zcut')
         [mat_tran(:,1),seas_Zcut(season-1).*ones(size(mat_tran,1),1)]
        end
    end
end

  if season == 1
    N_match = size(mat_tran,1);
  end
  
  if t==2  
    seas_Zcut_lag = seas_Zcut;
    seas_tran_lag = seas_tran;
  end      % (will want to discard observations from first year)
  
      seas_tran{1,season} = [[t,season,year].*ones(size(mat_tran,1),1),mat_tran,ship_cur,age_vec];
      seas_Zcut(season)   = drop_Zcut;

    N_match = N_match + size(mat_tran,1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%% construct annualized variables  
    if season == pd_per_yr           
        
         [mat_yr_sales,firm_yr_sales] =...
           season_merge(seas_tran,N_match,N_firms,pd_per_yr);

% mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age w/in yr, firm age] 
% firm_yr_sales:[firm ID, total dom. sales, total dom. shipments, firm age in domestic market]
        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the following matrices accumulate annualized values over time and firm types
          theta_h_match = theta_h(mat_yr_sales(:,1));
          tt =  ones(size(mat_yr_sales,1),1).*[t,type];
          agg_mat_yr_sales  = [agg_mat_yr_sales;[tt,mat_yr_sales]];
       % agg_mat_yr_sales: [t,type,firm ID, match sales, shipments, boy Z, eoy Z, match age, firm age] 
          theta_h_firm = theta_h(firm_yr_sales(:,1));
          ttt = ones(size(firm_yr_sales,1),1).*[t,type]; 
          agg_firm_yr_sales = [agg_firm_yr_sales;[ttt,firm_yr_sales]]; 
          agg_theta_h_firm  = [agg_theta_h_firm;theta_h_firm]; % keep track of domestic thetas for each firm
       % agg_firm_yr_sales: [t,type,firm ID, total sales, total shipments,firm age]
       
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% update lagged variables for next year's loop         
         seas_tran_lag = seas_tran;
         seas_Zcut_lag = seas_Zcut;
         N_match_lag   = N_match;
         seas_Zcut = zeros(1,pd_per_yr);     

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%% Construct and cumulate moments        
  if year > mm.burn  % cosntruct moments for firm domestic sales regressions   
    % autoregressions and degree distribution       
          
      [x,y,fmoms_xx,fmoms_xy,fysum,fn_obs] = firm_reg_h_moms(firm_yr_sales,firm_yr_sales_lag,N_firms);
         
         agg_x_fsales_h   = [agg_x_fsales_h;x];
         agg_y_fsales_h   = [agg_y_fsales_h;y];
         agg_fmoms_h_xx = agg_fmoms_h_xx + fmoms_xx; % cumulate moments for home sales AR1
         agg_fmoms_h_xy = agg_fmoms_h_xy + fmoms_xy; % cumulate moments for home sales AR1
         agg_fysum_h    = agg_fysum_h + fysum;
         agg_fnobs_h    = agg_fnobs_h + fn_obs ;          
         
  end   % year > mm.burn if statement    
         mat_yr_sales_lag = mat_yr_sales;   % stack data for match regression
         firm_yr_sales_lag = firm_yr_sales; % stack data for firm regression
         
  end   % season == pd_per_yr if statement 
    
      season = season + 1;

%% load lagged client state matrix and re-initialize objects

      lag_cli_zst  = cur_cli_zst;
      new_cli_zst  = zeros(N_firms,N_Z);
      die_cli_zst  = zeros(N_firms,N_Z);  
      trans_zst    = zeros(N_firms,N_Z); 
      trans_count  = zeros(N_Z+1,N_Z+1,N_firms);

      tlag = t;
%   end   % end of home-theta type loop (alternative to end above)

  if toc > 50
     disp("ERROR: simulations taking too long in matchdat_gen_h")
     my_flag = 1;
     return
  end
  
%   if N_firms>=20  && sum(sum(cur_cli_cnt))>100 && t>=540
%     t
%     N_firms
%     'pause here'
%     
% end

end

% This function is called from discrete_sim_parfor3.m

% It simulates panel data on foreign sales, matches, and shipments for all firms of a
% particular type. It also generates type-specific moments to be aggregated
% in discrete_sim_parfor3.m. Type is determined by productivity and foreign theta.

function [singletons, agg_time_gaps,agg_mat_yr_sales,agg_mat_yr_sales_adj,agg_firm_yr_sales,...
          agg_mat_ar1_x,agg_mat_ar1_y,... 
          agg_moms_xx,agg_moms_xy,agg_ysum,agg_nobs,agg_fmoms_xx,agg_fmoms_xy,... 
          agg_fysum,agg_fnobs,agg_exit_moms_xx,agg_exit_moms_xy,agg_sum_succ_rate,... 
          agg_exit_obs,agg_sum_exits,agg_mat_exit_moms_xx,agg_mat_exit_moms_xy,...
          agg_mat_exit_x,agg_mat_exit_y,...
          agg_mat_obs,agg_nmat_exit,agg_ship_obs,agg_ln_ships,agg_match_count,my_flag]...
    = matchdat_gen_f(N_firms,N_types,N_Z,N_Phi,N_mic_state,pd_per_yr,frac_of_year,periods,...
           max_ships,typemat,macro_state_f,theta_ndx,prod_ndx,mm,pmat_cum_f,...
           c_val_f_orig,succ_prob,prod_lvl,cum_pz,poisCDF,prob_mdeath,prob_fdeath,...
           haz_ship,pmat_cum_z,nn2,nn4,lambda_f,max_match)  

    
%% Initialize matrices

    Q_size = nn2*(nn2+1)/2;    % was nn1 instead of nn2--both have same value, but need to check arguments of lambda_f
    Q_index = zeros(Q_size,3); % [index,trials,successes]
    
    counter = 1;
    for i=1:1:nn2    % number of meetings, plus 1
        for ss=1:1:i % number of successes, plus 1
            Q_index(counter,:) = [counter,i,ss];
            counter = counter + 1;
        end
    end 
    
  seas_tran = cell(1,pd_per_yr); % cells will hold one year's worth of season- and match-specific outcomes for all firms w/in type
  seas_Zcut = zeros(1,pd_per_yr);    % elements will hold season-specifics Z cut-offs for endog. drops
% Each firm begins with zero trials zero successes, macro state at median position
  lag_row = ones(N_firms,1);

  cur_cli_cnt  = zeros(N_firms,periods,1); % clients active in the current period
  add_cli_cnt  = zeros(N_firms,periods,1); % gross additions to client count
  cum_meets    = zeros(N_firms,periods,1); % cumulative number of meetings
  cum_succ     = zeros(N_firms,periods,1); % cumulative number of successes
  new_firm     = zeros(N_firms,periods);   % will mark new firms that haven't made a match
  exit_firm    = zeros(N_firms,periods);   % will mark last period of exiting firm 
%  tot_ships    = zeros(N_firms,periods,1); % total number of shipments
  exog_deaths  = zeros(N_firms,periods,1); % number of exogenous match deaths
  micro_state  = ones(N_firms,periods,1);  % scalar indices for #success/#meetings 
  cur_cli_zst  = zeros(N_firms,N_Z);  % breaks down current client counts by z state
  lag_cli_zst  = zeros(N_firms,N_Z);  % breaks down lagged clients by z state
  new_cli_zst  = zeros(N_firms,N_Z);  % breaks down new client counts by z state
  die_cli_zst  = zeros(N_firms,N_Z);  % breaks down client death counts by z state
  surviv_zst   = zeros(N_firms,N_Z);  % breaks down surviving client counts by z state
  trans_zst    = zeros(N_firms,N_Z);  % counts survival types by firm after z innovations
  drop_cnt     = zeros(N_firms,1);    % number of clients endogenously dropped
  flrlag       = ones(N_firms,1);     % initializing vector for age debugging
  yr_flrlag    = ones(N_firms,1);     % initializing vector for age debugging
  cumage       = zeros(N_firms,1);    % initializing vector for age debugging 
%  cont_expr    = ones(N_firms,1);    % continuing exporter indicator
  
% create first observation on firm-year level aggregates (will concatenate below)
% max_match         = 50; % upper bound on number of matches to be counted 
agg_match_count      = zeros(max_match,1);
agg_mat_yr_sales     = zeros(0,9);
agg_mat_yr_sales_adj = zeros(0,9);
agg_mat_matur        = zeros(0,7);
mat_matur_dat        = zeros(0,7);
agg_firm_yr_sales    = zeros(0,6);
agg_time_gaps        = zeros(0,7);
agg_mat_ar1_x        = zeros(0,4);
agg_mat_ar1_y        = zeros(0,1);
mkt_exit             = zeros(1,3);
mat_yr_sales_adj     = zeros(0,9);

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
agg_fmoms_xx = zeros(4,4);
agg_fmoms_xy = zeros(4,1);
agg_fysum    = 0;
agg_fnobs    = 0; 

agg_exit_moms_xx = zeros(6,6);
agg_exit_moms_xy = zeros(6,1);
agg_sum_succ_rate = 0;
agg_exit_obs = 0;
agg_sum_exits = 0;

agg_mat_exit_moms_xx  = zeros(5,5);
agg_mat_exit_moms_xy  = zeros(5,1);
agg_mat_obs       = 0;
agg_nmat_exit     = 0;
agg_mat_exit_x   = zeros(0,5);
agg_mat_exit_y   = zeros(0,1);

trans_count    = zeros(N_Z+1,N_Z+1,N_firms); % counts transitions across buyer types, 
% for each seller type. New buyer types are considered type 0 at beginning of period, hence the +1. 
% Exiting firms are considered to move to type 0 at the the end of the period.
% Dimensions: (1) initial z-state (2) new z-state (3) firm index, given type    

%  trans_tot    = zeros(N_Z+1,N_Z+1,periods,N_types); % transition counts through time, aggregeted across firms

%% Initialize stuff

% Set period 1 values for keep_cli. It will be used to select matches thar are endogenously dropped.
  keep_cli      = ones(1,N_Z); % applies to clients existing in period 1
  keep_cli(1:5) =  zeros(1,5); % implying worst 5 client types from period 1 are dropped.
   
  year = 1;
  tlag = 1;
  season = 1;
  N_match = 0;
  N_match_lag = 0;
  firm_yr_sales_lag = zeros(N_firms,4); 
  % firm_yr_sales_lag will contain: [firmID,sales,#shipments,firm age]

%% TIME LOOP BEGINS HERE
   tic
   my_flag = 0;
   for t = 2:1:periods

% update year, type and pmat_cum_t
    year = floor((t-1)/pd_per_yr);
    mic_type = find(typemat(:,3) == theta_ndx & typemat(:,4) == prod_ndx,1,'first');  
    type     = find(typemat(:,2) == macro_state_f(t) & typemat(:,3) == theta_ndx & typemat(:,4) == prod_ndx,1,'first');  

% Load relevant transition probabilites for current micro type & macro state.
% pmat_cum_t holds cumulative transition probs across (#success, #meeting) pairs.
% It's created in inten_sim using the relevant Q matrix.   
    pmat_cum_t = pmat_cum_f{type}; % type-specific cum trans. probs., micro state  

% reset season when previous period completes a year
     if abs(floor((t-1)/pd_per_yr) - (t-1)/pd_per_yr) <= 0.0001 
        season = 1;  
     end

%% gross additions to clients, before drops and deaths between t-1 and t

%  Need to treat firms that have maxed out learning separately from others. 
%  To get to the no-learning state, firms must first land on or above n = mm.n_size-3. 
%  where mm.n_size is the maximum number of meetings that firms learn from.
%  The -3 keeps firms transiting to the no-learning states from first having to
%  land first on mm.n_size. They can land on mm.n_size,mm.n_size-1,mm.n_size-2, or mm.n_size-3 
   stay     = ones(N_firms,1);  % will flag no-learning firms that continue, t-1 to t
   no_learn = cum_meets(:,t-1) >= mm.n_size-3; % for picking off seasoned exporters
   learn    = cum_meets(:,t-1) <  mm.n_size-3; % for picking off learning exporters
     
 % First deal with firms that have not maxed out learning using pmat_cum, based on Q matrix. 
 
 % Find N_learn randomly selected micro states next period, given 
 % macro state (common to all firms), initial micro states, and pmat_cum.
   N_learn = sum(learn,1) ;
   if  N_learn > 0
     trans_rands = pmat_cum_t(micro_state(learn,t-1),:)> rand(N_learn,1)*ones(1,N_mic_state);
     micro_state(learn,t) = int16(N_mic_state + 1 - sum(trans_rands,2)); % drawn new micro states
     cum_meets(learn,t)   = Q_index(micro_state(learn ,t),2) - 1; % trials, new state, matrix for all firms (t+1)
     cum_succ(learn,t)    = Q_index(micro_state(learn ,t),3) - 1; % successes, new state, matrix for all firms starting (t+1)
   end   

   stay(learn)  = micro_state(learn,t-1) - micro_state(learn,1) ~= 0; % wasn't in initial state last period

 % Now deal with firms that have maxed-out learning. 
 
 % Update cumulative meetings and successes for firms with >= mm.n_size meets, if
 % any. These firms are presumed to know their thetas, hence nn1 = floor(succ_prob*nn2) 
 % These calculations don't involve the Q matrix or the associated pmat trans. probs. 
    N_no_learn  = sum(no_learn);
    if N_no_learn >0        
    % no exog. death, and shipments prev. year
    stay(no_learn) = (rand(N_no_learn,1) > prob_fdeath);        
    cum_meets(no_learn,t) = (cum_meets(no_learn,t-1) + poissinv(rand(N_no_learn,1),succ_prob ...
        * reshape(lambda_f(floor(succ_prob*nn2),nn2,1,min(nn4,cum_succ(no_learn,t-1)+1),...
          prod_ndx,macro_state_f(t-1)),N_no_learn,1))).*stay(no_learn); % resets to 0 if stay==0 or new firm this period               
      
    cum_succ(no_learn,t) = cum_succ(no_learn,t-1).*stay(no_learn)...  %.*(1-new_firm(no_learn)) ...
        + random('bino',cum_meets(no_learn,t)-cum_meets(no_learn,t-1).*stay(no_learn),succ_prob);  
    % Note: cum_succ includes single shipment matches (with bad z's) that are dropped after 1 period. 
    % cum_succ and cum_meets will always be 0 after a period t reset due to exit (stay==0).   
    end
    
    %(succ,trial,common succ rate (defunct), network size, prod of firm, macro shock) 
    
 %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   
 % Hereafter can treat learning and no-learning firms together
    
    add_cli_cnt(:,t) = max(cum_succ(:,t) - cum_succ(:,t-1),0); % 
    % max() resets count to 0 for neg. differences to deal with new exporters 
     
    % identify first period in which a new exporter is active. First
    % condition is for learning firms; second (stay==0) is for no-learning
    
    new_firm(:,t)    = max((cum_meets(:,t) - cum_meets(:,t-1) < 0),(stay == 0));      
    exit_firm(:,t-1) = cum_meets(:,t) - cum_meets(:,t-1) < 0 ; % last period before exit
    
    % NOTE: new_firm = 1 each period that stay = 0. 
    
    incumb        = ones(size(new_firm(:,t),1),1)- new_firm(:,t);     
    % Careful: incumb means no reset after previous period.  It is not the same
    % as cont_expr, which tracks whether an exporter had at least one
    % shipment during the previous year.
   
   
%% Deal with endogenous and exogenous drops

% identify z values at which exporters keep current matches from t-1 to t
  keep_cli = c_val_f_orig(:,prod_ndx,macro_state_f(t-1))' > 0; % = 1 if want to keep type for t 
  drop_Zcut = N_Z - sum(keep_cli); % cutoff: matches dropped at z value <= drop_Zcut
  
% count endogenous drops for all exporter hotel rooms (exporters): z too low to continue
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
    % break down new clients that occur between t-1 and t into e.o.p. z-types 
      new_cli_zst(i,:) = new_vec_C(add_cli_cnt(i,t),N_Z,cum_pz); % distribute gross additions 
      if sum(new_cli_zst(i,:)) ~= add_cli_cnt(i,t) % cheap patch--better to clean up C++ code
       'C routine for random draws failed in matchdat_gen_f. Trying Matlab version'
      new_cli_zst(i,:) = new_vec(add_cli_cnt(i,t),N_Z,cum_pz);
      end
    if exog_deaths(i,t-1) > 0
      % break down exogenous deaths that occur between t-1 and t down by b.o.p. z state: 
      die_cli_zst(i,:) = die_vec(lag_cli_zst(i,:).*keep_cli,exog_deaths(i,t-1),N_Z);            
    end   
    %trans_count_test = trans_count;
    trans_count(2:N_Z+1,1,i) = (lag_cli_zst(i,:).*(1-keep_cli))' + die_cli_zst(i,:)';
    % For each firm (i) of a particular type, column 1 of trans_count(:,:,i)    
    % now contains counts of all exiting matches (endog. and exog.), by buyer type (row).
    
    % Update surviving client counts by z type using transition matrix for z. 
    % Do this for those that don't die for endogenous or exogenous reasons. 
    surviv_zst(i,:) = lag_cli_zst(i,:).*keep_cli - die_cli_zst(i,:); 
    N_sur = sum(surviv_zst(i,:),2); % number of survivors from t-1 by b.o.p. type, firm i
    
     if N_sur > 0
        sur_typ = find(surviv_zst(i,:)); % addresses for z-states populated by at least one survivor, firm i
        for jj = sur_typ  % loop over initial (b.o.p.) states of surviving matches, exporter i
         draw = rand(surviv_zst(i,jj),1);  
       % identify destination z states for each surviving client (could be multiple survivors per initial type):
         trans_z = ones(surviv_zst(i,jj),1)*pmat_cum_z(jj,:) > draw;    
       % count # clients in each destination z state for each beginning z state. 
       % Record e.o.p. counts in cols 2:N_Z+1 of trans_count. Row indices are b.o.p. states, plus 1:
         trans_count(jj+1,2:N_Z+1,i) = sum(trans_z(:,1:N_Z) - [zeros(size(draw,1),1),trans_z(:,1:N_Z-1)],1); 
       % cumulate over b.o.p. z types to get row vector of surviving client e.o.p. types. Rows (i) index exporter hotel rooms:   
         trans_zst(i,:) = trans_zst(i,:) +  trans_count(jj+1,2:N_Z+1,i); 
        end        
     end
    if sum(new_cli_zst(i,:),2)>0
        trans_count(1,2:N_Z+1,i) = new_cli_zst(i,:); % load new clients for exporter i in first row of trans_count
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
%         [t,season,year,prod_ndx,theta_ndx,macro_state(t)]
%         warning('new client figures inconsistent: new_cli_zst versus add_cli_cnt')
%         adds
%       end
%       try
%          assert(sum((dies(1,:)-dies(2,:)).^2)==0)
%       catch
%          [t,season,year,prod_ndx,theta_ndx,macro_state(t)]
%          warning('dying client figures inconsistent: die_cli_zst versus exog_deaths')
%          dies
%       end
%       try
%          assert(sum((surv(1,:)-surv(2,:)).^2)==0)
%       catch
%          [t,season,year,prod_ndx,theta_ndx,macro_state(t)]
%          warning('survivor counts inconsistent: surviv_zst versus trans_zst')
%          surv
%       end
%       try
%         assert(sum((zst_test(1,:)-zst_test(2,:)).^2)==0)
%       catch
%          [t,season,year,prod_ndx,theta_ndx,macro_state(t)]
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
% 
% try
%   assert(min(min(cur_cli_zst))>=0)
% %  [cur_cli_cnt(:,t),reshape(sum(sum(trans_mat(:,2:N_Z+1,:,type),2),1),N_firms,1)]'
% catch
%   fprintf('\r Warning: minimum cur_cli_zst < 0 %.1f\n', min(min(cur_cli_zst))); 
% %  'pause here'
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
 
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 
%% Deal with dormant firms: if a firm has no clients for more than one year,   
%  kick it out of its slot and start a new firm

if t > pd_per_yr
  dmt = sum(cur_cli_cnt(:,t-pd_per_yr+1:t),2)==0; % identify dormant firms
  micro_state(dmt,t)  = 1;  % reset initial micro state to 1 (entrant)
  new_firm(dmt,t)     = 1;  % will mark new firms that haven't made a match
  exit_firm(dmt,t)    = 1;  % will mark last period of exiting firm (same thing)
  cum_meets(dmt,t)    = 0;  % cumulative number of meetings
  cum_succ(dmt,t)     = 0;  % cumulative number of successes
end

 %% End diagnostics
    
%% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Calculate time (in periods) in export market  
  flr       = max(flrlag,t*new_firm(:,t)); % floor resets to current year for new exporters. 
  age       = t*ones(N_firms,1) - flr;     % age in periods. age=0 all year for firms with no shipment previous year
  flrlag    = flr ; % carry floor forward for continuing matches
  cumage    = cat(2,cumage,age);

  % Calculate time (in years) in export market 
%   yr_age    = floor(age.*frac_of_year);
%   yr_flr    = max(yr_flrlag, year*new_firm(:,t));
%   yr_age    = year*ones(N_firms,1) - yr_flr;
%   yr_flrlag = yr_flr;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%% Construct period-specific variables

%  Load season to season transitions into mat_tran, which describes matches  
%  of all N_firms of a particular type for a particular transition (t-1 to t).

mat_tran_all_zeros = ~any(trans_count(:));
if mat_tran_all_zeros
    mat_tran = zeros(0,4);ship_cur = zeros(0,1); age_vec = zeros(0,1);
else
[mat_tran,ship_cur,age_vec] =...
    match_sales(mm.scale_f,mm.eta,trans_count(:,:,:),age,mm.X_f(macro_state_f(t)),t,poisCDF,...
    max_ships,N_Z,mm.Z,mm.Phi(prod_ndx));
% mat_tran:  [initial state, exporter id, ending state, match revenue]
% ship_cur:   match's number of shipments within the current period
% age_vec:    firm age (# periods)
end

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
  end % needed to get started (will be discarded with the burn-in)
  
 % Concatenate time index, season index and year index with match variables 
 % and collect results for all seasons in year in seas_tran
    seas_tran{1,season} = [[t,season,year].*ones(size(mat_tran,1),1),mat_tran,ship_cur,age_vec];
 % seas_tran: [t, season, year, initial state, exporter id, ending state, match revenue,
 %             # shipments, firm age (# periods)]    
    seas_Zcut(season)   = drop_Zcut;

    N_match = N_match + size(mat_tran,1); % cumulate match count within current year

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%% construct annualized variables  
   if season == pd_per_yr           
       
         [mat_yr_sales,firm_yr_sales] =...
           season_merge(seas_tran,N_match,N_firms,pd_per_yr);
% mat_yr_sales:  [firm ID, match-specific sales, shipments, boy Z, eoy Z, 
%                 match age in periods (w/in year), firm age in periods] 
% firm_yr_sales: [firmID,sales,#shipments,firm age]
        
% convert age in periods to age in years for first-year observations.
% After the first year, conversion to years handled by mat_yr_splice
if year==1
  mat_yr_sales(:,6)  = mat_yr_sales(:,6) > 0;  % set match age in years to 1 if match age in periods > 0 
  mat_yr_sales(:,7)  = mat_yr_sales(:,7) > 0;  % set firm age in years to 1 if year age in periods > 0 
  firm_yr_sales(:,4) = firm_yr_sales(:,4) > 0; % set firm age in years to 1 if year age in periods > 0 
  year_lag           = 1;
end 

% # unsuccessful meetings (duds) over previous year, by firm (needed for degree distribution later)
        yr_tlag = t-pd_per_yr;
        cum_duds  = cum_meets(:,yr_tlag:t) - cum_succ(:,yr_tlag:t); % previous pd_per_yr + 1 cumulative duds
        curr_duds = cum_duds(:,2:pd_per_yr+1)-(new_firm(:,yr_tlag+1:t)==0).*cum_duds(:,1:pd_per_yr);

        try
        assert(min(min(curr_duds))>=0) 
        catch
           warning('firm turnover not handled correctly')
           [row_temp,~] = find(curr_duds<0);
           curr_duds(row_temp,:)
           [t, mic_type, season, size(curr_duds)]
           curr_duds = max(curr_duds,zeros(size(curr_duds)));
        end
        
        singletons = sum(sum(curr_duds));
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  Generate observations on intervals between meetings, cumulative meetings, and cumulative succcesses

        if t>3*pd_per_yr
          [time_gap,mkt_exit] = time_gaps(t,exit_firm,pd_per_yr,cum_meets,cum_succ);
%           time_gap: (1) firm_ID, (2) periods into interval, (3) time gap, 
%           (4) # new meetings, (5) t, (6) cum. meetings, (7) cum succeses

          agg_time_gaps = [agg_time_gaps;time_gap];
        end
             
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% update lagged variables for next year's loop         
         seas_tran_lag = seas_tran;
         seas_Zcut_lag = seas_Zcut;
         N_match_lag   = N_match;
         seas_Zcut = zeros(1,pd_per_yr);     

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
  if year > 1  % construct moments for match regressions, firm export regressions   

  if size(mat_yr_sales,1)*size(mat_yr_sales_lag)>0
  % drops types without sales in both current and lagged years
      
 % mat_yr_splice splices consecutive obs. on annualized data and updates 
 % match ages and firm ages, converting them to years. It's really costly! 
 
    ncols = size(mat_yr_sales,2);  

    [mat_cont_2yr,mat_yr_sales,mat_yr_sales_adj,year_lag] =...
        mat_yr_splice_v2(mat_yr_sales,mat_yr_sales_lag,mm,year_lag,year);    
    
    %  mat_cont_2yr: [mat_yr_sales_lag(ff_cont_lag,:), mat_yr_sales(ff_cont,:)]
    %  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age in yrs, firm age in yrs] 
    %  mat_yr_sales_adj: same as lagged mat_yr_sales except eoy Z set to zero if no sales next year
    
      try
        assert(sum(abs(mat_cont_2yr(:,5)-mat_cont_2yr(:,ncols+4)))==0) % Z lines up  
        assert(sum(abs(mat_cont_2yr(:,1)-mat_cont_2yr(:,ncols+1)))==0) % firm ID lines up
%       assert(sum(mat_cont_2yr(:,14)- mat_cont_2yr(:,13)<0)==0) % firm age never less than match age
      catch
       'pause here'
        warning('splicing mismatch or firm age < match age')
      end 

%     if sum(mat_cont_2yr(:,14)- mat_cont_2yr(:,13)<0)>0 % print obs if firm age less than match age
%         fff = mat_cont_2yr(:,14)- mat_cont_2yr(:,13)<0;
%         'printing from line 526 in matchdat_gen_f'
%         [ones(sum(fff),1)*[t,macro_state_f(t),theta_ndx,prod_ndx], mat_cont_2yr(fff,6:7), mat_cont_2yr(fff,13:14)]
%     end  
    
    
    %% the following matrices accumulate annualized values over time and firm types
 
   mat_matur_dat =  [mat_yr_sales(:,2),mat_yr_sales(:,4:7),mat_yr_sales(:,1),ones(size(mat_yr_sales,1),1).*t/pd_per_yr] ; 
   % agg_mat_matur: [sales, boy Z, eoy Z, match age, firm age, firm_ID, yr] 
            
  if year > mm.burn  % don't start building simulated data set until burn-in finished              
      
       tt =  ones(size(mat_yr_sales,1),1).*[t,mic_type]; % add cols 1 and 2: t, firm type
       agg_mat_yr_sales  = [agg_mat_yr_sales;[tt,mat_yr_sales]];
       % agg_mat_yr_sales: [t,type,firm ID, match sales, shipments, boy Z, eoy Z, w/in yr. match age, firm age] 

   if year_lag == year % check that mat_yr_splice_v2 ran & updated mat_yr_sales_adj 
       tt2 =  ones(size(mat_yr_sales_adj,1),1).*[t-pd_per_yr,mic_type];
          agg_mat_yr_sales_adj = [agg_mat_yr_sales_adj;[tt2,mat_yr_sales_adj]]; % add cols 1 and 2: t, firm type 
       % agg_mat_yr_sales_adj: [t,type,firm ID, match sales, shipments, boy Z, adj_eoy Z, w/in yr. match age, firm age] 
          ttt = ones(size(firm_yr_sales,1),1).*[t,mic_type]; 
          agg_firm_yr_sales = [agg_firm_yr_sales;[ttt,firm_yr_sales]]; % add cols 1 and 2: t, firm type
       % agg_firm_yr_sales: [t,type,firm ID, total exports,total shipments,firm age]
   end
   
   agg_mat_matur =  [agg_mat_matur; mat_matur_dat]; 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  match regression moments

      [mat_ar1_x,mat_ar1_y,moms_xx,moms_xy,ysum,n_obs] = match_reg_moms(mat_cont_2yr,ncols);
          
         agg_moms_xx = agg_moms_xx + moms_xx; % cumulate moments for match regression
         agg_moms_xy = agg_moms_xy + moms_xy; % cumulate moments for match regression
         agg_ysum    = agg_ysum + ysum;
         agg_nobs    = agg_nobs + n_obs;  
         
         agg_mat_ar1_x = [agg_mat_ar1_x;mat_ar1_x];
         agg_mat_ar1_y = [agg_mat_ar1_y;mat_ar1_y];
  end 
  
  end       
         
    if year >= mm.burn  % don't start building simulated data set until burn-in finished 
  
    % autoregressions and degree distribution       
          
      [fmoms_xx,fmoms_xy,fysum,fn_obs] = firm_reg_moms(firm_yr_sales,firm_yr_sales_lag,N_firms);
         
         agg_fmoms_xx = agg_fmoms_xx + fmoms_xx; % cumulate moments for firm regression
         agg_fmoms_xy = agg_fmoms_xy + fmoms_xy; % cumulate moments for firm regression
         agg_fysum    = agg_fysum + fysum;
         agg_fnobs    = agg_fnobs + fn_obs ; 
         
    % foreign market exit regression moments
    
      ff_exit = mkt_exit(:,2)>0;  % column 2 of mkt_exit is number of meetings
      if sum(ff_exit,1)>0
      [exit_moms_xx,exit_moms_xy,sum_succ_rate,sum_exits,exit_obs] = mkt_exit_moms(mkt_exit);
      
         agg_exit_moms_xx = agg_exit_moms_xx + exit_moms_xx;
         agg_exit_moms_xy = agg_exit_moms_xy + exit_moms_xy;
         agg_sum_succ_rate = agg_sum_succ_rate + sum_succ_rate;
         agg_exit_obs = agg_exit_obs + exit_obs;
         agg_sum_exits = agg_sum_exits + sum_exits;
      end
    % match exit regression moments
      if year_lag == year
        ff_mexit = mat_yr_sales_adj(:,2)>0; 
        if sum(ff_mexit,1)>0 % positive sales for at least one match      
            [mat_exit_x,mat_exit_y,mat_exit_moms_xx,mat_exit_moms_xy,mat_obs,nmat_exit]...
                = match_exit_moms(mat_yr_sales_adj(ff_mexit,:),pd_per_yr);
        
% Notes on variables:

%         mat_exit_y  = matches(ff,5)==0;         % match dead by end of year         
%         mat_exit_x = [x0,x1,x2,x3,x4];
%           x0 = ones(size(ff,1),1);              % ff picks off non-missing obs.
%           x1 = matches(ff,4)==0;                % first year dummy
%           x2 = log(matches(ff,2));              % sales during year
%           x3 = log(1+matches(ff,6)./pd_per_yr); % age of match
%           x4 = log(1+matches(ff,7));            % age of exporter

            agg_mat_exit_moms_xx = agg_mat_exit_moms_xx + mat_exit_moms_xx;
            agg_mat_exit_moms_xy = agg_mat_exit_moms_xy + mat_exit_moms_xy;
            agg_mat_obs      = agg_mat_obs + mat_obs;
            agg_nmat_exit    = agg_nmat_exit + nmat_exit;
            agg_mat_exit_x   = [agg_mat_exit_x;mat_exit_x];
            agg_mat_exit_y   = [agg_mat_exit_y;mat_exit_y];
        end
      end      
    % shipment and match counter
    
      [nship_obs,ln_ships,match_count] = match_shpt_cntr(mat_yr_sales_adj,max_match);
    
        agg_ship_obs    = agg_ship_obs + nship_obs ;
        agg_ln_ships    = agg_ln_ships + ln_ships ;
        agg_match_count = agg_match_count + match_count ;
        % include all the matches that generated a single sample shipment:
        agg_match_count(1) = agg_match_count(1) + singletons;
         
  end   % year > mm.burn  if statement    
   end
         mat_yr_sales_lag = mat_yr_sales;   % update lags
         firm_yr_sales_lag = firm_yr_sales; 

         
  end   % season == pd_per_yr if statement 
    
      season = season + 1;

%% load lagged client state matrix and re-initialize objects

      lag_cli_zst  = cur_cli_zst;
      new_cli_zst  = zeros(N_firms,N_Z);
      die_cli_zst  = zeros(N_firms,N_Z);  
      trans_zst    = zeros(N_firms,N_Z); 
      trans_count  = zeros(N_Z+1,N_Z+1,N_firms);

      tlag = t;
      
      
  if toc > 50
     disp("ERROR: simulations taking too long in matchdat_gen_f")
     my_flag = 1;
     return
  end

  
   end     % end of time loop
end

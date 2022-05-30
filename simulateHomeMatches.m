function iter_out = simulateHomeMatches(pt_ndx,macro_state_h,mm,policy,iter_out)  

% iter_in = struct;

% create home theta draws
theta_df = (length(mm.theta1)+1)*ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1) - ...
    sum(ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1)*mm.th1_cdf >...
    rand(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1)*ones(1,size(mm.theta1,2)),2);

% list the non-zero values and their frequencies
[uv,~,idx]   = unique(theta_df);
theta1_cntr  = [uv,accumarray(idx(:),1)];
theta_h      = sortrows(theta_df);

%% Initialize matrices

seas_tran = cell(1,mm.pd_per_yr); % cells will hold one year's worth of season- and match-specific outcomes for all firms w/in type
seas_Zcut = zeros(1,mm.pd_per_yr);    % elements will hold season-specifics Z cut-offs for endog. drops

cur_cli_cnt  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % clients active in the current period
add_cli_cnt  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % gross additions to client count
cum_succ     = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % cumulative number of successes
new_firm     = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods);   % will mark first period of a new firm
exog_deaths  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % number of exogenous match deaths
micro_state  = ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % scalar indices for #success/#meetings
lag_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % breaks down lagged clients by z state
new_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % breaks down new client counts by z state
die_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % breaks down client death counts by z state
surviv_zst   = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % breaks down surviving client counts by z state
trans_zst    = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % counts survival types by firm after z innovations
flrlag       = ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1);     % initializing vector for age debugging
cumage       = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1);    % initializing vector for age debugging

% create first observation on firm-year level aggregates (will concatenate below)
iter_out.firm_h_yr_sales = double.empty(0,6);
iter_out.theta_h_firm  = double.empty(0,1);
iter_out.x_fsales_h    = zeros(0,2); % will contain rows of x matrix for regression
iter_out.y_fsales_h    = zeros(0,1); % will contain rows of y matrix for regression

% firm level moment aggregators
iter_out.fmoms_h_xx = zeros(2,2);
iter_out.fmoms_h_xy = zeros(2,1);
iter_out.fysum_h    = 0;
iter_out.fnobs_h    = 0;

trans_count    = zeros(length(mm.Z)+1,length(mm.Z)+1,mm.sim_firm_num_by_prod_succ_type(pt_ndx)); % counts transitions across buyer types,
% for each seller type. New buyer types are considered type 0 at beginning of period, hence the +1.
% Exiting firms are considered to move to type 0 at the the end of the period.
% columns: (1) initial z-state (2) new z-state (3) firm index, given type (4) firm type

%% Initialize stuff

% keep_cli will be used to select matches that are endogenously dropped.
%keep_cli      = ones(1,size(mm.Z,1)); % applies to clients existing in period 1
%keep_cli(1:5) =  zeros(1,5); % implying worst 5 client types from period 1 are dropped

%% TIME LOOP BEGINS HERE
season = 1;
N_match = 0;
firm_yr_sales_lag = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),4);
% firm_yr_sales_lag will contain: [firmID,sales,#shipments,firm age]

tic
iter_out.abort_flag_h = 0;
for t = 2:1:mm.periods

    % update year, type and pmat_cum_t
    year = floor((t-1)/mm.pd_per_yr);

    % reset season when previous period completes a year
    if abs(floor((t-1)/mm.pd_per_yr) - (t-1)/mm.pd_per_yr) <= 0.001
        season = 1;
    end

    %% gross additions to clients, before drops and deaths between t-1 and t

    if  mm.sim_firm_num_by_prod_succ_type(pt_ndx) > 0
        trans_rands = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.nn_h);
        cntr = 0;

        for ii = 1:size(theta1_cntr,1)
            cntr2 = cntr+theta1_cntr(ii,2);
            type = find(policy.firm_type_prod_succ_macro(:,2) == macro_state_h(t)...
                & policy.firm_type_prod_succ_macro(:,3) == theta1_cntr(ii,1)...
                & policy.firm_type_prod_succ_macro(:,4) == mm.pt_type(pt_ndx,1),1,'first');
            pmat_cum_ht = policy.pmat_cum_h{type};
            
            trans_rands(cntr+1:cntr2,:) = pmat_cum_ht(micro_state(cntr+1:cntr2,t-1),:)...
                > rand(theta1_cntr(ii,2),1)*ones(1,mm.nn_h);
            cntr = cntr2;
        end  % end of home-theta type loop (alternative to end below)
        micro_state(:,t) = int16(mm.nn_h + 1 - sum(trans_rands,2)); % drawn new micro states
        cum_succ(:,t)    =  micro_state(:,t) - 1; % cumulative successes, new state, matrix for all firms starting (t+1)
    end

    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    add_cli_cnt(:,t) = max(cum_succ(:,t) - cum_succ(:,t-1),0); %
    % max() resets count to 0 for neg. differences to deal with new exporters

    % identify first period in which a new exporter is active. First
    new_firm(:,t) = micro_state(:,t) == 1;

    %% Deal with endogenous and exogenous drops (not in policy.pmat_cum_h)

    % identify z values at which exporters keep current matches from t-1 to t
    keep_cli = policy.c_val_h(:,mm.pt_type(pt_ndx,1),macro_state_h(t-1))' > 0; % = 1 if want to keep type for t
    drop_Zcut = size(mm.Z,1) - sum(keep_cli); % matches dropped at z value <= drop_Zcut
    % count endogenous drops (z too low to continue)
    drop_cnt = sum(lag_cli_zst.*(1-keep_cli),2);

    % draw the number of exogenous deaths of remaining matches between t-1 and t
    ddum = find(cur_cli_cnt(:,t-1)-drop_cnt > 0);
    if sum(ddum)>0
        exog_deaths(ddum,t-1) =...
            random('bino',cur_cli_cnt(ddum,t-1)-drop_cnt(ddum),(1-exp(-mm.delta)));
    end

    %% update current count for new matches, drops, and exogenous deaths
    cur_cli_cnt(:,t) = add_cli_cnt(:,t) + cur_cli_cnt(:,t-1) ...
        - drop_cnt - exog_deaths(:,t-1) ;

    %% break down by buyer types (z)

    for i=1:mm.sim_firm_num_by_prod_succ_type(pt_ndx)
        % break down new clients that occur between t-1 and t into z-types
        %      new_cli_zst(i,:) = new_vec(add_cli_cnt(i,t),size(mm.Z,1),cumsum(mm.erg_pz)); % distribute gross additions

        new_cli_zst(i,:) = new_vec_C(add_cli_cnt(i,t),size(mm.Z,1),cumsum(mm.erg_pz)); % distribute gross additions
        % The following patch was added because occasionally new_vec_C doesn't work right. Need to clean up
        if sum(new_cli_zst(i,:)) ~= add_cli_cnt(i,t)
            'C routine for random draws failed in matchdat_gen_h. Trying Matlab version'
            new_cli_zst(i,:) = new_vec(add_cli_cnt(i,t),size(mm.Z,1),cumsum(mm.erg_pz));
        end

        if exog_deaths(i,t-1) > 0
            % break down exogenous deaths that occur between t-1 and t down by z state:
            die_cli_zst(i,:) = createDieVec(lag_cli_zst(i,:).*keep_cli,exog_deaths(i,t-1),size(mm.Z,1));
            % record number of endogenous plus exogenous exits by z state in 1st column of trans_count:
        end
        trans_count(2:size(mm.Z,1)+1,1,i) = (lag_cli_zst(i,:).*(1-keep_cli))' + die_cli_zst(i,:)';
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
                trans_z = ones(surviv_zst(i,jj),1)*policy.pmat_cum_z(jj,:) > draw;
                % count # clients in each destination z state for each beginning z state.
                % Record counts in cols 2:size(mm.Z,1)+1 of trans_mat. Rows are initial states, plus 1:
                trans_count(jj+1,2:size(mm.Z,1)+1,i) = sum(trans_z(:,1:size(mm.Z,1)) - [zeros(size(draw,1),1),trans_z(:,1:size(mm.Z,1)-1)],1);
                % cumulate over b.o.p. z types to get vector of surviving client e.o.p. types. Rows (i) are exporters:
                trans_zst(i,:) = trans_zst(i,:) +  trans_count(jj+1,2:size(mm.Z,1)+1,i);
            end
        end
        if sum(new_cli_zst(i,:),2)>0
            trans_count(1,2:size(mm.Z,1)+1,i) = new_cli_zst(i,:);
        end
    end
    cur_cli_zst = new_cli_zst + trans_zst;

    %% Deal with dormant firms: after two years with no clients, swap them out
    if t > 2*mm.pd_per_yr
        dmt = sum(cur_cli_cnt(:,t-2*mm.pd_per_yr+1:t),2)==0; % identify dormant firms
        micro_state(dmt,t)  = 1;  % reset initial micro state to 1 (entrant)
        new_firm(dmt,t)     = 1;  % will mark firms that haven't made a match in 2 yrs
        cum_succ(dmt,t)     = 0;  % cumulative number of successes
    end

    %% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    % calculate time in domestic market (can do this outside t loop--see earlier version of code)
    flr       = max(flrlag,t*new_firm(:,t));  % floor resets to current year for new firms. Hence
    age       = t*ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1) - flr;    % age=0 all year for firms with no shipment previous year
    flrlag    = flr ;
    cumage    = cat(2,cumage,age);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Construct period-specific variables

    %  First load season to season transitions into mat_tran, which describes
    %  matches of all mm.sim_firm_num_by_prod_succ_type(pt_ndx) of a particular type for a particular transition (t-1 to t).
    mat_tran_all_zeros = ~any(trans_count(:));
    if mat_tran_all_zeros
        mat_tran = zeros(0,4);ship_cur = zeros(0,1); age_vec = zeros(0,1);
    else

        mkt = 2; % =2 for domestic market
        [mat_tran,ship_cur,age_vec] = match_sales(mkt,mm,trans_count,age,pt_ndx,macro_state_h(t));

    end
    % mat_tran:  [initial state, exporter id, ending state, match revenue]

    if season == 1
        N_match = size(mat_tran,1);
    end

    seas_tran{1,season} = [[t,season,year].*ones(size(mat_tran,1),1),mat_tran,ship_cur,age_vec];
    seas_Zcut(season)   = drop_Zcut;

    N_match = N_match + size(mat_tran,1);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% construct annualized variables
    if season == mm.pd_per_yr

        [~,firm_yr_sales] =...
            season_merge(seas_tran,N_match,mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.pd_per_yr);

        % firm_yr_sales:[firm ID, total dom. sales, total dom. shipments, firm age in domestic market]

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% the following matrices accumulate annualized values over time and firm types
        theta_h_firm = theta_h(firm_yr_sales(:,1));
        ttt = ones(size(firm_yr_sales,1),1).*[t,type];
        iter_out.firm_h_yr_sales = [iter_out.firm_h_yr_sales;[ttt,firm_yr_sales]];
        iter_out.theta_h_firm  = [iter_out.theta_h_firm;theta_h_firm]; % keep track of domestic thetas for each firm
        % iter_out.firm_h_yr_sales: [t,type,firm ID, total sales, total shipments,firm age]

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        seas_Zcut = zeros(1,mm.pd_per_yr);

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Construct and cumulate moments
        if year > mm.burn  % cosntruct moments for firm domestic sales regressions
            % autoregressions and degree distribution

            [x,y,fmoms_xx,fmoms_xy,fysum,fn_obs] = firm_reg_h_moms(firm_yr_sales,firm_yr_sales_lag,mm.sim_firm_num_by_prod_succ_type(pt_ndx));

            iter_out.x_fsales_h   = [iter_out.x_fsales_h;x];
            iter_out.y_fsales_h   = [iter_out.y_fsales_h;y];
            iter_out.fmoms_h_xx = iter_out.fmoms_h_xx + fmoms_xx; % cumulate moments for home sales AR1
            iter_out.fmoms_h_xy = iter_out.fmoms_h_xy + fmoms_xy; % cumulate moments for home sales AR1
            iter_out.fysum_h    = iter_out.fysum_h + fysum;
            iter_out.fnobs_h    = iter_out.fnobs_h + fn_obs ;

        end   % year > mm.burn if statement
        firm_yr_sales_lag = firm_yr_sales; % stack data for firm regression

    end   % season == mm.pd_per_yr if statement

    season = season + 1;

    %% load lagged client state matrix and re-initialize objects

    lag_cli_zst  = cur_cli_zst;
    new_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    die_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    trans_zst    = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    trans_count  = zeros(size(mm.Z,1)+1,size(mm.Z,1)+1,mm.sim_firm_num_by_prod_succ_type(pt_ndx));

end

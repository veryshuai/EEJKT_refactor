function [iterH_in, iter_out] = simulateHomeInnerInitialize(mm, pt_ndx, macro_state_h, iter_out)

%% Initialize matrices

iterH_in = struct;

iterH_in.pt_ndx    = pt_ndx;
iterH_in.macro_state_h = macro_state_h;
iterH_in.seas_tran = cell(1,mm.pd_per_yr);  % cells will hold one year's worth of season- and match-specific outcomes for all firms w/in type
iterH_in.seas_Zcut = zeros(1,mm.pd_per_yr); % elements will hold season-specifics Z cut-offs for endog. drops

iterH_in.cur_cli_cnt  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % clients active in the current period
iterH_in.add_cli_cnt  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % gross additions to client count
iterH_in.cum_succ     = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % cumulative number of successes
iterH_in.new_firm     = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods);   % will mark first period of a new firm
iterH_in.exog_deaths  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % number of exogenous match deaths
iterH_in.micro_state  = ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % scalar indices for #success/#meetings
iterH_in.lag_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % breaks down lagged clients by z state
iterH_in.new_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % breaks down new client counts by z state
iterH_in.die_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % breaks down client death counts by z state
iterH_in.surviv_zst   = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % breaks down surviving client counts by z state
iterH_in.trans_zst    = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % counts survival types by firm after z innovations
iterH_in.flrlag       = ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1);     % initializing vector for age debugging
iterH_in.cumage       = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1);    % initializing vector for age debugging

iterH_in.trans_count    = zeros(size(mm.Z,1)+1,size(mm.Z,1)+1,mm.sim_firm_num_by_prod_succ_type(pt_ndx)); % counts transitions across buyer types,
% for each seller type. New buyer types are considered type 0 at beginning of period, hence the +1.
% Exiting firms are considered to move to type 0 at the the end of the period.
% columns: (1) initial z-state (2) new z-state (3) firm index, given type (4) firm type

% create home theta draws
theta_df = (size(mm.theta1,2)+1)*ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1) - ...
    sum(ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1)*mm.th1_cdf >...
    rand(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1)*ones(1,size(mm.theta1,2)),2);

% list the non-zero values and their frequencies
[uv,~,idx]   = unique(theta_df);
iterH_in.theta1_cntr  = [uv,accumarray(idx(:),1)];
iterH_in.theta_h      = sortrows(theta_df);

iterH_in.firm_h_yr_sales_lag = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),4);
% firm_yr_sales_lag will contain: [firmID,sales,#shipments,firm age]
iterH_in.mat_h_yr_sales_lag = double.empty(0,4);
iterH_in.mat_h_yr_sales = double.empty(0,4);

% initialize keep_cli for first period
iterH_in.keep_cli = ones(1,size(mm.Z,1)); % applies to clients existing in period 1
iterH_in.keep_cli(1:5) =  zeros(1,5); % implying worst 5 client types from period 1 are dropped
iterH_in.year = 1;
iterH_in.N_match = 0;
iterH_in.season = 1;

iter_out.transH = cell(mm.N_pt,5);
% to hold (1) firm_ID, (2) cur_cli_cnt, (3) cum_succ, (4) age, (5) new_firm 

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

iter_out.abort_flag_h = 0;
end
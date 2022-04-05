function [iter_in, iter_out] = simulateForeignInnerInitialize(mm, pt_ndx)

iter_in = struct;
iter_out = struct;

iter_in.pt_ndx       = pt_ndx;
iter_in.seas_tran = cell(1,mm.pd_per_yr);  % cells will hold one year's worth of season- and match-specific outcomes for all firms w/in type
iter_in.seas_Zcut = zeros(1,mm.pd_per_yr); % elements will hold season-specifics Z cut-offs for endog. drops
iter_in.cur_cli_cnt  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % clients active in the current period
iter_in.add_cli_cnt  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % gross additions to client count
iter_in.cum_meets    = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % cumulative number of meetings
iter_in.cum_succ     = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % cumulative number of successes
iter_in.new_firm     = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods);   % will mark new firms that haven't made a match
iter_in.exit_firm    = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods);   % will mark last period of exiting firm
iter_in.exog_deaths  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1); % number of exogenous match deaths
iter_in.micro_state  = ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.periods,1);  % scalar indices for #success/#meetings
iter_in.lag_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % breaks down lagged clients by z state
iter_in.new_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % breaks down new client counts by z state
iter_in.die_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % breaks down client death counts by z state
iter_in.surviv_zst   = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % breaks down surviving client counts by z state
iter_in.trans_zst    = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));  % counts survival types by firm after z innovations
iter_in.flrlag       = ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1);     % initializing vector for age debugging
iter_in.cumage       = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1);    % initializing vector for age debugging
iter_in.mkt_exit                  = zeros(1,3);
iter_in.mat_yr_sales_adj          = zeros(0,9);
iter_in.trans_count    = zeros(size(mm.Z,1)+1,size(mm.Z,1)+1,mm.sim_firm_num_by_prod_succ_type(pt_ndx)); % counts transitions across buyer types,
% for each seller type. New buyer types are considered type 0 at beginning of period, hence the +1.
% Exiting firms are considered to move to type 0 at the the end of the period.
% Dimensions: (1) initial z-state (2) new z-state (3) firm index, given type

iter_in.keep_cli      = ones(1,size(mm.Z,1)); % applies to clients existing in period 1
iter_in.keep_cli(1:5) =  zeros(1,5); % implying worst 5 client types from period 1 are dropped.

% create first observation on firm-year level aggregates (will concatenate below)
iter_out.match_count      = zeros(mm.max_match,1);
iter_out.mat_yr_sales     = zeros(0,9);
iter_out.mat_yr_sales_adj = zeros(0,9);
iter_out.firm_f_yr_sales  = zeros(0,6);
iter_out.time_gaps        = zeros(0,7);
iter_out.mat_ar1_x        = zeros(0,4);
iter_out.mat_ar1_y        = zeros(0,1);

% match level moment aggregators
iter_out.moms_xx   = zeros(4,4);
iter_out.moms_xy   = zeros(4,1);
iter_out.ysum      = 0;
iter_out.nobs      = 0;
iter_out.ship_obs  = 0;
iter_out.ln_ships = 0;

% firm level moment aggregators
iter_out.fmoms_xx = zeros(4,4);
iter_out.fmoms_xy = zeros(4,1);
iter_out.fysum    = 0;
iter_out.fnobs    = 0;

iter_out.exit_xx       = zeros(6,6);
iter_out.exit_xy       = zeros(1,6);
iter_out.sum_succ_rate = 0;
iter_out.exit_obs      = 0;
iter_out.sum_exits     = 0;

iter_out.mat_exit_moms_xx = zeros(5,5);
iter_out.mat_exit_moms_xy = zeros(5,1);
iter_out.mat_obs          = 0;
iter_out.nmat_exit        = 0;
iter_out.mat_exit_x       = zeros(0,5);
iter_out.mat_exit_y       = zeros(0,1);
iter_out.mat_matur        = zeros(0,7);
end
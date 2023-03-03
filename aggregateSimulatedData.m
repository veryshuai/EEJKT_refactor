function sim_cum = aggregateSimulatedData(sim_out,mm)

sim_cum = struct;

sim_cum.agg_mat_yr_sales    = zeros(0,9); % obs. for analysis of match dynamics
sim_cum.agg_dud_matches     = zeros(0,9); % obs. on dud matches
sim_cum.agg_mat_ar1_x       = zeros(0,4); % for match ar1 regression residuals
sim_cum.agg_mat_ar1_y       = zeros(0,1); % for match ar1 regression residuals
sim_cum.agg_mat_exit_x      = zeros(0,5); % for match exit regression residuals
sim_cum.agg_mat_exit_y      = zeros(0,1); % for match exit regression residuals
sim_cum.agg_x_hf            = zeros(0,2); % for home-foreign sales regression residuals
sim_cum.agg_y_hf            = zeros(0,1); % for home-foreign sales regression residuals
sim_cum.agg_x_fsales_h      = zeros(0,2); % for firm-level home sales AR1 residuals
sim_cum.agg_y_fsales_h      = zeros(0,1); % for firm-level home sales AR1 residuals
sim_cum.agg_time_gaps       = zeros(0,7); % for match hazard analysis
sim_cum.agg_moms_xx   = zeros(4,4);       % match autoregression
sim_cum.agg_moms_xy   = zeros(4,1);       % match autoregression
sim_cum.agg_ysum      = 0;                % match autoregression
sim_cum.agg_nobs      = 0;                % match autoregression
sim_cum.agg_ship_obs  = 0;                % number of shipments
sim_cum.agg_fmoms_h_xx = zeros(2,2); % home market firm sales autoreg.
sim_cum.agg_fmoms_h_xy = zeros(2,1); % home market firm sales autoreg.
sim_cum.agg_fysum_h  = 0; % cross-firm sum of domestic sales, dom. sales autoreg.
sim_cum.agg_fnobs_h  = 0;      % # firms in domestic sales autoreg.
sim_cum.agg_sum_succ_rate = 0; % firm-level ratios, cum. succ to cum. meetings 
sim_cum.agg_hfmoms_xx = zeros(2,2); % for firm-level home-foreign regression
sim_cum.agg_hfmoms_xy = zeros(2,1); % for firm-level home-foreign regression
sim_cum.agg_hfysum    = 0;          % for firm-level home-foreign regression
sim_cum.agg_hf_nobs    = 0;         % for firm-level home-foreign regression
sim_cum.agg_nfirm     = 0;          % for # active firms, pooling years
sim_cum.agg_nexptr    = 0;          % for # exporting firms, pooling years
sim_cum.agg_expt_rate = zeros(0,1); % for share of output exported among exporters
sim_cum.agg_mat_exit_moms_xx  = zeros(5,5); % match death regression
sim_cum.agg_mat_exit_moms_xy  = zeros(5,1); % match death regression
sim_cum.agg_mat_obs       = 0;  % # matches used, match exit regression
sim_cum.agg_nmat_exit     = 0;  % # match exits, match exit regression
sim_cum.agg_match_hist  = zeros(1,mm.max_match); % # matches for each firm, excl. duds
sim_cum.agg_match_histD = zeros(1,mm.max_match); % # matches for each firm, incl. duds

sim_cum.duds       = 0;
sim_cum.agg_ln_ships = 0;

%% Cumulate over firm types

for pt_ndx = 1:1:mm.N_pt
% for pt_ndx = 113
% mm.sim_firm_num_by_prod_succ_type(pt_ndx)
    if mm.sim_firm_num_by_prod_succ_type(pt_ndx) > 0
        sim_cum.agg_time_gaps    = [sim_cum.agg_time_gaps;sim_out{pt_ndx}.time_gaps];
        sim_cum.agg_mat_yr_sales = [sim_cum.agg_mat_yr_sales;sim_out{pt_ndx}.mat_yr_sales];
        sim_cum.agg_dud_matches  = [sim_cum.agg_dud_matches;sim_out{pt_ndx}.dud_matches];         
        sim_cum.agg_x_hf       = [sim_cum.agg_x_hf;sim_out{pt_ndx}.x_hf];
        sim_cum.agg_y_hf       = [sim_cum.agg_y_hf;sim_out{pt_ndx}.y_hf];
        sim_cum.agg_x_fsales_h = [sim_cum.agg_x_fsales_h;sim_out{pt_ndx}.x_fsales_h];
        sim_cum.agg_y_fsales_h = [sim_cum.agg_y_fsales_h;sim_out{pt_ndx}.y_fsales_h];
        sim_cum.agg_moms_xx    = sim_cum.agg_moms_xx + squeeze(sim_out{pt_ndx}.moms_xx); 
        sim_cum.agg_moms_xy    = sim_cum.agg_moms_xy + sim_out{pt_ndx}.moms_xy; %check direction
        sim_cum.agg_ysum       = sim_cum.agg_ysum    + sim_out{pt_ndx}.ysum;
        sim_cum.agg_nobs       = sim_cum.agg_nobs    + sim_out{pt_ndx}.nobs;
        sim_cum.agg_mat_ar1_x  = [sim_cum.agg_mat_ar1_x;sim_out{pt_ndx}.mat_ar1_x];
        sim_cum.agg_mat_ar1_y  = [sim_cum.agg_mat_ar1_y;sim_out{pt_ndx}.mat_ar1_y];
        sim_cum.agg_fmoms_h_xx = sim_cum.agg_fmoms_h_xx + squeeze(sim_out{pt_ndx}.fmoms_h_xx); 
        sim_cum.agg_fmoms_h_xy = sim_cum.agg_fmoms_h_xy + squeeze(sim_out{pt_ndx}.fmoms_h_xy); %check direction
        sim_cum.agg_fysum_h    = sim_cum.agg_fysum_h + sim_out{pt_ndx}.fysum_h;
        sim_cum.agg_fnobs_h    = sim_cum.agg_fnobs_h + sim_out{pt_ndx}.fnobs_h;
        sim_cum.agg_hfmoms_xx  = sim_cum.agg_hfmoms_xx + squeeze(sim_out{pt_ndx}.hfmoms_xx);

        sim_cum.agg_hfmoms_xy  = sim_cum.agg_hfmoms_xy + squeeze(sim_out{pt_ndx}.hfmoms_xy); %check direction       
        sim_cum.agg_hfysum     = sim_cum.agg_hfysum    + sim_out{pt_ndx}.hfysum;
        sim_cum.agg_hf_nobs    = sim_cum.agg_hf_nobs   + sim_out{pt_ndx}.hf_nobs;              
                    
        sim_cum.agg_nfirm      = sim_cum.agg_nfirm + sim_out{pt_ndx}.nfirm;   % number of firms
        sim_cum.agg_nexptr     = sim_cum.agg_nexptr + sim_out{pt_ndx}.nexptr; % number of exporters
        sim_cum.agg_expt_rate  = [sim_cum.agg_expt_rate;sim_out{pt_ndx}.expt_rate]; %vector of shares of output exported, exporters only
        
        sim_cum.agg_sum_succ_rate    = sim_cum.agg_sum_succ_rate + sim_out{pt_ndx}.sum_succ_rate;
        sim_cum.agg_mat_exit_moms_xx = sim_cum.agg_mat_exit_moms_xx + squeeze(sim_out{pt_ndx}.mat_exit_moms_xx);
        sim_cum.agg_mat_exit_moms_xy = sim_cum.agg_mat_exit_moms_xy + squeeze(sim_out{pt_ndx}.mat_exit_moms_xy);
        
        sim_cum.agg_mat_obs      = sim_cum.agg_mat_obs + sim_out{pt_ndx}.mat_obs;
        sim_cum.agg_nmat_exit    = sim_cum.agg_nmat_exit + sim_out{pt_ndx}.nmat_exit;
        sim_cum.agg_mat_exit_x   = [sim_cum.agg_mat_exit_x;sim_out{pt_ndx}.mat_exit_x];
        sim_cum.agg_mat_exit_y   = [sim_cum.agg_mat_exit_y;sim_out{pt_ndx}.mat_exit_y];
        sim_cum.agg_ship_obs     = sim_cum.agg_ship_obs    + sim_out{pt_ndx}.ship_obs;
        sim_cum.agg_ln_ships     = sim_cum.agg_ln_ships    + sim_out{pt_ndx}.ln_ships;
        sim_cum.agg_match_hist  = sim_cum.agg_match_hist + sim_out{pt_ndx}.match_hist; 
        sim_cum.agg_match_histD = sim_cum.agg_match_histD + sim_out{pt_ndx}.match_histD; 
        sim_cum.duds             = sim_cum.duds + sim_out{pt_ndx}.duds;
    end

end       
end
function sim_cum = aggregateSimulatedData(sim_out,mm)

sim_cum = struct;

sim_cum.agg_mat_yr_sales    = zeros(0,9); % for analysis of match dynamics
sim_cum.agg_dud_matches     = zeros(0,9);
sim_cum.agg_mat_ar1_x       = zeros(0,4); % for match ar1 regression residuals
sim_cum.agg_mat_ar1_y       = zeros(0,1); % for match ar1 regression residuals
sim_cum.agg_mat_exit_x      = zeros(0,5); % for match exit regression residuals
sim_cum.agg_mat_exit_y      = zeros(0,1); % for match exit regression residuals
% sim_cum.agg_mat_matur       = zeros(0,8); % for match maturation analysis
sim_cum.agg_x_hf            = zeros(0,2); % for home-foreign sales regression residuals
sim_cum.agg_y_hf            = zeros(0,1); % for home-foreign sales regression residuals
sim_cum.agg_x_fsales_h      = zeros(0,2); % for home sales AR1 residuals
sim_cum.agg_y_fsales_h      = zeros(0,1); % for home sales AR1 residuals
sim_cum.agg_time_gaps       = zeros(0,7); % for match hazard analysis
sim_cum.agg_moms_xx   = zeros(4,4);
sim_cum.agg_moms_xy   = zeros(4,1);
sim_cum.agg_ysum      = 0;
sim_cum.agg_nobs      = 0;
sim_cum.agg_ship_obs  = 0;
sim_cum.agg_fmoms_xx = zeros(4,4);
sim_cum.agg_fmoms_xy = zeros(4,1);
sim_cum.agg_fmoms_h_xx = zeros(2,2);
sim_cum.agg_fmoms_h_xy = zeros(2,1);
sim_cum.agg_fysum    = 0;
sim_cum.agg_fnobs    = 0;
sim_cum.agg_fysum_h  = 0;
sim_cum.agg_fnobs_h  = 0;
sim_cum.agg_exit_xx = zeros(6,6);
sim_cum.agg_exit_xy = zeros(6,1);
sim_cum.agg_sum_succ_rate = 0;
sim_cum.agg_exit_obs = 0;
sim_cum.agg_sum_exits = 0;
sim_cum.agg_hfmoms_xx = zeros(2,2);
sim_cum.agg_hfmoms_xy = zeros(2,1);
sim_cum.agg_hfysum    = 0;
sim_cum.agg_hf_nobs    = 0;
sim_cum.agg_nfirm     = 0;
sim_cum.agg_nexptr    = 0;
sim_cum.agg_expt_rate = zeros(0,1);
sim_cum.agg_mat_exit_moms_xx  = zeros(5,5);
sim_cum.agg_mat_exit_moms_xy  = zeros(5,1);
sim_cum.agg_mat_obs       = 0;
sim_cum.agg_nmat_exit     = 0;
sim_cum.agg_match_count  = zeros(1,mm.max_match);
sim_cum.agg_match_countD = zeros(1,mm.max_match);

sim_cum.singletons = 0;
sim_cum.duds       = 0;
sim_cum.agg_ln_ships = 0;

%% Cumulate over firm types

for pt_ndx = 1:1:mm.N_pt
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
        sim_cum.agg_fmoms_xx   = sim_cum.agg_fmoms_xx + squeeze(sim_out{pt_ndx}.fmoms_xx); 
        sim_cum.agg_fmoms_xy   = sim_cum.agg_fmoms_xy + squeeze(sim_out{pt_ndx}.fmoms_xy);  %check direction
        sim_cum.agg_fysum      = sim_cum.agg_fysum + sim_out{pt_ndx}.fysum;
        sim_cum.agg_fnobs      = sim_cum.agg_fnobs + sim_out{pt_ndx}.fnobs;
        sim_cum.agg_fmoms_h_xx = sim_cum.agg_fmoms_h_xx + squeeze(sim_out{pt_ndx}.fmoms_h_xx); 
        sim_cum.agg_fmoms_h_xy = sim_cum.agg_fmoms_h_xy + squeeze(sim_out{pt_ndx}.fmoms_h_xy); %check direction
        sim_cum.agg_fysum_h    = sim_cum.agg_fysum_h + sim_out{pt_ndx}.fysum_h;
        sim_cum.agg_fnobs_h    = sim_cum.agg_fnobs_h + sim_out{pt_ndx}.fnobs_h;
        sim_cum.agg_hfmoms_xx  = sim_cum.agg_hfmoms_xx + squeeze(sim_out{pt_ndx}.hfmoms_xx);
        sim_cum.agg_hfmoms_xy  = sim_cum.agg_hfmoms_xy + squeeze(sim_out{pt_ndx}.hfmoms_xy); %check direction
        sim_cum.agg_hfysum     = sim_cum.agg_hfysum    + sim_out{pt_ndx}.hfysum;
        sim_cum.agg_hf_nobs    = sim_cum.agg_hf_nobs   + sim_out{pt_ndx}.hf_nobs;
        sim_cum.agg_nfirm      = sim_cum.agg_nfirm + sim_out{pt_ndx}.nfirm;
        sim_cum.agg_nexptr     = sim_cum.agg_nexptr + sim_out{pt_ndx}.nexptr;
        sim_cum.agg_expt_rate  = [sim_cum.agg_expt_rate;sim_out{pt_ndx}.expt_rate];
%       sim_cum.agg_exit_xx    = sim_cum.agg_exit_xx + squeeze(sim_out{pt_ndx}.exit_xx);
        sim_cum.agg_exit_xy    = sim_cum.agg_exit_xy + squeeze(sim_out{pt_ndx}.exit_xy)'; %chekc direction
        
        sim_cum.agg_sum_succ_rate    = sim_cum.agg_sum_succ_rate + sim_out{pt_ndx}.sum_succ_rate;
        sim_cum.agg_exit_obs         = sim_cum.agg_exit_obs + sim_out{pt_ndx}.exit_obs;
        sim_cum.agg_sum_exits        = sim_cum.agg_sum_exits + sim_out{pt_ndx}.sum_exits;
        sim_cum.agg_mat_exit_moms_xx = sim_cum.agg_mat_exit_moms_xx + squeeze(sim_out{pt_ndx}.mat_exit_moms_xx);
        sim_cum.agg_mat_exit_moms_xy = sim_cum.agg_mat_exit_moms_xy + squeeze(sim_out{pt_ndx}.mat_exit_moms_xy);
        
        sim_cum.agg_mat_obs      = sim_cum.agg_mat_obs + sim_out{pt_ndx}.mat_obs;
        sim_cum.agg_nmat_exit    = sim_cum.agg_nmat_exit + sim_out{pt_ndx}.nmat_exit;
        sim_cum.agg_mat_exit_x   = [sim_cum.agg_mat_exit_x;sim_out{pt_ndx}.mat_exit_x];
        sim_cum.agg_mat_exit_y   = [sim_cum.agg_mat_exit_y;sim_out{pt_ndx}.mat_exit_y];
        sim_cum.agg_ship_obs     = sim_cum.agg_ship_obs    + sim_out{pt_ndx}.ship_obs;
        sim_cum.agg_ln_ships     = sim_cum.agg_ln_ships    + sim_out{pt_ndx}.ln_ships;
        sim_cum.agg_match_count  = sim_cum.agg_match_count + sim_out{pt_ndx}.match_count; 
        sim_cum.agg_match_countD = sim_cum.agg_match_countD + sim_out{pt_ndx}.match_countD; 
        sim_cum.singletons       = sim_cum.singletons + sim_out{pt_ndx}.singletons;
        sim_cum.duds             = sim_cum.duds + sim_out{pt_ndx}.duds;
    end

end       
end
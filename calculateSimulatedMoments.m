function simMoms = calculateSimulatedMoments(sim_cum)

%% Construct simulated statistics
simMoms = struct; %container for all simulated moments

simMoms.agg_nexptr = sim_cum.agg_nexptr;
simMoms.agg_nfirm = sim_cum.agg_nfirm;

% match-level autoregression
if rank(sim_cum.agg_moms_xx) == size(sim_cum.agg_moms_xx,2)
    inv_agg_moms_xx = inv(sim_cum.agg_moms_xx);
else
    rank_xx = rank(sim_cum.agg_moms_xx);
    fprintf('\r Warning: singular matrix for match-level AR1. Rank: %.1f\n', rank_xx);
    inv_agg_moms_xx = [inv(sim_cum.agg_moms_xx(1:2,1:2)), zeros(2,2); zeros(2,4)];
end
simMoms.beta_match = inv_agg_moms_xx*sim_cum.agg_moms_xy;
simMoms.ybar_match = sim_cum.agg_ysum/sim_cum.agg_nobs;
simMoms.mse_match_ar1 = (sim_cum.agg_mat_ar1_y - sim_cum.agg_mat_ar1_x*simMoms.beta_match)'*...
    (sim_cum.agg_mat_ar1_y - sim_cum.agg_mat_ar1_x*simMoms.beta_match)/size(sim_cum.agg_mat_ar1_x,1);


% firm-level autoregression
simMoms.beta_fsales_h = inv(sim_cum.agg_fmoms_h_xx)*sim_cum.agg_fmoms_h_xy;
simMoms.mse_h         = (sim_cum.agg_y_fsales_h - sim_cum.agg_x_fsales_h*simMoms.beta_fsales_h)'*...
    (sim_cum.agg_y_fsales_h - sim_cum.agg_x_fsales_h*simMoms.beta_fsales_h)/size(sim_cum.agg_y_fsales_h,1);
simMoms.ybar_fsales_h = sim_cum.agg_fysum_h/sim_cum.agg_fnobs_h ;

% home-foreign regression and export rates
simMoms.beta_hfsales = inv(sim_cum.agg_hfmoms_xx)*sim_cum.agg_hfmoms_xy;
simMoms.mse_hf       = (sim_cum.agg_y_hf - sim_cum.agg_x_hf*simMoms.beta_hfsales)'*...
    (sim_cum.agg_y_hf - sim_cum.agg_x_hf*simMoms.beta_hfsales)/size(sim_cum.agg_x_hf,1);
simMoms.ybar_hfsales  = sim_cum.agg_hfysum/sim_cum.agg_hf_nobs;
simMoms.avg_expt_rate = mean(sim_cum.agg_expt_rate);
simMoms.share_exptr   = sim_cum.agg_nexptr/sim_cum.agg_nfirm;

% market exit regression
if rank(sim_cum.agg_exit_xx) == size(sim_cum.agg_exit_xx,2)
    inv_agg_exit_xx = inv(sim_cum.agg_exit_xx);
else
    rank_xx = rank(sim_cum.agg_exit_xx);
    fprintf('\r Warning: singular matrix for market exit. Rank: %.1f\n', rank_xx);
    inv_agg_exit_xx = [inv(sim_cum.agg_exit_xx(1:3,1:3)), zeros(3,3); zeros(3,6)];
end
simMoms.beta_mkt_exit   = inv_agg_exit_xx*sim_cum.agg_exit_xy;
simMoms.mkt_exit_rate   = sim_cum.agg_sum_exits/sim_cum.agg_exit_obs;

% match exit regression
simMoms.beta_match_exit = inv(sim_cum.agg_mat_exit_moms_xx)*sim_cum.agg_mat_exit_moms_xy; 
simMoms.match_exit_rate = sim_cum.agg_nmat_exit/sim_cum.agg_mat_obs;

% average log #shipments
simMoms.avg_ln_ships = sim_cum.agg_ln_ships/sim_cum.agg_ship_obs;

% create variables for analysis of degree distribution
simMoms.ff_sim_max      = find(cumsum(sim_cum.agg_match_count)./sum(sim_cum.agg_match_count)<1);
log_compCDF     = log(1 - cumsum(sim_cum.agg_match_count(simMoms.ff_sim_max))./sum(sim_cum.agg_match_count));
log_matches     = log(1:1:size(simMoms.ff_sim_max,1))';
xmat            = [ones(size(simMoms.ff_sim_max)),log_matches,log_matches.^2];
% quadratic regression approximating degree distribution
simMoms.b_degree        = regress(log_compCDF,xmat);
% linear regression approximating degree distribution

% plot histogram of frequencies for meeting hazards
sim_cum.agg_time_gaps = sim_cum.agg_time_gaps(2:size(sim_cum.agg_time_gaps,1),:);

% create variables for hazard regressions
ln_haz = log(1./sim_cum.agg_time_gaps(:,3));
ln_csucc = log(1+sim_cum.agg_time_gaps(:,7));
ln_meet = log(sim_cum.agg_time_gaps(:,6));
const = ones(size(ln_haz,1),1);
ln_succ_rate = log(1+(sim_cum.agg_time_gaps(:,7)./sim_cum.agg_time_gaps(:,6)));

% success rate regression
succ_rate = sim_cum.agg_time_gaps(:,7)./sim_cum.agg_time_gaps(:,6);
[simMoms.b_succ_rate,~,uu] = regress(succ_rate,[const, ln_meet]);
usq_succ = uu.^2;
simMoms.b_usq_succ = regress(usq_succ,[const, ln_meet]);

% translog meeting hazard regression
X_haz = [const, ln_csucc, ln_csucc.^2, ln_succ_rate, ln_succ_rate.^2, ln_succ_rate.*ln_csucc];
simMoms.b_haz = regress(ln_haz,X_haz);

% means of log hazard rate, success rate, and squared residuals from succ rate
means_vec = mean([ln_haz,succ_rate,usq_succ]);
simMoms.mean_ln_haz    = means_vec(1);
simMoms.mean_succ_rate = means_vec(2);
simMoms.mean_usq_succ  = means_vec(3);

end
%function moms = simulateMomentsMain(mm,policy)

rng(80085,'twister');

lambda_h = policy.lambda_h;
lambda_f = policy.lambda_f;
c_val_h_orig = policy.c_val_h;
c_val_f_orig = policy.c_val_f;
value_f = policy.value_f;
value_h = policy.value_h;
nn_h = 70; %apparently this was the total possible number of clients at home, set deep in the code
typemat = policy.firm_type_prod_succ_macro;
pmat_cum_f = policy.pmat_cum_f;
pmat_cum_h = policy.pmat_cum_h;
pmat_cum_msf = policy.pmat_cum_msf;
pmat_cum_msh = policy.pmat_cum_msh;
pmat_cum_z = policy.pmat_cum_z;

[macro_state_f, macro_state_h] = simulateMacroTrajectories(mm, policy);

%% Create objects that will store firm-type specific simulated data

sim_out_f = cell(mm.N_pt,1);
sim_out_h = cell(mm.N_pt,1);
sim_out_hf = cell(mm.N_pt,1);
sim_out = cell(mm.N_pt,1);


%% Load storage objects by looping over firm types

seeds = randi(1e6,size(mm.Phi,1),2);

%parfor pt_ndx = 1:1:mm.N_pt
for pt_ndx = 1:1:mm.N_pt % use this for loop for debugging only
    
    rng(seeds(mm.pt_type(pt_ndx,1),1),'twister');
    seed_crand(seeds(mm.pt_type(pt_ndx,1),2));

    if mm.sim_firm_num_by_prod_succ_type(pt_ndx)>0

        sim_out_f{pt_ndx} = matchdat_gen_f(pt_ndx,macro_state_f, mm, policy);

        sim_out_h{pt_ndx} = matchdat_gen_h(pt_ndx,macro_state_h, mm, policy);

        sim_out_hf{pt_ndx} = splice_hf(sim_out_h{pt_ndx},sim_out_f{pt_ndx},policy,mm);

        sim_out{pt_ndx} = cell2struct([struct2cell(sim_out_h{pt_ndx});struct2cell(sim_out_f{pt_ndx});struct2cell(sim_out_hf{pt_ndx})],...
            [fieldnames(sim_out_h{pt_ndx});fieldnames(sim_out_f{pt_ndx});fieldnames(sim_out_hf{pt_ndx})]);

    end

end

%% Initialize objects that will hold cumulated values

sim_cum = aggregateSimulatedData(sim_out,mm);

%% Construct simulated statistics

% save ('match_count.mat','s_match_count','pi_tilda_f_new','pt_type','NN')
% active_client_value

simMoms = struct; %container for all simulated moments

% Some numbers from above
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
beta_fsales = inv(sim_cum.agg_fmoms_xx)*sim_cum.agg_fmoms_xy;
ybar_fsales = sim_cum.agg_fysum/sim_cum.agg_fnobs;

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
match_succ_rate = sim_cum.agg_sum_succ_rate/sim_cum.agg_exit_obs;

% match exit regression
simMoms.beta_match_exit = inv(sim_cum.agg_mat_exit_moms_xx)*sim_cum.agg_mat_exit_moms_xy; 
simMoms.match_exit_rate = sim_cum.agg_nmat_exit/sim_cum.agg_mat_obs;
mse_match_exit  = (sim_cum.agg_mat_exit_y - sim_cum.agg_mat_exit_x*simMoms.beta_match_exit)'*...
    (sim_cum.agg_mat_exit_y - sim_cum.agg_mat_exit_x*simMoms.beta_match_exit)/size(sim_cum.agg_mat_exit_y,1);


% average log #shipments
simMoms.avg_ln_ships = sim_cum.agg_ln_ships/sim_cum.agg_ship_obs;

%       % create variables for analysis of degree distribution
%
simMoms.ff_sim_max      = find(cumsum(sim_cum.agg_match_count)./sum(sim_cum.agg_match_count)<1);
log_compCDF     = log(1 - cumsum(sim_cum.agg_match_count(simMoms.ff_sim_max))./sum(sim_cum.agg_match_count));
log_matches     = log(1:1:size(simMoms.ff_sim_max,1))';
xmat            = [ones(size(simMoms.ff_sim_max)),log_matches,log_matches.^2];

% quadratic regression approximating degree distribution
simMoms.b_degree        = regress(log_compCDF,xmat);
% linear regression approximating degree distribution
xmat_linear     = [ones(size(simMoms.ff_sim_max)),log_matches];
b_degree_linear = regress(log_compCDF,xmat_linear);
% nonparametric plot of degree distribution
%       scatter(log_matches,log_compCDF)


% plot histogram of frequencies for meeting hazards
sim_cum.agg_time_gaps = sim_cum.agg_time_gaps(2:size(sim_cum.agg_time_gaps,1),:);
%         histogram(sim_cum.agg_time_gaps(:,3))

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
% regression explaining squared residuals of meeting hazard equation
usq_haz = (ln_haz - X_haz*simMoms.b_haz).^2;
b_haz_usq = regress(usq_haz,[const ln_meet]);
% means of log hazard rate, success rate, and squared residuals from succ rate and hazard regression
means_vec = mean([ln_haz,succ_rate,usq_succ,usq_haz]);
simMoms.mean_ln_haz    = means_vec(1);
simMoms.mean_succ_rate = means_vec(2);
simMoms.mean_usq_succ  = means_vec(3);
mean_usq_haz   = means_vec(4);

















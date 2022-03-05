%function moms = simulateMomentsMain(mm,policy)

rng(80085,'twister');

%seed = 1;
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

NN            = zeros(mm.N_pt,1);

%% Load storage objects by looping over firm types

N_theta1 = size(mm.theta1,2);    % number of thetas in home market
th1_cdf  = betacdf(mm.theta1,mm.ah,mm.bh); % cdf for home theta draws
th1_cdf(N_theta1) = 1; % to deal with rounding problem

seeds = randi(1e6,size(mm.Phi,1),2);
pt_type = [kron((1:size(mm.Phi,1))',ones(size(mm.theta2,2),1)),kron(ones(size(mm.Phi,1),1),(1:size(mm.theta2,2))')];

parfor pt_ndx = 1:1:mm.N_pt
%for pt_ndx = 1:1:mm.N_pt % use this for loop for debugging only

    prod_ndx  = pt_type(pt_ndx,1);
    theta_ndx = pt_type(pt_ndx,2);

    rng(seeds(prod_ndx,1),'twister');
    seed_crand(seeds(prod_ndx,2));

    prod_lvl  = mm.Phi(prod_ndx);
    succ_prob = mm.theta2(theta_ndx);

    % number of firms to simulate of this particular type
    N_firms = round(mm.erg_pp(prod_ndx)*policy.th2_pdf(theta_ndx)*mm.S);
    NN(pt_ndx) = N_firms;

    if N_firms>0

        sim_out_f{pt_ndx} = matchdat_gen_f(N_firms, macro_state_f, theta_ndx, prod_ndx, mm, policy);

        sim_out_h{pt_ndx} = matchdat_gen_h(N_firms,size(mm.Z,1),size(mm.Phi,1),nn_h,mm.pd_per_yr,1/mm.pd_per_yr,mm.periods,mm.max_ships,typemat,macro_state_h,...
            theta_ndx,prod_ndx,mm,pmat_cum_h,c_val_h_orig,cumsum(mm.erg_pz),mm.poisCDF_shipments,...
            1-exp(-mm.delta),1-exp(-mm.firm_death_haz),mm.L_b,pmat_cum_z,N_theta1,th1_cdf);

        sim_out_hf{pt_ndx} = splice_hf(sim_out_h{pt_ndx},sim_out_f{pt_ndx},policy,mm);

        sim_out{pt_ndx} = cell2struct([struct2cell(sim_out_h{pt_ndx});struct2cell(sim_out_f{pt_ndx});struct2cell(sim_out_hf{pt_ndx})],...
            [fieldnames(sim_out_h{pt_ndx});fieldnames(sim_out_f{pt_ndx});fieldnames(sim_out_hf{pt_ndx})]);

    end

end

rand

%% Initialize objects that will hold cumulated values

sim_cum = struct;

sim_cum.agg_mat_yr_sales    = zeros(0,9); % for analysis of match dynamics
sim_cum.agg_mat_yr_sales_adj= zeros(0,9); % for analysis of match exit
sim_cum.agg_mat_ar1_x       = zeros(0,4); % for match ar1 regression residuals
sim_cum.agg_mat_ar1_y       = zeros(0,1); % for match ar1 regression residuals
sim_cum.agg_mat_exit_x      = zeros(0,5); % for match exit regression residuals
sim_cum.agg_mat_exit_y      = zeros(0,1); % for match exit regression residuals
sim_cum.agg_mat_matur       = zeros(0,8); % for match maturation analysis
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
sim_cum.agg_match_count  = zeros(mm.max_match,1);
sim_cum.singletons = 0;
sim_cum.agg_ln_ships = 0;

%% Cumulate over firm types

for pt_ndx = 1:1:mm.N_pt

    if NN(pt_ndx) > 0
        sim_cum.agg_time_gaps = [sim_cum.agg_time_gaps;sim_out{pt_ndx}.time_gaps];
        sim_cum.agg_mat_yr_sales  = [sim_cum.agg_mat_yr_sales;sim_out{pt_ndx}.mat_yr_sales];
        sim_cum.agg_mat_yr_sales_adj  = [sim_cum.agg_mat_yr_sales_adj;sim_out{pt_ndx}.mat_yr_sales_adj];
        sim_cum.agg_x_hf = [sim_cum.agg_x_hf;sim_out{pt_ndx}.x_hf];
        sim_cum.agg_y_hf = [sim_cum.agg_y_hf;sim_out{pt_ndx}.y_hf];
        sim_cum.agg_x_fsales_h = [sim_cum.agg_x_fsales_h;sim_out{pt_ndx}.x_fsales_h];
        sim_cum.agg_y_fsales_h = [sim_cum.agg_y_fsales_h;sim_out{pt_ndx}.y_fsales_h];
        sim_cum.agg_moms_xx = sim_cum.agg_moms_xx + squeeze(sim_out{pt_ndx}.moms_xx); 
        sim_cum.agg_moms_xy = sim_cum.agg_moms_xy + sim_out{pt_ndx}.moms_xy; %check direction
        sim_cum.agg_ysum    = sim_cum.agg_ysum    + sim_out{pt_ndx}.ysum;
        sim_cum.agg_nobs    = sim_cum.agg_nobs    + sim_out{pt_ndx}.nobs;
        sim_cum.agg_mat_ar1_x = [sim_cum.agg_mat_ar1_x;sim_out{pt_ndx}.mat_ar1_x];
        sim_cum.agg_mat_ar1_y = [sim_cum.agg_mat_ar1_y;sim_out{pt_ndx}.mat_ar1_y];
        sim_cum.agg_fmoms_xx = sim_cum.agg_fmoms_xx + squeeze(sim_out{pt_ndx}.fmoms_xx); 
        sim_cum.agg_fmoms_xy = sim_cum.agg_fmoms_xy + squeeze(sim_out{pt_ndx}.fmoms_xy);  %check direction
        sim_cum.agg_fysum    = sim_cum.agg_fysum + sim_out{pt_ndx}.fysum;
        sim_cum.agg_fnobs    = sim_cum.agg_fnobs + sim_out{pt_ndx}.fnobs;
        sim_cum.agg_fmoms_h_xx = sim_cum.agg_fmoms_h_xx + squeeze(sim_out{pt_ndx}.fmoms_h_xx); 
        sim_cum.agg_fmoms_h_xy = sim_cum.agg_fmoms_h_xy + squeeze(sim_out{pt_ndx}.fmoms_h_xy); %check direction
        sim_cum.agg_fysum_h    = sim_cum.agg_fysum_h + sim_out{pt_ndx}.fysum_h;
        sim_cum.agg_fnobs_h    = sim_cum.agg_fnobs_h + sim_out{pt_ndx}.fnobs_h;
        sim_cum.agg_hfmoms_xx = sim_cum.agg_hfmoms_xx + squeeze(sim_out{pt_ndx}.hfmoms_xx);
        sim_cum.agg_hfmoms_xy = sim_cum.agg_hfmoms_xy + squeeze(sim_out{pt_ndx}.hfmoms_xy); %check direction
        sim_cum.agg_hfysum    = sim_cum.agg_hfysum    + sim_out{pt_ndx}.hfysum;
        sim_cum.agg_hf_nobs   = sim_cum.agg_hf_nobs   + sim_out{pt_ndx}.hf_nobs;
        sim_cum.agg_nfirm     = sim_cum.agg_nfirm + sim_out{pt_ndx}.nfirm;
        sim_cum.agg_nexptr    = sim_cum.agg_nexptr + sim_out{pt_ndx}.nexptr;
        sim_cum.agg_expt_rate = [sim_cum.agg_expt_rate;sim_out{pt_ndx}.expt_rate];
        sim_cum.agg_exit_xx = sim_cum.agg_exit_xx + squeeze(sim_out{pt_ndx}.exit_xx);
        sim_cum.agg_exit_xy = sim_cum.agg_exit_xy + squeeze(sim_out{pt_ndx}.exit_xy)'; %chekc direction
        sim_cum.agg_sum_succ_rate = sim_cum.agg_sum_succ_rate + sim_out{pt_ndx}.sum_succ_rate;
        sim_cum.agg_exit_obs = sim_cum.agg_exit_obs + sim_out{pt_ndx}.exit_obs;
        sim_cum.agg_sum_exits = sim_cum.agg_sum_exits + sim_out{pt_ndx}.sum_exits;
        sim_cum.agg_mat_exit_moms_xx = sim_cum.agg_mat_exit_moms_xx + squeeze(sim_out{pt_ndx}.mat_exit_moms_xx);
        sim_cum.agg_mat_exit_moms_xy = sim_cum.agg_mat_exit_moms_xy + squeeze(sim_out{pt_ndx}.mat_exit_moms_xy);
        sim_cum.agg_mat_obs          = sim_cum.agg_mat_obs + sim_out{pt_ndx}.mat_obs;
        sim_cum.agg_nmat_exit        = sim_cum.agg_nmat_exit + sim_out{pt_ndx}.nmat_exit;
        sim_cum.agg_mat_exit_x       = [sim_cum.agg_mat_exit_x;sim_out{pt_ndx}.mat_exit_x];
        sim_cum.agg_mat_exit_y       = [sim_cum.agg_mat_exit_y;sim_out{pt_ndx}.mat_exit_y];
        sim_cum.agg_ship_obs    = sim_cum.agg_ship_obs    + sim_out{pt_ndx}.ship_obs;
        sim_cum.agg_ln_ships    = sim_cum.agg_ln_ships    + sim_out{pt_ndx}.ln_ships;
        sim_cum.agg_match_count = sim_cum.agg_match_count + squeeze(sim_out{pt_ndx}.match_count); %potential flip
        sim_cum.singletons      = sim_cum.singletons + sim_out{pt_ndx}.singletons;
    end

end        % end of prod_ndx loop

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

















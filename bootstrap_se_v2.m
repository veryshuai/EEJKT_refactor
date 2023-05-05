% %This script bootstraps standard errors for EEJKT
% %
% % Three steps:
% % 1. Simulate moments to recover covariance
% % 2. Find dM/dP, the dependence of moments on parameters
% % 3. Combine to recover standard erros
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Simulate moments many times
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = [-3.60529  -3.87941  0.23156  2.46487  1.88176  15.42634...
     0.38246  11.72211  1.38618  -1.21819  13.00238  -6.13506];

 
% %% First bootstrap the covariance matrix for the moments
% 
% N_boot_sims = 100;
% 
% X2params;
% SetParams;
% inten_sim_v1;
% 
% % simulate N_boot_sim sets of moments, holding X fixed
% % want different rands each time, set seed once here
% seed = 0;
% rng(80085,'twister');
% 
% Model_holder = zeros(38,N_boot_sims);
% 
% for k=1:N_boot_sims
% 
%     display(k);
%     
%     discrete_sim_parfor3;
%     
%     match_death_coefsSIM = [match_exit_rate;beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
%     match_ar1_coefsSIM   = [ybar_match;beta_match(2:4);mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]   
%     loglog_coefsSIM      = [b_degree]; % [intercept, slope, quadratic term]
%     mavshipSIM           = [avg_ln_ships]; % average ln(# shipments) 
%     exp_dom_coefsSIM     = [ybar_hfsales;beta_hfsales(2);mse_hf]; % [mean dep var.,coef,MSE]  
%     dom_ar1_coefsSIM     = [ybar_fsales_h;beta_fsales_h(2);mse_h]; % [mean dep var.,coef,MSE] 
%     match_lag_coefsSIM   = [mean_ln_haz;b_haz(2:6)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)] 
%     last_match_coefsSIM  = [mkt_exit_rate;beta_mkt_exit(2:6)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)]            
%     succ_rate_coefsSIM   = [mean_succ_rate;b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
%     sr_var_coefsSIM      = [mean_usq_succ;b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
%     for_sales_shrSIM     = [avg_expt_rate]; % mean share of exports to U.S. in total sales 
%     exp_fracSIM          = [share_exptr]; % fraction of firms exporting to U.S.  
% 
%     Model = cat(1,match_death_coefsSIM,match_ar1_coefsSIM,loglog_coefsSIM,...
%         mavshipSIM,exp_dom_coefsSIM,dom_ar1_coefsSIM,match_lag_coefsSIM,...   
%         last_match_coefsSIM,succ_rate_coefsSIM,sr_var_coefsSIM,for_sales_shrSIM,...    
%         exp_fracSIM);   
% 
%     Model_holder(:,k) = Model;
%     
% end
% 
% save results/boot_data_temp
% %% Now recover the gradient of the moments with respect to parameter changes
% 
load results/boot_data_temp

X2params;

param_vec    = [scale_h,F_h,F_f,scale_f,ah,bh,D_z,L_b,gam,cs_h,cs_f,sig_p];
param_names  = {'scale_h','F_h','F_f','scale_f','ah','bh','D_z','L_b','gam','cs_h','cs_f','sig_p'};

param_vec_ln   = [scale_h,F_h,F_f,scale_f,ah,bh,D_z,L_b,gam,log(cs_h),log(cs_f),sig_p];
param_names_ln = {'scale_h','F_h','F_f','scale_f','ah','bh','D_z','L_b','gam','ln_cs_h','ln_cs_f','sig_p'};

param_moments = zeros(size(Model,1),size(param_vec,2));
pct_del  = 3e-1;   % percentage change in parameters used to calcuate gradients
grad_del = pct_del * param_vec; % absolute changes in the size of parameters
logged_par = [10 11];  % flag coefficients that will be converted to logs for reporting

dMdP = zeros(size(param_moments));
dMdP_ln = zeros(size(param_moments));

for loop_ind = 0:size(X,2)

    %make a small change to param_vec
    param_vec_new = param_vec;
    if loop_ind ~= 0
        param_vec_new(loop_ind) = param_vec(loop_ind) + grad_del(loop_ind);
        
    end
    pCell = num2cell(param_vec_new); 
    [scale_h,F_h,F_f,scale_f,ah,bh,D_z,L_b,gam,cs_h,cs_f,sig_p] = pCell{:}; %insane that there is no better way to assign vector values to variables
    %sig_p = param_vec_new(end);
    
%   cs_h = exp(ln_cs_h);
%   cs_f = exp(ln_cs_f);

    %don't want randomness driving results
    %NOTE THIS IS HARD CODED IN DISCRETE SIM AT THE MOMENT, MUST MANUALLY
    %CHECK
    seed = 1;
    rng(80085,'twister');
    seed_crand(80085);
    
    %generate moments
    SetParams;
    inten_sim_v1;
    discrete_sim_parfor3;
    
    match_death_coefsSIM = [match_exit_rate;beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
    match_ar1_coefsSIM   = [ybar_match;beta_match(2:4);mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]   
    loglog_coefsSIM      = [b_degree]; % [intercept, slope, quadratic term]
    mavshipSIM           = [avg_ln_ships]; % average ln(# shipments) 
    exp_dom_coefsSIM     = [ybar_hfsales;beta_hfsales(2);mse_hf]; % [mean dep var.,coef,MSE]  
    dom_ar1_coefsSIM     = [ybar_fsales_h;beta_fsales_h(2);mse_h]; % [mean dep var.,coef,MSE] 
    match_lag_coefsSIM   = [mean_ln_haz;b_haz(2:6)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)] 
    last_match_coefsSIM  = [mkt_exit_rate;beta_mkt_exit(2:6)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)]            
    succ_rate_coefsSIM   = [mean_succ_rate;b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
    sr_var_coefsSIM      = [mean_usq_succ;b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
    for_sales_shrSIM     = [avg_expt_rate]; % mean share of exports to U.S. in total sales 
    exp_fracSIM          = [share_exptr]; % fraction of firms exporting to U.S.  

    Model = cat(1,match_death_coefsSIM,match_ar1_coefsSIM,loglog_coefsSIM,...
        mavshipSIM,exp_dom_coefsSIM,dom_ar1_coefsSIM,match_lag_coefsSIM,...   
        last_match_coefsSIM,succ_rate_coefsSIM,sr_var_coefsSIM,for_sales_shrSIM,...    
        exp_fracSIM);
    
    if loop_ind == 0
        base_moments = Model;
    else
        param_moments(:,loop_ind) = Model;

        dMdP(:,loop_ind) = (param_moments(:,loop_ind) - base_moments) / grad_del(loop_ind);
        % calculate derivatives of parameters expressed in logs
    end
end

% construct gradients with respect to logs of cost function parameters
dMdP_ln = dMdP;
dMdP_ln(:,logged_par) = dMdP(:,logged_par).* abs(param_vec(logged_par));

Mcov = cov(Model_holder');

% deal with moments excluded from estimation in distance.m
kp = ones(38,1);
kp(11) = 0;
kp(15) = 0;
kp(17) = 0;
kp(18) = 0;
kp(20) = 0;
kp(27:32)=0;

% throw out excluded moments and moments that don't move
keep = kp>0 & mean(Mcov)'~=0; 
momlist = 1:1:size(Model,1); % display moments kept
momlist(keep) % display moments kept
% strip rows and cols for unused moments out of covariance matrix
Mcov2 = Mcov(keep,:)';
Mcov3 = Mcov2(keep,:)';

G = -dMdP; %just change the notation to fit Joris' note 
G3 = -dMdP(keep,:);

G3_ln = -dMdP_ln(keep,:); % version with cost scalars in log form
[Data, W_inv] = target_stats(); %load in the weighting matrix
W = inv(W_inv);
W2 = W(keep,:)';
W3 = W2(keep,:)';
% 

V_std_ln3  = inv(G3_ln' * W3 * G3_ln) * G3_ln' * W3 * Mcov3 * W3 * G3_ln * inv(G3_ln' * W3 * G3_ln);
V_std3     = inv(G3' * W3 * G3) * G3' * W3 * Mcov3 * W3 * G3 * inv(G3' * W3 * G3);
% V_simple3  = inv(G3' * W3  * G3);
% V_pure_bs3 = inv(G3' * inv(Mcov3) * G3);

stderr3 = sqrt(diag(V_std3));
z_stat = param_vec'./stderr3;
se_table4 = table(param_names',param_vec',stderr3,z_stat);
% 
stderr3_ln = sqrt(diag(V_std_ln3));
z_stat_ln = param_vec_ln'./stderr3_ln;
se_table5 = table(param_names_ln',param_vec_ln',stderr3_ln,z_stat_ln);
% 
save results/bootstrap_results_benchmark
%%
% load results/bootstrap_results_benchmark
format shortG
se_table4
se_table5


AGS_sens_ln  = -inv(G3_ln' * W3 * G3_ln) * G3_ln' * W3;
AGS_elas_ln  = AGS_sens_ln .* repmat(base_moments(keep)',size(AGS_sens_ln,1),1) ./ repmat(param_vec_ln',1,size(AGS_sens_ln,2));
AGS_selas_ln  = AGS_sens_ln .* repmat(base_moments(keep)',size(AGS_sens,1),1);  % semi-elasticity wrt %change in moments


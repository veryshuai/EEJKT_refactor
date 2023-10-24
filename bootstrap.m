% This script bootstraps standard errors for EEJKT
% Basic tasks are
% (1) Generate bootstrapped moments
% (2) Generate finite differences
% (3) Generate covariance matrix for moments
% (4) Generate gradients for moments on parameters
% (5) Use delta method to get variance for parameters


clear all;

rng(80085,'twister');% set the seed and the rng (default parallel rand generator)
seed_crand(80085);

runs = 20; %20 %set number of runs to use in bootstrap
     
X =[-5.73144888007375	-16.1844310171305    0.109738437680579...
     0.242942844550223    0.605350029387283 11.0996895678257...
     0.0467130185363216   3.07428217391970	 2.20454924542622...
    11.1042769087445]; % fit = 10.9807
%  separate cost parameters, lambda_h = 2*lambda_f, same profit function scalars.

mm = setModelParameters(X);
mm.check_type = 108;

%solve the model 
policy = generatePolicyAndValueFunctions(mm);

boot_rand_list = floor(rand(runs,1) * 1e6); %set the random seed differently for each parallel run

%generate a cell for each moment
    match_death_coefs_boot_holder = cell(runs,1);% [match_exit_rate;beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
    match_ar1_coefs_boot_holder   = cell(runs,1);% [ybar_match;beta_match(2:4);mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]   
    mavship_boot_holder           = cell(runs,1);% [avg_ln_ships]; % average ln(# shipments) 
    exp_dom_coefs_boot_holder     = cell(runs,1);% [beta_hfsales(2)]; % [coef]  
    dom_ar1_coefs_boot_holder     = cell(runs,1);% [beta_fsales_h(2)]; % [coef] 
    match_lag_coefs_boot_holder   = cell(runs,1);% [mean_ln_haz;b_hazNew(2:2)]; % [mean dep. var, ln(1+a)] 
    succ_rate_coefs_boot_holder   = cell(runs,1);% [mean_succ_rate;b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
    sr_var_coefs_boot_holder      = cell(runs,1);% [mean_usq_succ;b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
    for_sales_shr_boot_holder     = cell(runs,1);% [avg_expt_rate]; % mean share of exports to U.S. in total sales 
    exp_frac_boot_holder          = cell(runs,1);% [share_exptr]; % fraction of firms exporting to U.S.  
    
%data moments (use these for sizing later)
[Data, W] = target_stats();

% Adjust moment vector and covariance matrix for exclusions
temp            = ones(1,size(Data,2));
[Data, temp, ~] = do_not_match_some_moments(Data, temp, temp);
MomsUsed        = temp > 0;
Data            = Data(:,MomsUsed);
temp            = W(MomsUsed,:);
W               = temp(:,MomsUsed);

for iter = 1:runs
    
    display(iter)
    
    rng(boot_rand_list(iter),'combRecursive');

    %%calculated moments to match from model (DONT FORGET TO COMMENT THE SEED AT THE TOP)
    simMoms = simulateMomentsMain(policy,mm);

    %Put the _boot_holder runs into their particular cells
    match_death_coefs_boot_holder{iter} = [simMoms.match_exit_rate;simMoms.beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
    match_ar1_coefs_boot_holder{iter}   = [simMoms.ybar_match;simMoms.beta_match(2:4);simMoms.mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]   
    mavship_boot_holder{iter}           = [simMoms.avg_ln_ships]; % average ln(# shipments) 
    exp_dom_coefs_boot_holder{iter}     = [simMoms.beta_hfsales(2)]; % [coef]  
    dom_ar1_coefs_boot_holder{iter}     = [simMoms.beta_fsales_h(2)]; % [mean dep var.,coef,MSE] 
    match_lag_coefs_boot_holder{iter}   = [simMoms.mean_ln_haz;simMoms.b_hazNew(2:end)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)] 
    succ_rate_coefs_boot_holder{iter}   = [simMoms.mean_succ_rate;simMoms.b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
    sr_var_coefs_boot_holder{iter}      = [simMoms.mean_usq_succ;simMoms.b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
    for_sales_shr_boot_holder{iter}     = [simMoms.avg_expt_rate]; % mean share of exports to U.S. in total sales 
    exp_frac_boot_holder{iter}          = [simMoms.share_exptr]; % fraction of firms exporting to U.S.  

end
  save results/moment_var

% (2) Now we need finite differences of parameters on moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fin_diff_size = (1 + 0.08); %percentage of parameter value 

param_vec = [mm.F_h, mm.scale_h, mm.ah, mm.bh, mm.D_z, mm.L_bF, mm.gam, mm.cs_h, mm.sig_p, mm.cs_f]';
pv_siz    = size(param_vec,1);

%each column is a finite differenced version of the parameter set
fin_diff_param_mat = repmat(param_vec,1,pv_siz);
fin_diff_param_mat(1:size(fin_diff_param_mat,1)+1:end) = (fin_diff_param_mat(1:size(fin_diff_param_mat,1)+1:end)) * fin_diff_size;
fin_diff_param_mat = [param_vec,fin_diff_param_mat]; %append "baseline" parameter set to difference with 

% %generate a cell for each moment
match_death_coefs_fd = cell(pv_siz + 1,1);% [match_exit_rate;beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
match_ar1_coefs_fd   = cell(pv_siz + 1,1);% [ybar_match;beta_match(2:4);mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]
loglog_coefs_fd      = cell(pv_siz + 1,1);% [b_degree]; % [intercept, slope, quadratic term]
mavship_fd           = cell(pv_siz + 1,1);% [avg_ln_ships]; % average ln(# shipments)
exp_dom_coefs_fd     = cell(pv_siz + 1,1);% [ybar_hfsales;beta_hfsales(2);mse_hf]; % [mean dep var.,coef,MSE]
dom_ar1_coefs_fd     = cell(pv_siz + 1,1);% [ybar_fsales_h;beta_fsales_h(2);mse_h]; % [mean dep var.,coef,MSE]
match_lag_coefs_fd   = cell(pv_siz + 1,1);% [mean_ln_haz;b_haz(2:6)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)]
succ_rate_coefs_fd   = cell(pv_siz + 1,1);% [mean_succ_rate;b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
sr_var_coefs_fd      = cell(pv_siz + 1,1);% [mean_usq_succ;b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
for_sales_shr_fd     = cell(pv_siz + 1,1);% [avg_expt_rate]; % mean share of exports to U.S. in total sales
exp_frac_fd          = cell(pv_siz + 1,1);% [share_exptr]; % fraction of firms exporting to U.S.
spb_fd               = cell(pv_siz + 1,1);% [simMoms.model_shareD]; % fraction of firms with 1 buyer, 2 buyers, etc.

for iter =1:pv_siz + 1

    display(iter)

    %point estimates
    mm = struct();
    cell_for_assignment = num2cell(fin_diff_param_mat(:,iter)); %need cell to unpack into multiple variables, matlab quirk
    [mm.F_h, mm.scale_h, mm.ah, mm.bh, mm.D_z, mm.L_bF, mm.gam, mm.cs_h, mm.sig_p, mm.cs_f] = cell_for_assignment{:};
    
%   set constrained parameters equal to their perturbed counterpart values    
    mm.F_f = mm.F_h;
    mm.scale_f =  mm.scale_h;
    mm.optimism = 0;
 
    rng(80085,'twister');% set the seed and the rng (default parallel rand generator)
    seed_crand(80085);
    mm = bootstrap_setModelParametersNoHead(X,mm);
    mm.check_type = 108;
    policy = generatePolicyAndValueFunctions(mm);
    rng(80085,'twister');% set the seed and the rng (default parallel rand generator)
    simMoms = simulateMomentsMain(policy,mm);

    %Put the _boot_holder runs into their particular cells
    match_death_coefs_fd{iter} = [simMoms.match_exit_rate;simMoms.beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
    match_ar1_coefs_fd{iter}   = [simMoms.ybar_match;simMoms.beta_match(2:4);simMoms.mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]   
    mavship_fd{iter}           = [simMoms.avg_ln_ships]; % average ln(# shipments) 
    exp_dom_coefs_fd{iter}     = [simMoms.beta_hfsales(2)]; % [coef]  
    dom_ar1_coefs_fd{iter}     = [simMoms.beta_fsales_h(2)]; % [coef] 
    match_lag_coefs_fd{iter}   = [simMoms.mean_ln_haz;simMoms.b_hazNew(2:end)]; % [mean dep. var, ln(1+a)] 
    succ_rate_coefs_fd{iter}   = [simMoms.mean_succ_rate;simMoms.b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
    sr_var_coefs_fd{iter}      = [simMoms.mean_usq_succ;simMoms.b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
    for_sales_shr_fd{iter}     = [simMoms.avg_expt_rate]; % mean share of exports to U.S. in total sales 
    exp_frac_fd{iter}          = [simMoms.share_exptr]; % fraction of firms exporting to U.S.  

end

%Reset to baseline values
cell_for_assignment = num2cell(fin_diff_param_mat(:,1)); %need cell to unpack into multiple variables
[mm.F_h, mm.scale_h, mm.ah, mm.bh, mm.D_z, mm.L_bF, mm.gam, mm.cs_h, mm.sig_p, mm.cs_f] = cell_for_assignment{:};


% (Probably don't need this block since we're not constructing std. errors for these parameters)
 mm.F_f = mm.F_h;
 mm.scale_f= mm.scale_h;
 mm.optimism = 0;

% (3) Construct moment covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Form moment vector, need to have vector to construct full covarience
%matrix
moments_boot_mat = zeros(size(Data,2),runs);
for iter = 1:runs
    moments_boot_mat(:,iter) = [... 
    match_death_coefs_boot_holder{iter};...
    match_ar1_coefs_boot_holder{iter}; ...
    mavship_boot_holder{iter};           ...
    exp_dom_coefs_boot_holder{iter};     ...
    dom_ar1_coefs_boot_holder{iter};     ...
    match_lag_coefs_boot_holder{iter};  ...
    succ_rate_coefs_boot_holder{iter};   ...
    sr_var_coefs_boot_holder{iter};      ...
    for_sales_shr_boot_holder{iter};     ...
    exp_frac_boot_holder{iter};          ...
    ];
end
Mcov = cov(moments_boot_mat'); %this is the estimated covariance matrix of the moments

% (4) Construct derivative of parameters on moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Finite differences
fin_diff_mat = zeros(size(Data,2),pv_siz + 1);
for iter = 1:pv_siz + 1
    fin_diff_mat(:,iter) = [...
        match_death_coefs_fd{iter}(:);...
        match_ar1_coefs_fd{iter}(:);...
        mavship_fd{iter}(:);...
        exp_dom_coefs_fd{iter}(:);...
        dom_ar1_coefs_fd{iter}(:);...
        match_lag_coefs_fd{iter}(:);...
        succ_rate_coefs_fd{iter}(:);...
        sr_var_coefs_fd{iter}(:);...
        for_sales_shr_fd{iter}(:);...
        exp_frac_fd{iter}(:);...
        ];
end

dMdP = zeros(size(Data,2),pv_siz);
for iter = 1:pv_siz
    dMdP(:,iter) = (fin_diff_mat(:,iter+1) - fin_diff_mat(:,1)) / ((fin_diff_size - 1) * param_vec(iter)); %this is our derivative approximation
end

% (5) Delta method to get each parameters variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use Joris lecture note 10.3

G = -dMdP; %just change the notation to fit Joris' note 

V = inv(G' * W * G) * G' * W * Mcov * W * G * inv(G' * W * G);

AGS_sens = -inv(G' * W * G) * G' * W;

AGS_elas = AGS_sens .* repmat((W * fin_diff_mat(:,1))',size(AGS_sens,1),1) ./ repmat(param_vec,1,size(AGS_sens,2));

% construct table for easy understanding of Gentzkow Shapiro weighting method

% AGS_param_names = {'F', 'scale_h', 'ah', 'bh', 'D_z', 'L_bF', 'gam', 'cs_h', 'sig_p', 'cs_f', 'scale_f'};
AGS_param_names = {'F', 'scale_h', 'ah', 'bh', 'D_z', 'L_bF', 'gam', 'cs_h', 'sig_p', 'cs_f'};

varnames = {'match_death_coefs1','match_death_coefs2','match_death_coefs3','match_death_coefs4','match_death_coefs5','match_ar1_coefs1','match_ar1_coefs2','match_ar1_coefs3','match_ar1_coefs4','match_ar1_coefs5','mavship','exp_dom_coefs1','dom_ar1_coefs2','match_lag_coefs1','match_lag_coefs2','succ_rate_coefs1','succ_rate_coefs2','sr_var_coefs1','sr_var_coefs2','for_sales_shr','exp_frac'};
AGS_cell = mat2cell(AGS_elas,ones(pv_siz,1),ones(size(Data,2),1));
AGS_table = cell2table(AGS_cell);
AGS_table.Properties.VariableNames = varnames;
AGS_table = [AGS_param_names',AGS_table];
writetable(AGS_table,'results/AGS_table.csv');

% Odds and ends, standard deviation and putting things into a readable
% formatfotex 
Psd = diag(V).^0.5;
std_error = Psd;
parameter = param_vec;
z_ratio = parameter./std_error;
format short
param_est = table(parameter,std_error,z_ratio,'Rownames',AGS_param_names)

save results/bootstrap_results

%% The following block handles the case where F is essentially 0 and shows no variation 

G2 = -dMdP(:,2:end); % drop first column, since no variation in F 

V2 = inv(G2' * W_alt * G2) * G2' * W_alt * Mcov * W_alt * G2 * inv(G2' * W_alt * G2);

AGS_param_names2 = {'scale_h', 'ah', 'bh', 'D_z', 'L_bF', 'gam', 'cs_h', 'sig_p', 'cs_f'};
Psd2 = diag(V2).^0.5;
std_error2 = Psd2;
param_vec2 = param_vec(2:end);
parameter2 = param_vec2;
z_ratio2 = parameter2./std_error2;
format short
param_est2 = table(parameter2,std_error2,z_ratio2,'Rownames',AGS_param_names2)


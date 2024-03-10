% This script bootstraps standard errors for EEJKT
% Basic tasks are
% (1) Replicate reported (benchmark) results
% (2) Generate bootstrapped covariance matrix for moments
% (3) Generate finite parameter differences
% (4) Generate gradients for moments with respect to parameters
% (5) Construct Andrews et al. sensitivity matrix 
% (6) Contruct covariance matrix for parameter estimates 

clear all;

runs = 20; %20 %set number of runs to use in bootstrap
SE_calc = 0; % set to 1 to get std. errors; to 0 for Andrews et. al calcs.
 
X = [-3.87377400411704,-19.6350579745564,0.141131266607483,0.224048944147957,...
    0.527541658446327,12.1673847502223,0.0462420845261331,5.13239605341601,...
    2.38418049835969,15.1614775287665]; % fit unix: 11.6935 PC: 11.6912

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% (1) Generate benchmark results
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%  Get baseline estimates and objects to be used. (This is basically the
%  distance.m function.)
    format long; 
    rng(80085,'twister');
    seed_crand(80085);    
    mm = setModelParameters(X);
    % choose the firm type to use for spot checking
    mm.check_type = 108;  
    %solve the model 
    policy = generatePolicyAndValueFunctions(mm);
    simMoms = simulateMomentsMain(policy,mm);
    testVec0 = [mm.F_h, mm.scale_h, mm.ah, mm.bh, mm.D_z, mm.L_bF, mm.gam, mm.cs_h, mm.sig_p, mm.cs_f]';
    [D,real_moms_and_sim_moms] = calculateDistanceAndPrint(simMoms,mm,X);


%data moments (use these for sizing later)
   [~, W] = target_stats();

% Adjust moment vector and covariance matrix for exclusions
    Data            = real_moms_and_sim_moms(:,1);
    Model           = real_moms_and_sim_moms(:,2);
    temp            = ones(size(Data,1),1);
    [Data, temp, ~] = do_not_match_some_moments(Data, temp, temp);
    MomsUsed        = temp > 0;
    Data            = Data(MomsUsed,1);
    Model           = Model(MomsUsed,1);
    temp            = W(MomsUsed,:);
    W               = temp(:,MomsUsed);
    
% Confirm transformed objects yield correct fit metric
    error           = (Data - Model);
    D2 = log((error'/W)*error);
    
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
% (2) Generate bootstrapped covariance matrix for moments  
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% Generate a cell for each moment
    match_death_coefs_boot_holder = cell(runs,1);% [match_exit_rate;beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
    match_ar1_coefs_boot_holder   = cell(runs,1);% [ybar_match;beta_match(2:4);mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]   
    mavship_boot_holder           = cell(runs,1);% [avg_ln_ships]; % average ln(# shipments) 
    exp_dom_coefs_boot_holder     = cell(runs,1);% [beta_hfsales(2)]; % [coef]  
    dom_ar1_coefs_boot_holder     = cell(runs,1);% [beta_fsales_h(2)]; % [coef] 
    match_haz_coefs_boot_holder   = cell(runs,1);% [mean_ln_haz;b_hazNew(2:2)]; % [mean dep. var, ln(1+a)] 
    succ_rate_coefs_boot_holder   = cell(runs,1);% [mean_succ_rate;b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
    sr_var_coefs_boot_holder      = cell(runs,1);% [mean_usq_succ;b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
    for_sales_shr_boot_holder     = cell(runs,1);% [avg_expt_rate]; % mean share of exports to U.S. in total sales 
    exp_frac_boot_holder          = cell(runs,1);% [share_exptr]; % fraction of firms exporting to U.S.  
 
    boot_rand_list = floor(rand(runs,1) * 1e6);  % set the random seed differently for each parallel run

parfor iter = 1:runs
    
    display(iter)
    
    rng(boot_rand_list(iter),'combRecursive');

    % Calculate moments to match from model using different seeds each iteration
    %(DONT FORGET TO COMMENT THE SEED AT THE TOP)
    simMoms = simulateMomentsMain(policy,mm);

    %Put the _boot_holder runs into their particular cells
    match_death_coefs_boot_holder{iter} = [simMoms.match_exit_rate;simMoms.beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
    match_ar1_coefs_boot_holder{iter}   = [simMoms.ybar_match;simMoms.beta_match(2:4);simMoms.mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]   
    mavship_boot_holder{iter}           = [simMoms.avg_ln_ships]; % average ln(# shipments) 
    exp_dom_coefs_boot_holder{iter}     = [simMoms.beta_hfsales(2)]; % [coef]  
    dom_ar1_coefs_boot_holder{iter}     = [simMoms.beta_fsales_h(2)]; % [mean dep var.,coef,MSE] 
    match_haz_coefs_boot_holder{iter}   = [simMoms.mean_ln_haz;simMoms.b_hazNew(2:end)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)] 
    succ_rate_coefs_boot_holder{iter}   = [simMoms.mean_succ_rate;simMoms.b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
    sr_var_coefs_boot_holder{iter}      = [simMoms.mean_usq_succ;simMoms.b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
    for_sales_shr_boot_holder{iter}     = [simMoms.avg_expt_rate]; % mean share of exports to U.S. in total sales 
    exp_frac_boot_holder{iter}          = [simMoms.share_exptr]; % fraction of firms exporting to U.S.  

end

% collect bootstrapped moment vectors in moments_boot_mat
moments_boot_mat = zeros(size(Data,1),runs);
for iter = 1:runs
    moments_boot_mat(:,iter) = [... 
    match_death_coefs_boot_holder{iter};...
    match_ar1_coefs_boot_holder{iter}; ...
    mavship_boot_holder{iter};           ...
    exp_dom_coefs_boot_holder{iter};     ...
    dom_ar1_coefs_boot_holder{iter};     ...
    match_haz_coefs_boot_holder{iter};  ...
    succ_rate_coefs_boot_holder{iter};   ...
    sr_var_coefs_boot_holder{iter};      ...
    for_sales_shr_boot_holder{iter};     ...
    exp_frac_boot_holder{iter};          ...
    ];
end
Mcov = cov(moments_boot_mat'); %this is the estimated covariance matrix of the moments
save results/moment_var

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
% (3) Generate set of perturbed parameter vectors (up and down)
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
load results/moment_var

if SE_calc == 1
% Some parameters are reported in logs. To get their std. errors,
% express fixed cost and cost function scalars in logs before generating differences
param_vec = [log(mm.F_h), mm.scale_h, mm.ah, mm.bh, mm.D_z, mm.L_bF, mm.gam, log(mm.cs_h), mm.sig_p, log(mm.cs_f)]';
else
param_vec = [mm.F_h, mm.scale_h, mm.ah, mm.bh, mm.D_z, mm.L_bF, mm.gam, mm.cs_h, mm.sig_p, mm.cs_f]';
end

pv_siz    = size(param_vec,1);
fin_diff_size1 = (1 + 0.05); % increase parameter values 

fin_diff_param_mat1 = repmat(param_vec,1,pv_siz);
fin_diff_param_mat1(1:size(fin_diff_param_mat1,1)+1:end) =...
            (fin_diff_param_mat1(1:size(fin_diff_param_mat1,1)+1:end)) * fin_diff_size1;        

fin_diff_size2 = (1 - 0.05); % decrease parameter value 
fin_diff_param_mat2 = repmat(param_vec,1,pv_siz);
fin_diff_param_mat2(1:size(fin_diff_param_mat2,1)+1:end) =...
            (fin_diff_param_mat2(1:size(fin_diff_param_mat2,1)+1:end)) * fin_diff_size2;        
     
% Concatenate "baseline" parameter vector with set of perturbed vectors         
fin_diff_param_mat = [param_vec,fin_diff_param_mat1,fin_diff_param_mat2]; 

if SE_calc == 1
% Restoring fixed cost and cost function scalars to levels
fin_diff_param_mat(1,:)  = exp(fin_diff_param_mat(1,:)); 
fin_diff_param_mat(8,:)  = exp(fin_diff_param_mat(8,:)); 
fin_diff_param_mat(10,:) = exp(fin_diff_param_mat(10,:)); 
end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
% (4) Generate gradients of moments and fit metric with respect to parameters
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 

% generate a cell for each moment vector
match_death_coefs_fd = cell(2*pv_siz + 1,1);% [match_exit_rate;beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
match_ar1_coefs_fd   = cell(2*pv_siz + 1,1);% [ybar_match;beta_match(2:4);mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]
mavship_fd           = cell(2*pv_siz + 1,1);% [avg_ln_ships]; % average ln(# shipments)
exp_dom_coefs_fd     = cell(2*pv_siz + 1,1);% [ybar_hfsales;beta_hfsales(2);mse_hf]; % [mean dep var.,coef,MSE]
dom_ar1_coefs_fd     = cell(2*pv_siz + 1,1);% [ybar_fsales_h;beta_fsales_h(2);mse_h]; % [mean dep var.,coef,MSE]
match_haz_coefs_fd   = cell(2*pv_siz + 1,1);% [mean_ln_haz;b_haz(2:6)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)]
succ_rate_coefs_fd   = cell(2*pv_siz + 1,1);% [mean_succ_rate;b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
sr_var_coefs_fd      = cell(2*pv_siz + 1,1);% [mean_usq_succ;b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
for_sales_shr_fd     = cell(2*pv_siz + 1,1);% [avg_expt_rate]; % mean share of exports to U.S. in total sales
exp_frac_fd          = cell(2*pv_siz + 1,1);% [share_exptr]; % fraction of firms exporting to U.S.
spb_fd               = cell(2*pv_siz + 1,1);% [simMoms.model_shareD]; % fraction of firms with 1 buyer, 2 buyers, etc.

Ddiff = zeros(2*pv_siz,1);
for iter =1:2*pv_siz + 1

    display(iter)

  % load parameter vectors
    cell_for_assignment = num2cell(fin_diff_param_mat(:,iter)); %need cell to unpack into multiple variables, matlab quirk
    [mm.F_h, mm.scale_h, mm.ah, mm.bh, mm.D_z, mm.L_bF, mm.gam, mm.cs_h, mm.sig_p, mm.cs_f] = cell_for_assignment{:};
    
  % set constrained parameters equal to their perturbed counterpart values    
    mm.F_f = mm.F_h;
    mm.scale_f =  mm.scale_h + 1;
    mm.optimism = 0;
 
    rng(80085,'twister');% set the seed and the rng (default parallel rand generator)
    seed_crand(80085);
    mm = bootstrap_setModelParametersNoHead(X,mm);
    testVec1 = [mm.F_h, mm.scale_h, mm.ah, mm.bh, mm.D_z, mm.L_bF, mm.gam, mm.cs_h, mm.sig_p, mm.cs_f]';
    policy1 = generatePolicyAndValueFunctions(mm);
    simMoms = simulateMomentsMain(policy1,mm);
  % retreive fit metric and moments for current parameter vector
    [Ddiff(iter),real_moms_and_sim_moms] = calculateDistanceAndPrint(simMoms,mm,X);
    
  % confirm that first parameter vector replicates benchmark results:
    if iter == 1
      Model = real_moms_and_sim_moms(MomsUsed,2);
      error = (Data - Model);
      D3 = log((error'/W)*error);
      test0 = sum(abs(testVec0 - testVec1));
      test1 = sum(sum(policy.c_val_f(:,:,1) - policy1.c_val_f(:,:,1)));
        try
        assert(D3-D==0);
        assert(test0==0);
        assert(test1==0);
        catch
        fprintf('Warning: Benchmark fit not replicated \n')
        end
    end
    
  % Put the simulated moments for the current vector into the relevant cells
    match_death_coefs_fd{iter} = [simMoms.match_exit_rate;simMoms.beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
    match_ar1_coefs_fd{iter}   = [simMoms.ybar_match;simMoms.beta_match(2:4);simMoms.mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]   
    mavship_fd{iter}           = [simMoms.avg_ln_ships]; % average ln(# shipments) 
    exp_dom_coefs_fd{iter}     = [simMoms.beta_hfsales(2)]; % [coef]  
    dom_ar1_coefs_fd{iter}     = [simMoms.beta_fsales_h(2)]; % [coef] 
    match_haz_coefs_fd{iter}   = [simMoms.mean_ln_haz;simMoms.b_hazNew(2:end)]; % [mean dep. var, ln(1+a)] 
    succ_rate_coefs_fd{iter}   = [simMoms.mean_succ_rate;simMoms.b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
    sr_var_coefs_fd{iter}      = [simMoms.mean_usq_succ;simMoms.b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
    for_sales_shr_fd{iter}     = [simMoms.avg_expt_rate]; % mean share of exports to U.S. in total sales 
    exp_frac_fd{iter}          = [simMoms.share_exptr]; % fraction of firms exporting to U.S.  

end

% average the vector of changes in the fit metric with respect to parameter upshocks and downshocks.
  AvgDdiff = (Ddiff(2:pv_siz+1)+Ddiff(pv_siz+2:end))./2 - Ddiff(1);

% Reset parameter vector to benchmark values
cell_for_assignment = num2cell(fin_diff_param_mat(:,1)); %need cell to unpack into multiple variables
[mm.F_h, mm.scale_h, mm.ah, mm.bh, mm.D_z, mm.L_bF, mm.gam, mm.log_cs_h, mm.sig_p, mm.log_cs_f] = cell_for_assignment{:};

% (Don't need this block since we're not constructing std. errors for these parameters)
%  mm.F_f = mm.F_h;
%  mm.scale_f= mm.scale_h +1;
%  mm.optimism = 0;
 
% Load moments at each perturbed parameter vector into fin_diff_mat
fin_diff_mat = zeros(size(Data,1),2*pv_siz + 1);
for iter = 1:2*pv_siz + 1
    fin_diff_mat(:,iter) = [...
        match_death_coefs_fd{iter}(:);...
        match_ar1_coefs_fd{iter}(:);...
        mavship_fd{iter}(:);...
        exp_dom_coefs_fd{iter}(:);...
        dom_ar1_coefs_fd{iter}(:);...
        match_haz_coefs_fd{iter}(:);...
        succ_rate_coefs_fd{iter}(:);...
        sr_var_coefs_fd{iter}(:);...
        for_sales_shr_fd{iter}(:);...
        exp_frac_fd{iter}(:);...
        ];
end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
% (5) Construct Andrews sensitivity matrix
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
% approximate moment derivatives using finite differences
  dMdP_up = zeros(size(Data,1),pv_siz);
  dMdP_down = zeros(size(Data,1),pv_siz);
  for iter = 1:pv_siz
    dMdP_up(:,iter) =...
    (fin_diff_mat(:,iter+1) - fin_diff_mat(:,1)) / ((fin_diff_size1 - 1) * param_vec(iter)); 
    dMdP_down(:,iter) =...
    (fin_diff_mat(:,pv_siz+iter+1) - fin_diff_mat(:,1)) / ((fin_diff_size2 - 1) * param_vec(iter)); 
  end

% average derivatives based on upshocks and downshocks
  dMdP = (dMdP_up + dMdP_down)./2 ;
  G = -dMdP; % change the notation to fit Joris' note 

% Construct Andrews et al. sensitivity matrix
  AGS_sens = -inv(G' * W * G) * G' * W;
% Convert matrix to alternative elasticity form
  AGS_elas = AGS_sens .* repmat((W * fin_diff_mat(:,1))',size(AGS_sens,1),1) ./ repmat(param_vec,1,size(AGS_sens,2));

% construct table for easy understanding of Gentzkow Shapiro weighting method
  AGS_param_names = {'F', 'log(Pi)', 'ah', 'bh', 'D_z', 'L_bF', 'gam', 'log(cs_h)', 'sig_p', 'log(cs_f)'};
  varnames = {'match_death_coefs1','match_death_coefs2','match_death_coefs3','match_death_coefs4','match_death_coefs5','match_ar1_coefs1','match_ar1_coefs2','match_ar1_coefs3','match_ar1_coefs4','match_ar1_coefs5','mavship','exp_dom_coefs1','dom_ar1_coefs2','match_lag_coefs1','match_lag_coefs2','succ_rate_coefs1','succ_rate_coefs2','sr_var_coefs1','sr_var_coefs2','for_sales_shr','exp_frac'};
  AGS_cell = mat2cell(AGS_elas,ones(pv_siz,1),ones(size(Data,1),1));
  AGS_table = cell2table(AGS_cell);
  AGS_table.Properties.VariableNames = varnames;
  AGS_table = [AGS_param_names',AGS_table];
  writetable(AGS_table,'results/AGS_table.csv');

  AvgDdiffcell = mat2cell(AvgDdiff,ones(pv_siz,1),1);
  AvgDiff_table = cell2table(AvgDdiffcell);
  AvgDiff_table = [AGS_param_names',AvgDiff_table];
  writetable(AvgDiff_table,'results/AvgDiff_table.csv');

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
% (6) Construct covariance matrix for estimated parameters 
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 

% use delta method to get covariance matrix
  V = inv(G' * W * G) * G' * W * Mcov * W * G * inv(G' * W * G);
  Psd = diag(V).^0.5;
  std_error = Psd;
  parameter = param_vec;
  z_ratio = parameter./std_error;
  format short
  param_est = table(parameter,std_error,z_ratio,'Rownames',AGS_param_names)

  
  
  save results/bootstrap_results


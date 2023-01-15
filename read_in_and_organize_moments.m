function [Data, W, Model] = read_in_and_organize_moments(simMoms)
[Data, W] = target_stats();

    match_death_coefsSIM = [simMoms.match_exit_rate;simMoms.beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
    match_ar1_coefsSIM   = [simMoms.ybar_match;simMoms.beta_match(2:4);simMoms.mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]   
    loglog_coefsSIM      = [simMoms.b_degree]; % [intercept, slope, quadratic term]
    mavshipSIM           = [simMoms.avg_ln_ships]; % average ln(# shipments) 
    exp_dom_coefsSIM     = [simMoms.ybar_hfsales;simMoms.beta_hfsales(2);simMoms.mse_hf]; % [mean dep var.,coef,MSE]  
    dom_ar1_coefsSIM     = [simMoms.ybar_fsales_h;simMoms.beta_fsales_h(2);simMoms.mse_h]; % [mean dep var.,coef,MSE] 
    ln_haz_coefsSIM      = [simMoms.mean_ln_haz;simMoms.b_haz(2:6)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)] 
    succ_rate_coefsSIM   = [simMoms.mean_succ_rate;simMoms.b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
    sr_var_coefsSIM      = [simMoms.mean_usq_succ;simMoms.b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
    for_sales_shrSIM     = [simMoms.avg_expt_rate]; % mean share of exports to U.S. in total sales 
    exp_fracSIM          = [simMoms.share_exptr]; % fraction of firms exporting to U.S. 
    
    degree_distSIM_D     = simMoms.model_shareD; % fraction of firms with 1 buyer, 2 buyers, etc.

 Model = cat(1,match_death_coefsSIM,match_ar1_coefsSIM,loglog_coefsSIM,...
    mavshipSIM,exp_dom_coefsSIM,dom_ar1_coefsSIM,ln_haz_coefsSIM,...   
    succ_rate_coefsSIM,sr_var_coefsSIM,for_sales_shrSIM,...    
    exp_fracSIM);
end
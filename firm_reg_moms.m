function  [moms_xx,moms_xy,ysum,n_obs] = firm_reg_moms(iter_in,mm)

% JT: This is a match level regression. Why is it called firm_reg_moms?
    
obs2use = (iter_in.mat_cont_2yr(:,2).*iter_in.mat_cont_2yr(:,9))>0;
% matches: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age, firm age]  
% mat_cont_2yr: [matches_lagged, matches] spliced, continuing matches only

data_mat = iter_in.mat_cont_2yr(obs2use,:);
n_obs = size(data_mat,1);

if n_obs>0
    
    ln_export     = log(data_mat(:,9));
    ln_export_lag = log(data_mat(:,2));
    first_yr      = data_mat(:,6) < mm.pd_per_yr;
    ln_age        = log(data_mat(:,6)./12);
    const         = ones(n_obs,1);

    y = ln_export ;

    x = [const,first_yr,ln_export_lag,ln_age];
    moms_xx = x'*x;
    moms_xy = x'*y;
    ysum = sum(y);

else
        
    kk = 4; % columns of x matrix
    x = double.empty(0,4);
    moms_xx = zeros(4,4);
    moms_xy = zeros(4,1);
    ysum    = 0;
    n_obs = 0;
    
end



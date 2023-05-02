function  [x,y,moms_xx,moms_xy,ysum,n_obs] = match_reg_moms(mat_cont,ncols,mm)

     % matches: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age in periods, firm age in periods]  
     % mat_cont: [matches_lagged, matches] spliced, continuing matches only
     
    ff_obs = find(mat_cont(:,2).*mat_cont(:,ncols+2)>0); % obs. with positive current and lagged sales
    dat = mat_cont(ff_obs,:);
    n_obs = size(ff_obs,1);
    y = log(dat(:,ncols+2));         % log current year sales
    ln_age = log(dat(:,ncols)/mm.pd_per_yr + 1); % log firm age last year (treating data-based age as # years)
    x = [ones(n_obs,1),log(dat(:,2)),dat(:,4)==0,ln_age];  % constant, log(lagged sales), 1st yr. match dummy, firm age
    moms_xx = x'*x;
    moms_xy = x'*y;
    ysum = sum(y);
    
  
end



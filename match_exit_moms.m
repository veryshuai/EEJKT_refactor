% Called from matchdat_gen_f.m, this function constructs the t- and type-specific
% moments needed for match exit regressions. It passes them back to
% matchdat_gen_f.m, which aggregates to the year level and passes them back 
% to discrete_sim_parfor3.m for aggregation across types.

function [x,y,mat_exit_moms_xx,mat_exit_moms_xy,nobs,nmatch_exit] = match_exit_moms(matches,pd_per_yr)
 % matches: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age, firm age] 

ff = find(matches(:,2)>0);
y  = matches(ff,5)==0;                % match dead by end of year
x1 = matches(ff,4)==0;                % new match this year 
x2 = log(matches(ff,2));              % sales during year
x3 = log(1+matches(ff,6)./pd_per_yr); % log age of match, in years
x4 = log(1+matches(ff,7)./pd_per_yr); % log age of exporter, in years
x0 = ones(size(ff,1),1);
x  = [x0,x1,x2,x3,x4];
nobs = size(matches,1);
nmatch_exit = sum(y);
mat_exit_moms_xx = x'*x;
mat_exit_moms_xy = x'*y;

end
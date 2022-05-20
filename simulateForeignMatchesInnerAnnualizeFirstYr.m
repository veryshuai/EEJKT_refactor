function iter_in = simulateForeignMatchesInnerAnnualizeFirstYr(iter_in)
% convert age in periods to age in years for first-year observations.
% After the first year, conversion to years handled by mat_yr_splice
if iter_in.year==1
    iter_in.mat_yr_sales(:,6)  = iter_in.mat_yr_sales(:,6) > 0;  % set match age in years to 1 if match age in periods > 0
    iter_in.mat_yr_sales(:,7)  = iter_in.mat_yr_sales(:,7) > 0;  % set firm age in years to 1 if year age in periods > 0
    iter_in.firm_yr_sales(:,4) = iter_in.firm_yr_sales(:,4) > 0; % set firm age in years to 1 if year age in periods > 0
    iter_in.year_lag           = 1;
end
end
function [iter_in,iter_out] = simulateForeignMatchesInnerAnnualize(iter_in,iter_out,mm)

[iter_in.mat_yr_sales,iter_in.firm_yr_sales] =  season_merge(iter_in,mm);

% mat_yr_sales:  [firm ID, match-specific sales, shipments, boy Z, eoy Z,
%                 match age in periods (w/in year), firm age in periods]

% firm_yr_sales: [firmID,sales,#shipments,firm age]

iter_in = simulateForeignMatchesInnerAnnualizeFirstYr(iter_in); %since no lag, first year is different
iter_out = simulateForeignMatchesInnerAnnualizeDuds(iter_in, mm, iter_out); 
% # unsuccessful meetings (duds) over previous year, by firm (needed for degree distribution later)
[iter_in, iter_out] = simulateForeignMatchesInnerAnnualizeMeetingGaps(iter_in, mm, iter_out);

if iter_in.year > 2
    if size(iter_in.mat_yr_sales,1)*size(iter_in.mat_yr_sales_lag)>0 % drops types without sales in both current and lagged years

      iter_in.ncols = size(iter_in.mat_yr_sales,2);
      
      [iter_in.mat_cont_2yr,iter_in.mat_yr_sales,iter_in.mat_yr_sales_adj,iter_in.year_lag]...
          = mat_yr_splice_v2(iter_in.mat_yr_sales,iter_in.mat_yr_sales_lag,mm,iter_in.year_lag,iter_in.year);
      %  mat_cont_2yr: [mat_yr_sales_lag(ff_cont_lag,:), mat_yr_sales(ff_cont,:)]
      %  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age in yrs, firm age in yrs]
      %  mat_yr_sales_adj: same as lagged mat_yr_sales except eoy Z set to zero if no sales next year

      iter_in.mat_matur_dat...
          = [iter_in.mat_yr_sales(:,2),iter_in.mat_yr_sales(:,4:7),...
           iter_in.mat_yr_sales(:,1),ones(size(iter_in.mat_yr_sales,1),1).*iter_in.t/mm.pd_per_yr] ;
      % agg_mat_matur: [sales, boy Z, eoy Z, match age, fmic_typeirm age, firm_ID, yr]
    end
end

end
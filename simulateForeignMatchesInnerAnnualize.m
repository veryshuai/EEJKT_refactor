function [iter_in,iter_out] = simulateForeignMatchesInnerAnnualize(iter_in,iter_out,mm)

[iter_in.mat_yr_sales,iter_in.firm_yr_sales,iter_in.Zcut_eoy] =  season_merge(iter_in,mm);

% mat_yr_sales:  [firm ID, match-specific sales, shipments, boy Z, eoy Z,
%                 match age in periods (w/in year), firm age in periods]

% firm_yr_sales: [firmID,sales,#shipments,firm age]

iter_in = simulateForeignMatchesInnerAnnualizeFirstYr(iter_in); 
% since no lag, first year is different

iter_out = simulateForeignMatchesInnerAnnualizeDuds(iter_in, mm, iter_out); 
% # unsuccessful meetings (duds) over previous year, by firm (needed for degree distribution later)

iter_in.cur_duds(:,iter_in.t-mm.pd_per_yr+1:iter_in.t) = iter_out.cur_yr_duds;
% load this year's dud counts into the cur_duds series for all firm_IDs

[iter_in, iter_out] = simulateForeignMatchesInnerAnnualizeMeetingGaps(iter_in, mm, iter_out);

if iter_in.year > 2
    if size(iter_in.mat_yr_sales,1)*size(iter_in.mat_yr_sales_lag,1)>0 % drops types without sales in both current and lagged years

      iter_in.ncols = size(iter_in.mat_yr_sales,2);
            
      [iter_in.mat_cont_2yr,iter_in.mat_yr_sales,iter_in.mat_yr_sales_lag,iter_in.year_lag]...
          = mat_yr_splice_v2(iter_in,mm,iter_in.year);
      %  mat_cont_2yr: [mat_yr_sales_lag(ff_cont_lag,:), mat_yr_sales(ff_cont,:)]
      %  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age in yrs, firm age in yrs]

    end
end

end
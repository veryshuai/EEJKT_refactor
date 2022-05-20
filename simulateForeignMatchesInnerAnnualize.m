function [iter_in,iter_out] = simulateForeignMatchesInnerAnnualize(iter_in,iter_out,mm)

[iter_in.mat_yr_sales,iter_in.firm_yr_sales] =...
    season_merge(iter_in.seas_tran,iter_in.N_match,mm.sim_firm_num_by_prod_succ_type(iter_in.pt_ndx),mm.pd_per_yr);
% mat_yr_sales:  [firm ID, match-specific sales, shipments, boy Z, eoy Z,
%                 match age in periods (w/in year), firm age in periods]
% firm_yr_sales: [firmID,sales,#shipments,firm age]

% convert age in periods to age in years for first-year observations.
% After the first year, conversion to years handled by mat_yr_splice
if iter_in.year==1
    iter_in.mat_yr_sales(:,6)  = iter_in.mat_yr_sales(:,6) > 0;  % set match age in years to 1 if match age in periods > 0
    iter_in.mat_yr_sales(:,7)  = iter_in.mat_yr_sales(:,7) > 0;  % set firm age in years to 1 if year age in periods > 0
    iter_in.firm_yr_sales(:,4) = iter_in.firm_yr_sales(:,4) > 0; % set firm age in years to 1 if year age in periods > 0
    iter_in.year_lag           = 1;
end

% # unsuccessful meetings (duds) over previous year, by firm (needed for degree distribution later)
yr_tlag = iter_in.t-mm.pd_per_yr;
cum_duds  = iter_in.cum_meets(:,yr_tlag:iter_in.t) - iter_in.cum_succ(:,yr_tlag:iter_in.t); % previous mm.pd_per_yr + 1 cumulative duds
curr_duds = cum_duds(:,2:mm.pd_per_yr+1)-(iter_in.new_firm(:,yr_tlag+1:iter_in.t)==0).*cum_duds(:,1:mm.pd_per_yr);

iter_out.singletons = sum(sum(curr_duds));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Generate observations on intervals between meetings, cumulative meetings, and cumulative succcesses

if iter_in.t>3*mm.pd_per_yr
    [time_gap,iter_in.mkt_exit] = time_gaps(iter_in.t,iter_in.exit_firm,mm.pd_per_yr,iter_in.cum_meets,iter_in.cum_succ);
    %           time_gap: (1) firm_ID, (2) periods into interval, (3) time gap,
    %           (4) # new meetings, (5) iter_in.t, (6) cum. meetings, (7) cum succeses

    iter_out.time_gaps = [iter_out.time_gaps;time_gap];
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iter_in.year > 1  % construct moments for match regressions, firm export regressions

    if size(iter_in.mat_yr_sales,1)*size(iter_in.mat_yr_sales_lag)>0
        % drops types without sales in both current and lagged years

        % mat_yr_splice splices consecutive obs. on annualized data and updates
        % match ages and firm ages, converting them to years. It's really costly!

        iter_in.ncols = size(iter_in.mat_yr_sales,2);

        [iter_in.mat_cont_2yr,iter_in.mat_yr_sales,iter_in.mat_yr_sales_adj,iter_in.year_lag] =...
            mat_yr_splice_v2(iter_in.mat_yr_sales,iter_in.mat_yr_sales_lag,mm,iter_in.year_lag,iter_in.year);

        %  mat_cont_2yr: [mat_yr_sales_lag(ff_cont_lag,:), mat_yr_sales(ff_cont,:)]
        %  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age in yrs, firm age in yrs]
        %  iter_in.mat_yr_sales_adj: same as lagged mat_yr_sales except eoy Z set to zero if no sales next year

        %% the following matrices accumulate annualized values over time and firm types

        iter_in.mat_matur_dat =  [iter_in.mat_yr_sales(:,2),iter_in.mat_yr_sales(:,4:7),iter_in.mat_yr_sales(:,1),ones(size(iter_in.mat_yr_sales,1),1).*iter_in.t/mm.pd_per_yr] ;
        % agg_mat_matur: [sales, boy Z, eoy Z, match age, fmic_typeirm age, firm_ID, yr]
    end
end



end
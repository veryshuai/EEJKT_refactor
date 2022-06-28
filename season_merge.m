function [mat_yr_sales,firm_yr_sales] = season_merge(iterX_in,mm)

%  This function takes a year's worth of season-specific match
%  outcomes for a particular type of firm and organizes the information
%  into annual aggregates at the match and firm level. 

%  iterX_in.seas_tran: [t,season,year,initial state,firm_ID,new state,rev,shipments,firm age (in periods)];

%% build within-yr trajectories for all matches--continuing, new, and dying

[mat_cols, all_seas, som_seas] = season_mergeWithinYrSequence(mm, iterX_in);

%%  Package up the match info. for use in regressions

[mat_yr_sales, firm_yr_sales] = season_mergeAnnualizeDat(all_seas, som_seas, mm, mat_cols);

end


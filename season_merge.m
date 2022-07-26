function [mat_yr_sales,firm_yr_sales,iterX_in] = season_merge(iterX_in,mm)

%  This function takes a year's worth of season-specific match
%  outcomes for a particular type of firm and organizes the information
%  into annual aggregates at the match and firm level. 

%  iterX_in.seas_tran: [t,season,year,initial state,firm_ID,new state,rev,shipments,firm age (in periods)];

%% build within-yr monthly trajectories for all matches--continuing, new, and dying

[mat_cols, all_seas, som_seas] = season_mergeWithinYrSequence(mm, iterX_in);

 % all_seas and som_seas 
 %  (1) t, (2) season, (3) year, (4) initial state, (5) exporter id, (6) ending state,
 %  (7) match revenue,(8) #shipments,(9) exporter age (#periods), (10) match age w/in year


%%  Package up the match info. for use in regressions

[mat_yr_sales, firm_yr_sales,iterX.in] = season_mergeAnnualizeDat(all_seas, som_seas, mm, mat_cols,iterX_in);

% mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z,
%                match age in periods (w/in year), firm age in periods] 

% firm_yr_sales: [firmID,sales,#shipments,firm age]
 
end


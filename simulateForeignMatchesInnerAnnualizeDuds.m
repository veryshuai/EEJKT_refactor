function iter_out = simulateForeignMatchesInnerAnnualizeDuds(iter_in, mm, iter_out)

% This function is called, firm type by firm type, at the end of each year.
% It infers "duds" within the year from the difference between increments to 
% the cumulative number of meetings and increments to the cumulative number
% of successes. 

yr_tlag   = iter_in.t-mm.pd_per_yr;

% cumulative number of duds, by month
cum_duds  = iter_in.cum_meets(:,yr_tlag:iter_in.t) - iter_in.cum_succ(:,yr_tlag:iter_in.t); 

% current number of duds, by month, resetting cumulative number for new entrants
iter_out.cur_yr_duds  = ...
  cum_duds(:,2:mm.pd_per_yr+1)-(iter_in.new_firm(:,yr_tlag+1:iter_in.t)==0).*cum_duds(:,1:mm.pd_per_yr);

% iter_out.singletons = sum(sum(iter_out.cur_yr_duds));
end
function iter_out = simulateForeignMatchesInnerAnnualizeDuds(iter_in, mm, iter_out)
yr_tlag = iter_in.t-mm.pd_per_yr;
cum_duds  = iter_in.cum_meets(:,yr_tlag:iter_in.t) - iter_in.cum_succ(:,yr_tlag:iter_in.t); % previous mm.pd_per_yr + 1 cumulative duds
curr_duds = cum_duds(:,2:mm.pd_per_yr+1)-(iter_in.new_firm(:,yr_tlag+1:iter_in.t)==0).*cum_duds(:,1:mm.pd_per_yr);

iter_out.singletons = sum(sum(curr_duds));
end
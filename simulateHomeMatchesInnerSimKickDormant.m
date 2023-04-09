function [iterH_in] = simulateHomeMatchesInnerSimKickDormant(iterH_in, mm)

%%              simulateHomeMatchesInnerSimKickDormant
%    Deal with dormant firms: after one year with no clients, swap them out

   t = iterH_in.t;

   if t > 2*mm.pd_per_yr
        dormant = sum(iterH_in.cur_cli_cnt(:,t - mm.pd_per_yr+1:t),2)==0; % no clients for 1 year
        iterH_in.micro_state(dormant,t)  = 1;  % reset initial micro state to 1 (entrant)
        iterH_in.new_firm(dormant,t)     = 1;  % will mark firms that haven't had an active match in 2 yrs
        iterH_in.cum_succ(dormant,t)     = 0;  % cumulative number of successes
        iterH_in.cur_cli_zst(dormant,:)  = 0;  % NEW: client counts by Z type
    end
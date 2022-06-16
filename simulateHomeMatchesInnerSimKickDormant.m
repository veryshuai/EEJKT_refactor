function [iterH_in] = simulateHomeMatchesInnerSimKickDormant(iterH_in, mm)

%%              simulateHomeMatchesInnerSimKickDormant
%    Deal with dormant firms: after two years with no clients, swap them out

   t = iterH_in.t;

   if t > 2*mm.pd_per_yr
        dmt = sum(iterH_in.cur_cli_cnt(:,t-2*mm.pd_per_yr+1:t),2)==0; % identify dormant firms
        iterH_in.micro_state(dmt,t)  = 1;  % reset initial micro state to 1 (entrant)
        iterH_in.new_firm(dmt,t)     = 1;  % will mark firms that haven't made a match in 2 yrs
        iterH_in.cum_succ(dmt,t)     = 0;  % cumulative number of successes
    end
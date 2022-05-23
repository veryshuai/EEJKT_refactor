function iter_in = simulateForeignMatchesInnerSimKickDormant(iter_in, mm)

    % Dormant firms have no clients for at least one year
    if iter_in.t > mm.pd_per_yr
        dormant = sum(iter_in.cur_cli_cnt(:,iter_in.t-mm.pd_per_yr+1:iter_in.t),2)==0; % identify dormant firms
        iter_in.micro_state(dormant,iter_in.t)  = 1;  % reset initial micro state to 1 (entrant)
        iter_in.new_firm(dormant,iter_in.t)     = 1;  % will mark new firms that haven'iter_in.t made a match
%       iter_in.exit_firm(dormant,iter_in.t)    = 1;  % will mark last period of exiting firm (same thing)
        iter_in.cum_meets(dormant,iter_in.t)    = 0;  % cumulative number of meetings
        iter_in.cum_succ(dormant,iter_in.t)     = 0;  % cumulative number of successes
    end
end
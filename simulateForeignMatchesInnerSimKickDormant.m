function iter_in = simulateForeignMatchesInnerSimKickDormant(iter_in, mm)

    if iter_in.t > mm.pd_per_yr
        dormant = sum(iter_in.cur_cli_cnt(:,iter_in.t-mm.pd_per_yr+1:iter_in.t),2) + iter_in.add_cli_cnt(:,iter_in.t)==0  ; % no clients for one year
        iter_in.micro_state(dormant,iter_in.t)   = 1;  % reset initial micro state to 1 (entrant)
        iter_in.new_firm(dormant,iter_in.t)      = 1;  % will mark new firms that haven't made a match yet
        iter_in.exit_firm(dormant,iter_in.t-1)   = 1;
        iter_in.cum_meets(dormant,iter_in.t)     = 0;  % cumulative number of meetings
        iter_in.cum_succ(dormant,iter_in.t)      = 0;  % cumulative number of successe        
        iter_in.cur_cli_zst(dormant,:)           = 0;  % NEW: client counts by Z type
    end


end
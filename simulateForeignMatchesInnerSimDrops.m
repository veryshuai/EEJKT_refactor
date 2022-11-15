function [iter_in, drop_Zcut, drop_cnt] = simulateForeignMatchesInnerSimDrops(iter_in, policy, mm)

    % identify z values at which exporters keep current matches from t-1 to t
    iter_in.keep_cli = policy.c_val_f(:,mm.pt_type(iter_in.pt_ndx,1),iter_in.macro_state_f(iter_in.t))' > 0; % = 1 if want to keep type for t
%    iter_in.keep_cli = policy.c_val_f(:,mm.pt_type(iter_in.pt_ndx,1),iter_in.macro_state_f(iter_in.t-1))' > 0; % = 1 if want to keep type for t
    drop_Zcut = size(mm.Z,1) - sum(iter_in.keep_cli); % cutoff: matches dropped at z value <= drop_Zcut

    % count endogenous drops for all exporter hotel rooms (exporters): z too low to continue
    drop_cnt = sum(iter_in.lag_cli_zst.*(1-iter_in.keep_cli),2);
    % NOTE: drops based on lagged Z because current Z not generated yet at
    % this point in the time loop.

    % draw the number of exogenous deaths of remaining matches between t-1 and t
    ddum = find(iter_in.cur_cli_cnt(:,iter_in.t-1)-drop_cnt > 0);
    if sum(ddum)>0
        iter_in.exog_deaths(ddum,iter_in.t-1) =...
            random('bino',iter_in.cur_cli_cnt(ddum,iter_in.t-1)-drop_cnt(ddum),1-exp(-mm.delta));
    end
end
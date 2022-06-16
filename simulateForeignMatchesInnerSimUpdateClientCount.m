function [iter_in, drop_Zcut] = simulateForeignMatchesInnerSimUpdateClientCount(iter_in, mm, policy)
    
    no_learn = iter_in.cum_meets(:,iter_in.t-1) >= mm.n_size-3; % for picking off seasoned exporters (up to 3 trials before maximum learning)
    learn    = iter_in.cum_meets(:,iter_in.t-1) <  mm.n_size-3; % for picking off exporters who are still learning
    stay     = ones(mm.sim_firm_num_by_prod_succ_type(iter_in.pt_ndx),1);  % will flag no-learning firms that continue, t-1 to t
    [stay, iter_in] = simulateForeignMatchesInnerSimLearners(learn, iter_in, policy, stay);
    iter_in.add_cli_cnt(:,iter_in.t) = max(iter_in.cum_succ(:,iter_in.t) - iter_in.cum_succ(:,iter_in.t-1),0); %max() resets count to 0 for neg. differences to deal with new exporter
    iter_in.new_firm(:,iter_in.t)    = max((iter_in.cum_meets(:,iter_in.t) - iter_in.cum_meets(:,iter_in.t-1) < 0),(stay == 0)); %first period for new exporter
    iter_in.exit_firm(:,iter_in.t-1) = iter_in.cum_meets(:,iter_in.t) - iter_in.cum_meets(:,iter_in.t-1) < 0 ; % last period before exit
    % DAVID QUESTION: Does exit firm pick up no learning firms that exit?
    % If so, why is stay==0 condition needed?
    % JIM ANSWER: stay==0 throws out no-learning firms that get an exogenous exit shock. (This isn't needed for the learning firms
    % becasue exit shocks are built into the Q matrix, and exit always means cum meets(t)<=cum_meets(t-1).) 

    [iter_in, drop_Zcut, drop_cnt] = simulateForeignMatchesInnerSimDrops(iter_in, policy, mm);
    iter_in.cur_cli_cnt(:,iter_in.t) = iter_in.add_cli_cnt(:,iter_in.t) + iter_in.cur_cli_cnt(:,iter_in.t-1) ...
        - drop_cnt - iter_in.exog_deaths(:,iter_in.t-1) ;
end
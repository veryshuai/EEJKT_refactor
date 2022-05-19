function [iter_in, stay] = simulateForeignMatchesInnerSimLearners(learn, iter_in, policy, stay)
    
    % Find N_learn randomly selected micro states next period, given
    % macro state (common to all firms), initial micro states, and pmat_cum.
    N_learn = sum(learn,1) ;
    if  N_learn > 0
        trans_rands = iter_in.pmat_cum_t(iter_in.micro_state(learn,iter_in.t-1),:)> rand(N_learn,1)*ones(1,size(policy.pmat_cum_f{1},2));
        iter_in.micro_state(learn,iter_in.t) = int16(size(policy.pmat_cum_f{1},2) + 1 - sum(trans_rands,2)); % drawn new micro states
        iter_in.cum_meets(learn,iter_in.t)   = policy.pmat_to_meets_succs(iter_in.micro_state(learn ,iter_in.t),2) - 1; % trials, new state, matrix for all firms (iter_in.t+1)
        iter_in.cum_succ(learn,iter_in.t)    = policy.pmat_to_meets_succs(iter_in.micro_state(learn ,iter_in.t),3) - 1; % successes, new state, matrix for all firms starting (t+1)
    end

    stay(learn)  = iter_in.micro_state(learn,iter_in.t-1) - iter_in.micro_state(learn,1) ~= 0; % wasn't in initial state last period
end
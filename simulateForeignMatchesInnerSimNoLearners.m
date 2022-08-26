function [stay, iter_in] = simulateForeignMatchesInnerSimNoLearners(no_learn, stay, mm, iter_in, policy)
  
    % Update cumulative meetings and successes for firms with >= mm.n_size meets, if
    % any. These firms are presumed to know their thetas, hence nn1 = floor(succ_prob*nn2)
    % These calculations don't involve the Q matrix or the associated pmat trans. probs.
    N_no_learn  = sum(no_learn);
    if N_no_learn >0
        % no exog. death, and shipments prev. year
        stay(no_learn) = (rand(N_no_learn,1) > 1-exp(-mm.firm_death_haz));

        %(succ,trial,common succ rate (defunct), network size, prod of firm, macro shock)
        iter_in.cum_meets(no_learn,iter_in.t) = (iter_in.cum_meets(no_learn,iter_in.t-1) + poissinv(rand(N_no_learn,1),mm.theta2(mm.pt_type(iter_in.pt_ndx,2)) ...
            * reshape(policy.lambda_f(floor(mm.theta2(mm.pt_type(iter_in.pt_ndx,2))*(mm.n_size+1)),(mm.n_size+1),1,min((mm.net_size+1),max(iter_in.cum_succ(no_learn,iter_in.t-1)+1,mm.n_size+1)),...
            mm.pt_type(iter_in.pt_ndx,1),iter_in.macro_state_f(iter_in.t-1)),N_no_learn,1))).*stay(no_learn); % resets to 0 if stay==0 or new firm this period

        iter_in.cum_succ(no_learn,iter_in.t) = iter_in.cum_succ(no_learn,iter_in.t-1).*stay(no_learn)...  %.*(1-iter_in.new_firm(no_learn)) ...
            + random('bino',iter_in.cum_meets(no_learn,iter_in.t)-iter_in.cum_meets(no_learn,iter_in.t-1).*stay(no_learn),mm.theta2(mm.pt_type(iter_in.pt_ndx,2)));
        % Note: cum_succ includes single shipment matches (with bad z's) that are dropped after 1 period, but not meetings that were rejected after the sample
        % shipment because of buyer-seller incompatibility.
        % cum_succ and cum_meets will always be 0 after a period t reset due to exit (stay==0).
    end

    
end
function [iterH_in] = simulateHomeMatchesInnerSimClientCounts(iterH_in,mm,policy)

%%                   simulateHomeMatchesInnerSimClientCounts

% This script creates a time series on micro states (#current clients,    
% cumulative successes) in the home market for firms of pt_ndx, given the 
% home market macro time series
%
%% gross additions to clients, before drops and deaths between t-1 and t

t = iterH_in.t;
pt_ndx = iterH_in.pt_ndx;

    if  mm.sim_firm_num_by_prod_succ_type(pt_ndx) > 0
        trans_rands = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.nn_h);
        cntr = 0;

        for ii = 1:size(iterH_in.theta1_cntr,1) % looping over firms with different home-thetas
            cntr2 = cntr+iterH_in.theta1_cntr(ii,2);
            iterH_in.ptm_type = find(policy.firm_type_prod_succ_macro(:,2) == iterH_in.macro_state_h(t)...
                & policy.firm_type_prod_succ_macro(:,3) == iterH_in.theta1_cntr(ii,1)...
                & policy.firm_type_prod_succ_macro(:,4) == mm.pt_type(pt_ndx,1),1,'first');
            pmat_cum_ht = policy.pmat_cum_h{iterH_in.ptm_type};
            
            trans_rands(cntr+1:cntr2,:) = pmat_cum_ht(iterH_in.micro_state(cntr+1:cntr+iterH_in.theta1_cntr(ii,2),t-1),:)...
                > rand(iterH_in.theta1_cntr(ii,2),1)*ones(1,mm.nn_h);
            cntr = cntr2;
        end  % end of home-theta type loop 

        iterH_in.micro_state(:,t) = int16(mm.nn_h + 1 - sum(trans_rands,2)); % drawn new micro states
        iterH_in.cum_succ(:,t)    =  iterH_in.micro_state(:,t) - 1; % cumulative successes, new state, matrix for all firms starting (t+1)

    end

    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    iterH_in.add_cli_cnt(:,t) = max(iterH_in.cum_succ(:,t) - iterH_in.cum_succ(:,t-1),0); %
    % max() resets count to 0 for neg. differences to deal with new exporters

    % identify first period in which a new exporter is active. First
    iterH_in.new_firm(:,t) = iterH_in.micro_state(:,t) == 1;

    %% Deal with endogenous and exogenous drops (not in policy.pmat_cum_h)

    % identify z values at which exporters keep current matches from t-1 to t
    iterH_in.keep.cli = policy.c_val_h(:,mm.pt_type(pt_ndx,1),iterH_in.macro_state_h(t-1))' > 0; % = 1 if want to keep type for t
    iterH_in.drop_Zcut = size(mm.Z,1) - sum(iterH_in.keep.cli); % matches dropped at z value <= drop_Zcut
    % count endogenous drops (z too low to continue)
    iterH_in.drop_cnt = sum(iterH_in.lag_cli_zst.*(1-iterH_in.keep.cli),2);

    % draw the number of exogenous deaths of remaining matches between t-1 and t
    ddum = find(iterH_in.cur_cli_cnt(:,t-1)-iterH_in.drop_cnt > 0);
    if sum(ddum)>0
        iterH_in.exog_deaths(ddum,t-1) =...
            random('bino',iterH_in.cur_cli_cnt(ddum,t-1)-iterH_in.drop_cnt(ddum),(1-exp(-mm.delta)));
    end

    %% update current count for new matches, drops, and exogenous deaths
    iterH_in.cur_cli_cnt(:,t) = iterH_in.add_cli_cnt(:,t) + iterH_in.cur_cli_cnt(:,t-1) ...
        - iterH_in.drop_cnt - iterH_in.exog_deaths(:,t-1) ;
    iterH_in.trans_rands = trans_rands;
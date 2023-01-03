function [iterH_in] = simulateHomeMatchesInnerSimClientCounts(iterH_in,mm,policy)

%%                   simulateHomeMatchesInnerSimClientCounts

% This script creates a time series on micro states (#current clients,    
% cumulative successes) in the home market for firms of pt_ndx, given the 
% home market macro time series
%
% It differs from simulateForeign matchesInnserSimUpdateClientCount because,
% without home market learning, all cases can be handled using the policy
% function implied by the home market intensity matrix. Also, firms with
% the same foreign theta can have different home market thetas. Finally, we
% don't need to keep track of number of meetings in the home market, only
% successes

%% gross additions to clients, before drops and deaths between t-1 and t

t      = iterH_in.t;
pt_ndx = iterH_in.pt_ndx;

    if  mm.sim_firm_num_by_prod_succ_type(pt_ndx) > 0
        trans_rands = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.max_match_h);
        cntr = 0;
        search_inten = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1); %holds search intensity, in case we need it
        theta_record = size(search_inten);

        % load policy matrices and draw states, looping over firms with different home-thetas
        for ii = 1:size(iterH_in.theta1_cntr,1) 
            cntr2 = cntr+iterH_in.theta1_cntr(ii,2);
            iterH_in.ptm_type = find(policy.firm_type_macro_succ_prod(:,2) == iterH_in.macro_state_h(t)...
                & policy.firm_type_macro_succ_prod(:,3) == iterH_in.theta1_cntr(ii,1)...
                & policy.firm_type_macro_succ_prod(:,4) == mm.pt_type(pt_ndx,1),1,'first');
            pmat_cum_ht = policy.pmat_cum_h{iterH_in.ptm_type};
            % draw next period states for given thetaH type
            trans_rands(cntr+1:cntr2,:) = pmat_cum_ht(iterH_in.micro_state(cntr+1:cntr+iterH_in.theta1_cntr(ii,2),t-1),:)...
                > rand(iterH_in.theta1_cntr(ii,2),1)*ones(1,mm.max_match_h);
            search_inten(cntr+1:cntr2,1) = policy.lambda_h(1,iterH_in.theta1_cntr(ii,1),min(iterH_in.cum_succ(cntr+1:cntr+iterH_in.theta1_cntr(ii,2),t-1)+1,mm.net_size+1),mm.pt_type(pt_ndx,1),iterH_in.macro_state_h(t-1));
            theta_record(cntr+1:cntr2,1) = iterH_in.theta1_cntr(ii,1);
            cntr = cntr2;
        end 

        % use trans_rand to infer new micro states (cumulative successes)
        iterH_in.micro_state(:,t) = int16(mm.max_match_h + 1 - sum(trans_rands,2)); 
        % cumulative successes, new state, matrix for all firms starting (t+1)
        iterH_in.cum_succ(:,t)    =  iterH_in.micro_state(:,t) - 1; 

    end

    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    % calculate number of number of new matches
    % how we do this depends on whether the matches are maxed out
    % until that point, change in cum successes are the same as new clients
    % once cum successes are maxed out, we need to use the policy function
    
    maxed_out = (iterH_in.cum_succ(:,t) >= mm.max_match_h-1) & (iterH_in.cum_succ(:,t-1) >= mm.max_match_h-1); 
    not_maxed = ~maxed_out;
    iterH_in.add_cli_cnt(not_maxed,t) = max(iterH_in.cum_succ(not_maxed,t) - iterH_in.cum_succ(not_maxed,t-1),0); %
    % max() resets count to 0 for neg. differences to deal with new exporters
    
    if any(maxed_out)
        maxed_meets = poissrnd(search_inten(maxed_out));
        maxed_succs = binornd(maxed_meets,mm.theta1(theta_record(maxed_out))');
        iterH_in.add_cli_cnt(maxed_out,t) = maxed_succs;
        %iterH_in.add_cli_cnt(maxed_out,t) = min((mm.max_match_h-1) - iterH_in.cur_cli_cnt(maxed_out,t-1), iterH_in.add_cli_cnt(maxed_out,t)); %keep client count under limit
    end

    % identify periods in which a new firm is active in the home mkt., 
    % but hasn't made matches yet 
    iterH_in.new_firm(:,t) = iterH_in.micro_state(:,t) == 1;
    % also identify firms which died and must start over from the next
    % period
    firm_dies = (rand(size(iterH_in.cum_succ(:,t))) < 1-exp(-mm.firm_death_haz));
    iterH_in.cum_succ(firm_dies,t) = 0;
    iterH_in.add_cli_cnt(firm_dies,t) = 0;
    iterH_in.new_firm(firm_dies,t) = 1;

    %% Deal with endogenous and exogenous drops (not in policy.pmat_cum_h)

    % identify z values at which exporters keep current matches from t-1 to t
    iterH_in.keep_cli = policy.c_val_h(:,mm.pt_type(pt_ndx,1),iterH_in.macro_state_h(t-1))' > 0; % = 1 if want to keep type for t
    % matches dropped at z value <= drop_Zcut:
    iterH_in.drop_Zcut = size(mm.Z,1) - sum(iterH_in.keep_cli);
    
%     fprintf('\r\n t=%4.0f, drop_Zcut =%2.0f\r',[t,iterH_in.drop_Zcut]);
    
    % count endogenous drops (z too low to continue)
    iterH_in.drop_cnt = sum(iterH_in.lag_cli_zst.*(1-iterH_in.keep_cli_lag),2);
   
   %% DIAGNOSTIC PRINT: Display Zcuts for firms types with interior values
%      Zcut_lag_temp = size(mm.Z,1) - sum(iterH_in.keep_cli_lag);
%      if Zcut_lag_temp > 0 && Zcut_lag_temp < 15
%        fprintf('\r\n pt_ndx=%4.0f, t=%4.0f, Zcut_lagH =%2.0f',[iterH_in.pt_ndx,t,Zcut_lag_temp]);
%      end
    %% 
     
    % draw the number of exogenous deaths of remaining matches between t-1 and t
    ddum = find(iterH_in.cur_cli_cnt(:,t-1)-iterH_in.drop_cnt > 0);
    if sum(ddum)>0
        iterH_in.exog_deaths(ddum,t-1) =...
            random('bino',iterH_in.cur_cli_cnt(ddum,t-1)-iterH_in.drop_cnt(ddum),(1-exp(-mm.delta)));
    end
   
%%  update current count for new matches, drops, and exogenous deaths 
    
        iterH_in.cur_cli_cnt(:,t) = iterH_in.add_cli_cnt(:,t) + (iterH_in.cur_cli_cnt(:,t-1) ...
        - iterH_in.drop_cnt - iterH_in.exog_deaths(:,t-1)  );
    
       
       iterH_in.cur_cli_cnt(:,t) =...
           (1 - iterH_in.new_firm(:,t).*(1-iterH_in.new_firm(:,t-1)))...
            .*iterH_in.cur_cli_cnt(:,t) ; % to reset client counts when exits occur
     
    
    iterH_in.trans_rands = trans_rands;

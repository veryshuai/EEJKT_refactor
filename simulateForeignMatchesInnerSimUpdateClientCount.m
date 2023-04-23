function [iter_in, drop_Zcut] = simulateForeignMatchesInnerSimUpdateClientCount(iter_in, mm, policy)
    
% if iter_in.t >= 530
%     'pause in simulateForeignMatchesInnerSimUpdateClientCount'
% end

    no_learn = iter_in.cum_meets(:,iter_in.t-1) >= mm.n_size-3; % for picking off seasoned exporters (up to 3 trials before maximum learning)
    learn    = iter_in.cum_meets(:,iter_in.t-1) <  mm.n_size-3; % for picking off exporters who are still learning
    %stay     = ones(mm.sim_firm_num_by_prod_succ_type(iter_in.pt_ndx),1);  % will flag no-learning firms that continue, t-1 to t
    stay = (rand(size(iter_in.cum_meets(:,iter_in.t-1))) > 1-exp(-mm.firm_death_haz));
    

%    [stay, iter_in] = simulateForeignMatchesInnerSimLearners(learn, iter_in, policy, stay);
%    [stay, iter_in] = simulateForeignMatchesInnerSimNoLearners(no_learn, stay, mm, iter_in, policy);
    [iter_in] = simulateForeignMatchesInnerSimLearners(learn, iter_in, policy, stay);
    [iter_in] = simulateForeignMatchesInnerSimNoLearners(no_learn, stay, mm, iter_in, policy);
    
    iter_in.add_cli_cnt(:,iter_in.t) = max(iter_in.cum_succ(:,iter_in.t) - iter_in.cum_succ(:,iter_in.t-1),0); 
    %max() resets count to 0 for neg. differences to deal with new exporter

    % first period for new exporter
%   iter_in.exit_firm(:,iter_in.t-1) = iter_in.cum_meets(:,iter_in.t) - iter_in.cum_meets(:,iter_in.t-1) < 0 ; 
    % = 1 if t-1 is last period before exit
    %DJ: This way of doing exit firm should preserve the idea that a firm
    %only exits if it has ever had positive meetings.
    [iter_in, drop_Zcut, drop_cnt] = simulateForeignMatchesInnerSimDrops(iter_in, policy, mm);
    iter_in.cur_cli_cnt(:,iter_in.t) = iter_in.add_cli_cnt(:,iter_in.t) + iter_in.cur_cli_cnt(:,iter_in.t-1) ...
        - drop_cnt - iter_in.exog_deaths(:,iter_in.t-1) ;
    
    % iter_in.cur_cli_cnt(:,iter_in.t) =...
    %    (1-iter_in.new_firm(:,iter_in.t-1)))...
    %    .*iter_in.cur_cli_cnt(:,iter_in.t) ; % to reset client counts when exits occur 

    %if death, reset everything (should already be correct for non
    %learners, but learners need to be set to zero)
    
  % exog. firm death can't apply to firms adding clients--this messes   
  % up reset of match records in seas_tran
    firm_died = ~stay & iter_in.add_cli_cnt(:,iter_in.t) == 0;
%   firm_died = ~stay;
    iter_in.new_firm(:,iter_in.t)    = firm_died ; 
    iter_in.exit_firm(:,iter_in.t-1) = firm_died & iter_in.cum_meets(:,iter_in.t-1) ~= 0;  
    iter_in.cur_cli_cnt(firm_died,iter_in.t) = 0;
    iter_in.cum_meets(firm_died,iter_in.t) = 0;
    iter_in.cum_succ(firm_died,iter_in.t) = 0;

end



function iter_out = simulateForeignMatches(pt_ndx,macro_state_f,mm,policy)

[iter_in, iter_out] = simulateForeignInnerInitialize(mm, pt_ndx, macro_state_f);

for t = 2:1:mm.periods

    iter_in.t = t;
    %NOTE TO SELF: THIS IS MORE CLEARLY DONE WITH MOD FUNCTION
    if abs(floor((iter_in.t-1)/mm.pd_per_yr) - (iter_in.t-1)/mm.pd_per_yr) <= 0.0001
        iter_in.season = 1; % reset season when previous period completes a year
    end

    % update year, productivity type, and transition probability matrix
    iter_in.year = floor((iter_in.t-1)/mm.pd_per_yr);
    iter_in.mic_type = find(policy.firm_type_prod_succ_macro(:,3) == mm.pt_type(pt_ndx,2) & policy.firm_type_prod_succ_macro(:,4) == mm.pt_type(pt_ndx,1),1,'first');
    type     = find(policy.firm_type_prod_succ_macro(:,2) == macro_state_f(iter_in.t) & policy.firm_type_prod_succ_macro(:,3) == mm.pt_type(pt_ndx,2) & policy.firm_type_prod_succ_macro(:,4) == mm.pt_type(pt_ndx,1),1,'first');
    iter_in.pmat_cum_t = policy.pmat_cum_f{type}; % holds cumulative transition probs across (#success, #meeting) pairs.

    iter_in = simulateForeignMatchesInnerSim(iter_in,mm,policy);

    if iter_in.season == mm.pd_per_yr
        [iter_in,iter_out] = simulateForeignInnerAnnualize(iter_in,iter_out,mm);
    end
    iter_in.season = iter_in.season + 1;

    iter_in.lag_cli_zst  = iter_in.cur_cli_zst;
    iter_in.new_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iter_in.die_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iter_in.trans_zst    = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iter_in.trans_count  = zeros(size(mm.Z,1)+1,size(mm.Z,1)+1,mm.sim_firm_num_by_prod_succ_type(pt_ndx));

end     
end

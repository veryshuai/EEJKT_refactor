%%                   simulateHomeMatchesInnerSimMatchLevelData 

%  First load season to season transitions into mat_tran, which describes
    %  matches of all mm.sim_firm_num_by_prod_succ_type(pt_ndx) of a particular type for a particular transition (t-1 to t).
    mat_tran_all_zeros = ~any(iterH_in.trans_count(:));
    if mat_tran_all_zeros
        mat_tran = zeros(0,4);ship_cur = zeros(0,1); age_vec = zeros(0,1);
    else
        mkt = 2; % =2 for domestic market
 %      [mat_tran,ship_cur,age_vec] = match_sales(mkt,mm,iterH_in.trans_count,age,pt_ndx,macro_state_h(t));
       
    trans_count = iterH_in.trans_count;
        [mat_tran,ship_cur,age_vec] = simulateMatchesInnerSimMatchSales(mkt,mm,iterH_in,age);

    end
    % mat_tran:  [initial state, exporter id, ending state, match revenue]

    if iterH_in.season == 1
        iterH_in.N_match = size(mat_tran,1);
    end

    iterH_in.seas_tran{1,iterH_in.season} = [[t,iterH_in.season,iterH_in.year].*ones(size(mat_tran,1),1),mat_tran,ship_cur,age_vec];
    iterH_in.seas_Zcut(iterH_in.season)   = iterH_in.drop_Zcut;
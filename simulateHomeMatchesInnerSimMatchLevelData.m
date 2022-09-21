function [iterH_in] = simulateHomeMatchesInnerSimMatchLevelData(iterH_in, mm)

%%                   simulateHomeMatchesInnerSimMatchLevelData 

t =iterH_in.t;

%  First load season to season transitions into mat_tran, which describes
    %  matches of all mm.sim_firm_num_by_prod_succ_type(pt_ndx) of a particular type for a particular transition (t-1 to t).
    mat_tran_all_zeros = ~any(iterH_in.trans_count(:));
    if mat_tran_all_zeros
        mat_tran = zeros(0,4);ship_cur = zeros(0,1); age_vec = zeros(0,1);
    else
        mm.mkt = 2; % =2 for domestic market
         [mat_tran,ship_cur,age_vec] = simulateMatchesInnerSimMatchSales(mm.mkt,mm,iterH_in,iterH_in.age);

   iterH_in.ship_cur = ship_cur;
   iterH_in.age_vec = age_vec;

    end
    % mat_tran:  [initial state, exporter id, ending state, match revenue]

    if iterH_in.season == 1
        iterH_in.N_match = size(mat_tran,1);
    end

    iterH_in.seas_tran{1,iterH_in.season} = [[t,iterH_in.season,iterH_in.year].*ones(size(mat_tran,1),1),mat_tran,ship_cur,age_vec];
    iterH_in.seas_Zcut(iterH_in.season)   = iterH_in.drop_Zcut;
    iterH_in.mat_tran = mat_tran;
    iterH_in.mkt = 2;
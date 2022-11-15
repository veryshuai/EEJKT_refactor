%%                   simulateForeignMatchesInnerSimMatchLevelData 

% This function appends sales, shipments and age variables to match records
%%
function [iter_in, mat_tran] = simulateForeignMatchesInnerSimMatchLevelData(iter_in, mm, age, drop_Zcut)

    mat_tran_all_zeros = ~any(iter_in.trans_count(:));
    if mat_tran_all_zeros
        mat_tran = zeros(0,4);ship_cur = zeros(0,1); age_vec = zeros(0,1);
    else
%         macro_shk = iter_in.macro_state_f(iter_in.t);
        mm.mkt = 1; % =1 for foreign market
        [mat_tran,ship_cur,age_vec] = simulateMatchesInnerSimMatchSales(mm.mkt,mm,iter_in,age);
        % mat_tran:  [initial state, exporter id, ending state, match revenue]
        % ship_cur:   match's number of shipments within the current period
        % age_vec:    firm age (# periods)
    end

    if iter_in.season == 1
        iter_in.N_match = size(mat_tran,1);
    end
    
%      if iter_in.year >= 16 && iter_in.season == 12
%           'pause here'
%       end

    iter_in.seas_tran{1,iter_in.season} = [[iter_in.t,iter_in.season,iter_in.year].*ones(size(mat_tran,1),1),mat_tran,ship_cur,age_vec];
    % iter_in.seas_tran: [t, season, year, initial state, exporter id, ending state, match revenue,
    %             # shipments, firm age (# periods)]
    iter_in.seas_Zcut(iter_in.season)   = drop_Zcut;
    
    iter_in.mkt = 1;

end
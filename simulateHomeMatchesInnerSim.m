function iter_out = simulateHomeMatchesInnerSim(iter_out, mm, iterH_in, pt_ndx, policy)

tic

theta1_cntr = iterH_in.theta1_cntr;
macro_state_h = iterH_in.macro_state_h;
firm_yr_sales_lag = iterH_in.firm_yr_sales_lag;

for t = 2:1:mm.periods
    iterH_in.t = t;
    if mod(iterH_in.t-1,mm.pd_per_yr) == 0
        iterH_in.season = 1; % reset season when previous period completes a year
    end

    iterH_in.year = floor((iterH_in.t-1)/mm.pd_per_yr);

    simulateHomeMatchesInnerSimClientCounts
    simulateHomeMatchesInnerSimUpdZHotel
    simulateHomeMatchesInnerSimKickDormant
    simulateHomeMatchesInnerSimFirmAge

    simulateHomeMatchesInnerSimMatchLevelData
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Construct period-specific variables

%     %  First load season to season transitions into mat_tran, which describes
%     %  matches of all mm.sim_firm_num_by_prod_succ_type(pt_ndx) of a particular type for a particular transition (t-1 to t).
%     mat_tran_all_zeros = ~any(iterH_in.trans_count(:));
%     if mat_tran_all_zeros
%         mat_tran = zeros(0,4);ship_cur = zeros(0,1); age_vec = zeros(0,1);
%     else
%         mkt = 2; % =2 for domestic market
%         [mat_tran,ship_cur,age_vec] = match_sales(mkt,mm,iterH_in.trans_count,age,pt_ndx,macro_state_h(t));
%   %     [mat_tran,ship_cur,age_vec] = simulateMatchesInnerSimMatchSales(mkt,mm,iter_in,trans_count,age);
% 
%     end
%     % mat_tran:  [initial state, exporter id, ending state, match revenue]
% 
%     if iterH_in.season == 1
%         iterH_in.N_match = size(mat_tran,1);
%     end
% 
%     iterH_in.seas_tran{1,iterH_in.season} = [[t,iterH_in.season,iterH_in.year].*ones(size(mat_tran,1),1),mat_tran,ship_cur,age_vec];
%     iterH_in.seas_Zcut(iterH_in.season)   = drop_Zcut;

    iterH_in.N_match = iterH_in.N_match + size(mat_tran,1);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% construct annualized variables
    if iterH_in.season == mm.pd_per_yr

        [~,firm_yr_sales] =...
            season_merge(iterH_in.seas_tran,iterH_in.N_match,mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.pd_per_yr);

        % firm_yr_sales:[firm ID, total dom. sales, total dom. shipments, firm age in domestic market]

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% the following matrices accumulate annualized values over time and firm types
        theta_h_firm = iterH_in.theta_h(firm_yr_sales(:,1));
        ttt = ones(size(firm_yr_sales,1),1).*[t,ptm_type];
        iter_out.firm_h_yr_sales = [iter_out.firm_h_yr_sales;[ttt,firm_yr_sales]];
        iter_out.theta_h_firm  = [iter_out.theta_h_firm;theta_h_firm]; % keep track of domestic thetas for each firm
        % iter_out.firm_h_yr_sales: [t,type,firm ID, total sales, total shipments,firm age]

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         seas_Zcut = zeros(1,mm.pd_per_yr);

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Construct and cumulate moments
        if iterH_in.year > mm.burn  % cosntruct moments for firm domestic sales regressions
            % autoregressions and degree distribution

            [x,y,fmoms_xx,fmoms_xy,fysum,fn_obs] = firm_reg_h_moms(firm_yr_sales,firm_yr_sales_lag,mm.sim_firm_num_by_prod_succ_type(pt_ndx));

            iter_out.x_fsales_h   = [iter_out.x_fsales_h;x];
            iter_out.y_fsales_h   = [iter_out.y_fsales_h;y];
            iter_out.fmoms_h_xx = iter_out.fmoms_h_xx + fmoms_xx; % cumulate moments for home sales AR1
            iter_out.fmoms_h_xy = iter_out.fmoms_h_xy + fmoms_xy; % cumulate moments for home sales AR1
            iter_out.fysum_h    = iter_out.fysum_h + fysum;
            iter_out.fnobs_h    = iter_out.fnobs_h + fn_obs ;

        end   % year > mm.burn if statement
        firm_yr_sales_lag = firm_yr_sales; % stack data for firm regression

    end   % season == mm.pd_per_yr if statement

    iterH_in.season = iterH_in.season + 1;

    %% load lagged client state matrix and re-initialize objects

    iterH_in.lag_cli_zst  = iterH_in.cur_cli_zst;
    iterH_in.new_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iterH_in.die_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iterH_in.trans_zst    = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iterH_in.trans_count  = zeros(size(mm.Z,1)+1,size(mm.Z,1)+1,mm.sim_firm_num_by_prod_succ_type(pt_ndx));

if t==mm.periods
    'pause here'
end

end
end
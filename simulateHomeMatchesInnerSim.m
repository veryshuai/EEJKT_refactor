function iter_out = simulateHomeMatchesInnerSim(iter_out, mm, iterH_in, pt_ndx, policy)

tic
iterH_in.pt_ndx = pt_ndx;

for t = 2:1:mm.periods
    iterH_in.t = t;
    if mod(iterH_in.t-1,mm.pd_per_yr) == 0
        iterH_in.season = 1; % reset season when previous period completes a year
    end

    iterH_in.year = floor((iterH_in.t-1)/mm.pd_per_yr);

[iterH_in] = simulateHomeMatchesInnerSimClientCounts(iterH_in, mm, policy);
[iterH_in] = simulateHomeMatchesInnerSimUpdZHotel(iterH_in, mm, policy);
[iterH_in] = simulateHomeMatchesInnerSimKickDormant(iterH_in, mm);
[iterH_in] = simulateHomeMatchesInnerSimFirmAge(iterH_in, mm);
[iterH_in] = simulateHomeMatchesInnerSimMatchLevelData(iterH_in, mm);

 mat_tran = iterH_in.mat_tran;
 iterH_in.N_match = iterH_in.N_match + size(mat_tran,1);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% construct annualized variables
    if iterH_in.season == mm.pd_per_yr

  [iterH_in.mat_h_yr_sales,iterH_in.firm_h_yr_sales] = season_merge(iterH_in,mm);

     
        % firm_h_yr_sales:[firm ID, total dom. sales, total dom. shipments, firm age in domestic market]

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% the following matrices accumulate annualized values over time and firm types

        theta_h_firm = iterH_in.theta_h(floor(iterH_in.firm_h_yr_sales(:,1)));
        ttt = ones(size(iterH_in.firm_h_yr_sales,1),1).*[t,iterH_in.ptm_type];
        iter_out.firm_h_yr_sales = [iter_out.firm_h_yr_sales;[ttt,iterH_in.firm_h_yr_sales]];
        iter_out.theta_h_firm  = [iter_out.theta_h_firm;theta_h_firm]; % keep track of domestic thetas for each firm
        % iter_out.firm_h_yr_sales: [t,type,firm ID, total sales, total shipments,firm age]

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Construct and cumulate moments
           
        
        if iterH_in.year > mm.burn  % cosntruct moments for firm domestic sales regressions
            % autoregressions and degree distribution
            
        [iterH_in.mat_h_cont_2yr,iterH_in.mat_h_yr_sales,iterH_in.mat_h_yr_sales_lag,~] = ...
            mat_yr_splice_v2(iterH_in.mat_h_yr_sales,iterH_in.mat_h_yr_sales_lag,mm,iterH_in.year);


            [x,y,fmoms_xx,fmoms_xy,fysum,fn_obs] = firm_reg_h_moms(iterH_in,mm);

            iter_out.x_fsales_h   = [iter_out.x_fsales_h;x];
            iter_out.y_fsales_h   = [iter_out.y_fsales_h;y];
            iter_out.fmoms_h_xx = iter_out.fmoms_h_xx + fmoms_xx; % cumulate moments for home sales AR1
            iter_out.fmoms_h_xy = iter_out.fmoms_h_xy + fmoms_xy; % cumulate moments for home sales AR1
            iter_out.fysum_h    = iter_out.fysum_h + fysum;
            iter_out.fnobs_h    = iter_out.fnobs_h + fn_obs ;

        end   % year > mm.burn if statement
        iterH_in.firm_h_yr_sales_lag = iterH_in.firm_h_yr_sales; % stack data for firm regression
        iterH_in.mat_h_yr_sales_lag = iterH_in.mat_h_yr_sales;   
    end   % season == mm.pd_per_yr if statement

    iterH_in.season = iterH_in.season + 1;

    %% load lagged client state matrix and re-initialize objects

    iterH_in.lag_cli_zst  = iterH_in.cur_cli_zst;
    iterH_in.new_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iterH_in.die_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iterH_in.trans_zst    = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iterH_in.trans_count  = zeros(size(mm.Z,1)+1,size(mm.Z,1)+1,mm.sim_firm_num_by_prod_succ_type(pt_ndx));

% if t==mm.periods
%     'pause here'
% end

end
end
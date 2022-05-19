function [iter_in,iter_out] = simulateForeignMatchesInnerAnnualize(iter_in,iter_out,mm)

        [mat_yr_sales,firm_yr_sales] =...
            season_merge(iter_in.seas_tran,iter_in.N_match,mm.sim_firm_num_by_prod_succ_type(iter_in.pt_ndx),mm.pd_per_yr);
        % mat_yr_sales:  [firm ID, match-specific sales, shipments, boy Z, eoy Z,
        %                 match age in periods (w/in year), firm age in periods]
        % firm_yr_sales: [firmID,sales,#shipments,firm age]

        % convert age in periods to age in years for first-year observations.
        % After the first year, conversion to years handled by mat_yr_splice
        if iter_in.year==1
            mat_yr_sales(:,6)  = mat_yr_sales(:,6) > 0;  % set match age in years to 1 if match age in periods > 0
            mat_yr_sales(:,7)  = mat_yr_sales(:,7) > 0;  % set firm age in years to 1 if year age in periods > 0
            firm_yr_sales(:,4) = firm_yr_sales(:,4) > 0; % set firm age in years to 1 if year age in periods > 0
            iter_in.year_lag           = 1;
        end

        % # unsuccessful meetings (duds) over previous year, by firm (needed for degree distribution later)
        yr_tlag = iter_in.t-mm.pd_per_yr;
        cum_duds  = iter_in.cum_meets(:,yr_tlag:iter_in.t) - iter_in.cum_succ(:,yr_tlag:iter_in.t); % previous mm.pd_per_yr + 1 cumulative duds
        curr_duds = cum_duds(:,2:mm.pd_per_yr+1)-(iter_in.new_firm(:,yr_tlag+1:iter_in.t)==0).*cum_duds(:,1:mm.pd_per_yr);

        iter_out.singletons = sum(sum(curr_duds));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%  Generate observations on intervals between meetings, cumulative meetings, and cumulative succcesses

        if iter_in.t>3*mm.pd_per_yr
            [time_gap,iter_in.mkt_exit] = time_gaps(iter_in.t,iter_in.exit_firm,mm.pd_per_yr,iter_in.cum_meets,iter_in.cum_succ);
            %           time_gap: (1) firm_ID, (2) periods into interval, (3) time gap,
            %           (4) # new meetings, (5) iter_in.t, (6) cum. meetings, (7) cum succeses

            iter_out.time_gaps = [iter_out.time_gaps;time_gap];
        end

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if iter_in.year > 1  % construct moments for match regressions, firm export regressions

            if size(mat_yr_sales,1)*size(iter_in.mat_yr_sales_lag)>0
                % drops types without sales in both current and lagged years

                % mat_yr_splice splices consecutive obs. on annualized data and updates
                % match ages and firm ages, converting them to years. It's really costly!

                ncols = size(mat_yr_sales,2);

                [mat_cont_2yr,mat_yr_sales,iter_in.mat_yr_sales_adj,iter_in.year_lag] =...
                    mat_yr_splice_v2(mat_yr_sales,iter_in.mat_yr_sales_lag,mm,iter_in.year_lag,iter_in.year);

                %  mat_cont_2yr: [mat_yr_sales_lag(ff_cont_lag,:), mat_yr_sales(ff_cont,:)]
                %  mat_yr_sales: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age in yrs, firm age in yrs]
                %  iter_in.mat_yr_sales_adj: same as lagged mat_yr_sales except eoy Z set to zero if no sales next year

                %% the following matrices accumulate annualized values over time and firm types

                mat_matur_dat =  [mat_yr_sales(:,2),mat_yr_sales(:,4:7),mat_yr_sales(:,1),ones(size(mat_yr_sales,1),1).*iter_in.t/mm.pd_per_yr] ;
                % agg_mat_matur: [sales, boy Z, eoy Z, match age, fmic_typeirm age, firm_ID, yr]

                if iter_in.year > mm.burn  % don't start building simulated data set until burn-in finished

                    tt =  ones(size(mat_yr_sales,1),1).*[iter_in.t,iter_in.mic_type]; % add cols 1 and 2: t, firm type
                    iter_out.mat_yr_sales  = [iter_out.mat_yr_sales;[tt,mat_yr_sales]];
                    % agg_mat_yr_sales: [t,type,firm ID, match sales, shipments, boy Z, eoy Z, w/in yr. match age, firm age]

                    if iter_in.year_lag == iter_in.year % check that mat_yr_splice_v2 ran & updated iter_in.mat_yr_sales_adj
                        tt2 =  ones(size(iter_in.mat_yr_sales_adj,1),1).*[iter_in.t-mm.pd_per_yr,iter_in.mic_type];
                        iter_out.mat_yr_sales_adj = [iter_out.mat_yr_sales_adj;[tt2,iter_in.mat_yr_sales_adj]]; % add cols 1 and 2: t, firm type
                        % agg_mat_yr_sales_adj: [t,type,firm ID, match sales, shipments, boy Z, adj_eoy Z, w/in yr. match age, firm age]
                        ttt = ones(size(firm_yr_sales,1),1).*[iter_in.t,iter_in.mic_type];
                        iter_out.firm_f_yr_sales = [iter_out.firm_f_yr_sales;[ttt,firm_yr_sales]]; % add cols 1 and 2: t, firm type
                        % agg_firm_yr_sales: [t,type,firm ID, total exports,total shipments,firm age]
                    end

                    iter_out.mat_matur =  [iter_out.mat_matur; mat_matur_dat];

                    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %  match regression moments

                    [mat_ar1_x,mat_ar1_y,moms_xx,moms_xy,ysum,n_obs] = match_reg_moms(mat_cont_2yr,ncols);

                    iter_out.moms_xx = iter_out.moms_xx + moms_xx; % cumulate moments for match regression
                    iter_out.moms_xy = iter_out.moms_xy + moms_xy; % cumulate moments for match regression
                    iter_out.ysum    = iter_out.ysum + ysum;
                    iter_out.nobs    = iter_out.nobs + n_obs;

                    iter_out.mat_ar1_x = [iter_out.mat_ar1_x;mat_ar1_x];
                    iter_out.mat_ar1_y = [iter_out.mat_ar1_y;mat_ar1_y];
                end

            end

            if iter_in.year >= mm.burn  % don't start building simulated data set until burn-in finished

                % autoregressions and degree distribution

                [fmoms_xx,fmoms_xy,fysum,fn_obs] = firm_reg_moms(firm_yr_sales,iter_in.firm_yr_sales_lag,mm.sim_firm_num_by_prod_succ_type(iter_in.pt_ndx));

                iter_out.fmoms_xx = iter_out.fmoms_xx + fmoms_xx; % cumulate moments for firm regression
                iter_out.fmoms_xy = iter_out.fmoms_xy + fmoms_xy; % cumulate moments for firm regression
                iter_out.fysum    = iter_out.fysum + fysum;
                iter_out.fnobs    = iter_out.fnobs + fn_obs ;

                % foreign market exit regression moments

                ff_exit = iter_in.mkt_exit(:,2)>0;  % column 2 of mkt_exit is number of meetings
                if sum(ff_exit,1)>0
                    [exit_moms_xx,exit_moms_xy,sum_succ_rate,sum_exits,exit_obs] = mkt_exit_moms(iter_in.mkt_exit);

                    iter_out.exit_xx = iter_out.exit_xx + exit_moms_xx;
                    iter_out.exit_xy = iter_out.exit_xy + exit_moms_xy';
                    iter_out.sum_succ_rate = iter_out.sum_succ_rate + sum_succ_rate;
                    iter_out.exit_obs = iter_out.exit_obs + exit_obs;
                    iter_out.sum_exits = iter_out.sum_exits + sum_exits;
                end
                % match exit regression moments
                if iter_in.year_lag == iter_in.year
                    ff_mexit = iter_in.mat_yr_sales_adj(:,2)>0;
                    if sum(ff_mexit,1)>0 % positive sales for at least one match
                        [mat_exit_x,mat_exit_y,mat_exit_moms_xx,mat_exit_moms_xy,mat_obs,nmat_exit]...
                            = match_exit_moms(iter_in.mat_yr_sales_adj(ff_mexit,:),mm.pd_per_yr);

                        % Notes on variables:

                        %         mat_exit_y  = matches(ff,5)==0;         % match dead by end of year
                        %         mat_exit_x = [x0,x1,x2,x3,x4];
                        %           x0 = ones(size(ff,1),1);              % ff picks off non-missing obs.
                        %           x1 = matches(ff,4)==0;                % first year dummy
                        %           x2 = log(matches(ff,2));              % sales during year
                        %           x3 = log(1+matches(ff,6)./mm.pd_per_yr); % age of match
                        %           x4 = log(1+matches(ff,7));            % age of exporter

                        iter_out.mat_exit_moms_xx = iter_out.mat_exit_moms_xx + mat_exit_moms_xx;
                        iter_out.mat_exit_moms_xy = iter_out.mat_exit_moms_xy + mat_exit_moms_xy;
                        iter_out.mat_obs      = iter_out.mat_obs + mat_obs;
                        iter_out.nmat_exit    = iter_out.nmat_exit + nmat_exit;
                        iter_out.mat_exit_x   = [iter_out.mat_exit_x;mat_exit_x];
                        iter_out.mat_exit_y   = [iter_out.mat_exit_y;mat_exit_y];
                    end
                end
                % shipment and match counter

                [nship_obs,ln_ships,match_count] = match_shpt_cntr(iter_in.mat_yr_sales_adj,mm.max_match);

                iter_out.ship_obs    = iter_out.ship_obs + nship_obs ;
                iter_out.ln_ships    = iter_out.ln_ships + ln_ships ;
                iter_out.match_count = iter_out.match_count + match_count ;
                % include all the matches that generated a single sample shipment:
                iter_out.match_count(1) = iter_out.match_count(1) + iter_out.singletons;

            end   % year > mm.burn  if statement
        end
        iter_in.mat_yr_sales_lag = mat_yr_sales;   % update lags
        iter_in.firm_yr_sales_lag = firm_yr_sales;

end
function [iterH_in,iter_out] = simulateHomeMatchesInnerMoments(iterH_in,iter_out,mm)

        if iterH_in.year > mm.burn  % cosntruct moments for firm domestic sales regressions
            % autoregressions and degree distribution

            [x,y,fmoms_xx,fmoms_xy,fysum,fn_obs] = firm_reg_h_moms(iterH_in.firm_yr_sales,...
                iterH_in.firm_yr_sales_lag,mm.sim_firm_num_by_prod_succ_type(iterH_in.pt_ndx));

            iter_out.x_fsales_h   = [iter_out.x_fsales_h;x];
            iter_out.y_fsales_h   = [iter_out.y_fsales_h;y];
            iter_out.fmoms_h_xx = iter_out.fmoms_h_xx + fmoms_xx; % cumulate moments for home sales AR1
            iter_out.fmoms_h_xy = iter_out.fmoms_h_xy + fmoms_xy; % cumulate moments for home sales AR1
            iter_out.fysum_h    = iter_out.fysum_h + fysum;
            iter_out.fnobs_h    = iter_out.fnobs_h + fn_obs ;

        end   % year > mm.burn if statement
        iterH_in.firm_yr_sales_lag = iterH_in.firm_yr_sales; % stack data for firm regression
        iterH_in.Zcut_eoy_lag      = iterH_in.Zcut_eoy;
end
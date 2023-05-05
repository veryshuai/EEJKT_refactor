function [iterH_in,iter_out] = simulateHomeMatchesInnerAnnualize(iterH_in,iter_out,mm)

%         [~,iterH_in.firm_yr_sales] =...
%             season_merge(iterH_in.seas_tran,iterH_in.N_match,mm.sim_firm_num_by_prod_succ_type(pt_ndx),mm.pd_per_yr);

        [iterH_in.mat_yr_sales,iterH_in.firm_yr_sales] =...
            season_merge(iterH_in.seas_tran,iterH_in.N_match,mm.sim_firm_num_by_prod_succ_type(iterH_in.pt_ndx),mm.pd_per_yr);

        % iterH_in.firm_yr_sales:[firm ID, total dom. sales, total dom. shipments, firm age in domestic market]

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% the following matrices accumulate annualized values over time and firm types
        theta_h_firm = iterH_in.theta_h(iterH_in.firm_yr_sales(:,1));
        ttt = ones(size(iterH_in.firm_yr_sales,1),1).*[iterH_in.t,iterH_in.ptm_type];
        iter_out.firm_h_yr_sales = [iter_out.firm_h_yr_sales;[ttt,iterH_in.firm_yr_sales]];
        iter_out.theta_h_firm  = [iter_out.theta_h_firm;theta_h_firm]; % keep track of domestic thetas for each firm
        % iter_out.firm_h_yr_sales: [t,type,firm ID, total sales, total shipments,firm age]

end
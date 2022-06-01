function iter_out = simulateHomeMatches(pt_ndx,macro_state_h,mm,policy,iter_out)  

% iter_in = struct;

% create home theta draws
theta_df = (size(mm.theta1,2)+1)*ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1) - ...
    sum(ones(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1)*mm.th1_cdf >...
    rand(mm.sim_firm_num_by_prod_succ_type(pt_ndx),1)*ones(1,size(mm.theta1,2)),2);

% list the non-zero values and their frequencies
[uv,~,idx]   = unique(theta_df);
theta1_cntr  = [uv,accumarray(idx(:),1)];
theta_h      = sortrows(theta_df);

[iterH_in, iter_out] = simulateHomeInnerInitialize(mm, pt_ndx, iter_out);

firm_yr_sales_lag = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),4);
% firm_yr_sales_lag will contain: [firmID,sales,#shipments,firm age]
iter_out.abort_flag_h = 0;

iter_out = simulateHomeMatchesInnerSim(iter_out, mm, iterH_in, pt_ndx, theta1_cntr, policy, macro_state_h, theta_h,firm_yr_sales_lag);

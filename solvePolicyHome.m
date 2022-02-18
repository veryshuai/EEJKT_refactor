function [V,l_opt] = solvePolicyHome(policy,mm)

V = zeros(size(mm.Q0_h,1),1,mm.dim1,mm.net_size+1); %[macro-prod states, theta_g,theta_h, max network effect]
V = zeros(1,mm.dim1,mm.net_size+1,2*mm.phi_size+1,2*mm.x_size+1);
%lambda_h (common succ rate (defunct), known succ rate, network size, prod of firm, H macro shock)
l_opt = zeros(size(V)); %optimal search intensities

times_nostop = 0;
for k = 1:mm.dim1
    l_opt(:,1,k,mm.net_size+1) = mm.l_opt_func_h(mm.theta1(k),mm.net_size,policy.pi_h,0,0,0);
    den = mm.r+abs(diag(mm.Q0_h));
    perm_flow = -mm.cost_h(l_opt(:,1,k,mm.net_size+1),mm.net_size+1)+mm.theta1(k)*l_opt(:,1,k,mm.net_size+1).*policy.pi_h;
    Q0_diag = -mm.Q0_h_d;
    Q0_diag(1:size(mm.Q0_h_d,1)+1:end) = den;
    V(:,1,k,mm.net_size+1) = max(Q0_diag^-1*perm_flow,0);

    for l = 1:mm.net_size

        init = V(:,1,k,mm.net_size-l+2);
        [new_v,nostop,differ] = val_func_iter(init,mm.b,mm.theta1(k),policy.pi_h,mm.net_size+1-l,mm.Q0_h,mm.r,diag_Q,mm.Q0_h_d,V(:,1,k,mm.net_size-l+2),1,mm.gam,scl,mm.cs_h,mm.v_tolerance,mm.v_rel_tolerance,1,[],[],[],c,mm.l_opt_func_h,policy.pi_h);
        times_nostop = times_nostop + nostop;

        if differ/max(policy.pi_h) > worst_scaled_error(3)
            worst_scaled_error = [differ,max(policy.pi_h),differ/max(policy.pi_h)];
        end
        if differ > worst_error(1)
            worst_error = [differ,max(policy.pi_h),differ/max(policy.pi_h)];
        end

        V(:,1,k,mm.net_size-l+1) = new_v;
        [~,~,l_opt(:,1,k,mm.net_size+1-l)] = sim_solve_h(V(:,1,k,mm.net_size+1-l),mm.b,mm.theta1(k),policy.pi_h,mm.net_size+1-l,size(mm.Q0_h,1),mm.r,diag_Q,mm.Q0_h_d,V(:,1,k,mm.net_size-l+2),mm.gam,scl,mm.cs_h,c,mm.l_opt_func_h);

    end
end
end


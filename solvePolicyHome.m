function [V,l_opt] = solvePolicyHome(policy,mm)

V = zeros(size(mm.Q_px_h,1),1,mm.dim1,mm.net_size+1); %[macro-prod states, theta_g,theta_h, max network effect]
l_opt = zeros(size(V));

for k = 1:mm.dim1

    l_opt(:,1,k,mm.net_size+1) = mm.l_opt_func_h(mm.theta1(k),mm.net_size,policy.pi_h,0,0,0);
    den = mm.r+abs(diag(mm.Q_px_h));
    perm_flow = -mm.cost_h(l_opt(:,1,k,mm.net_size+1),mm.net_size+1)+mm.theta1(k)*l_opt(:,1,k,mm.net_size+1).*policy.pi_h;
    Q_diag = -mm.Q_px_h_d;
    Q_diag(1:size(mm.Q_px_h_d,1)+1:end) = den;
    V(:,1,k,mm.net_size+1) = max(Q_diag^-1*perm_flow,0);

    for N = 1:mm.net_size

        val_guess = V(:,1,k,mm.net_size-N+2);
        new_v = val_func_iter(val_guess,mm.b,mm.theta1(k),policy.pi_h,mm.net_size+1-N,mm.Q_px_h,mm.r,abs(diag(mm.Q_px_h)),mm.Q_px_h_d,V(:,1,k,mm.net_size-N+2),1,mm.gam,1,mm.cs_h,mm.v_tolerance,mm.v_rel_tolerance,1,[],[],[],mm.cost_h,mm.l_opt_func_h,policy.pi_h);
        V(:,1,k,mm.net_size-N+1) = new_v;
        [~,~,l_opt(:,1,k,mm.net_size+1-N)] = iterateValueFunction_h(V(:,1,k,mm.net_size+1-N),mm.theta1(k),policy.pi_h,mm.net_size+1-N,mm.r,abs(diag(mm.Q_px_h)),mm.Q_px_h_d,V(:,1,k,mm.net_size-N+2),mm.cost_h,mm.l_opt_func_h);

    end
end
end


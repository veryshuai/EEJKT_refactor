function [V,l_opt] = solvePolicyHome(policy,mm)

V = zeros(1,mm.dim1,mm.net_size+1,2*mm.phi_size+1,2*mm.x_size+1);
%lambda_h (common succ rate (defunct), known succ rate, network size, prod of firm, H macro shock)
l_opt = zeros(size(V)); 

for prod_ind = 1:size(mm.Phi)
    for k = 1:mm.dim1
        
        l_opt(1,k,mm.net_size+1,prod_ind,:) = mm.l_opt_func_h(mm.theta1(k),mm.net_size,policy.pi_h(prod_ind,:),0,0,0);
        den = mm.r+abs(diag(mm.Q_h));
        perm_flow = -mm.cost_h(squeeze(l_opt(1,k,mm.net_size+1,prod_ind,:)),mm.net_size+1)+mm.theta1(k)*squeeze(l_opt(1,k,mm.net_size+1,prod_ind,:)).*policy.pi_h(prod_ind,:)';
        Q_diag = -mm.Q_h_d;
        Q_diag(1:size(mm.Q_h_d,1)+1:end) = den;
        V(1,k,mm.net_size+1,prod_ind,:) = max(Q_diag^-1*perm_flow,0);

        for N = 1:mm.net_size

            val_guess = squeeze(V(1,k,mm.net_size-N+2,prod_ind,:));
            new_v = val_func_iter(val_guess,mm.b,mm.theta1(k),policy.pi_h(prod_ind,:)',mm.net_size+1-N,mm.Q_h,mm.r,abs(diag(mm.Q_h)),mm.Q_h_d,squeeze(V(1,k,mm.net_size-N+2,prod_ind,:)),1,mm.gam,1,mm.cs_h,mm.v_tolerance,mm.v_rel_tolerance,1,[],[],[],mm.cost_h,mm.l_opt_func_h,policy.pi_h(prod_ind,:)');
            V(1,k,mm.net_size-N+1,prod_ind,:) = new_v;
            [~,~,l_opt(1,k,mm.net_size+1-N,prod_ind,:)] = iterateValueFunction_h(squeeze(V(1,k,mm.net_size+1-N,prod_ind,:)),mm.theta1(k),policy.pi_h(prod_ind,:)',mm.net_size+1-N,mm.r,abs(diag(mm.Q_h)),mm.Q_h_d,squeeze(V(1,k,mm.net_size-N+2,prod_ind,:)),mm.cost_h,mm.l_opt_func_h);

        end
    end
end
end


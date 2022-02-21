function [V,l_opt] = solvePolicyForeign(policy,mm)

V = zeros(mm.n_size+1,mm.n_size+1,1,mm.net_size+1,size(mm.Phi,1),size(mm.X_f,1));
% lambda_f (succ, trial, common succ rate (defunct), network size, prod of firm, F macro shock)
l_opt = zeros(size(V));

for prod_ind = 1:size(mm.Phi,1)
    for succ = 1:mm.n_size+1 

        l_opt(succ,mm.n_size+1,1,mm.net_size+1,prod_ind,:) = mm.l_opt_func_f(policy.postSuccessProb_f(mm.n_size+1,succ),mm.net_size+1,policy.pi_f(prod_ind,:),0,0,0);
        Q0_diag = -mm.Q_f_d;
        Q0_diag(1:size(mm.Q_f_d,1)+1:end) = mm.r+abs(diag(mm.Q_f));
        perm_flow = policy.postSuccessProb_f(mm.n_size+1,succ)*squeeze(l_opt(succ,mm.n_size+1,1,mm.net_size+1,prod_ind,:)).*policy.pi_f(prod_ind,:)'-squeeze(mm.cost_f(l_opt(succ,mm.n_size+1,1,mm.net_size+1,prod_ind,:),mm.net_size+1));
        % next eq. follows from V = (mm.r + abs(diag(Q0)))^-1 * (perm_flow + Q0_d * V)
        V(succ,mm.n_size+1,1,mm.net_size+1,prod_ind,:) = max(Q0_diag^-1*squeeze(perm_flow),0);

        % learning maxed out
        for m = 1:mm.net_size+1-succ

            %value function iteration
            new_v = val_func_iter(squeeze(V(succ,mm.n_size+1,1,mm.net_size+2-m,prod_ind,:)),mm.b,policy.postSuccessProb_f(mm.n_size+1,succ),policy.pi_f(prod_ind,:)',mm.net_size+1-m,mm.Q_f,mm.r,abs(diag(mm.Q_f)),mm.Q_f_d,squeeze(V(succ,mm.n_size+1,1,mm.net_size+2-m,prod_ind,:)),1,mm.gam,1,mm.cs_f,mm.v_tolerance,mm.v_rel_tolerance,1,[],[],[],mm.cost_f,mm.l_opt_func_f,policy.pi_f(prod_ind,:)');
            V(succ,mm.n_size+1,1,mm.net_size+1-m,prod_ind,:) = new_v;

            % get policy function
            [~,~,l_opt(succ,mm.n_size+1,1,mm.net_size+1-m,prod_ind,:)] = iterateValueFunction_h(squeeze(V(succ,mm.n_size+1,1,mm.net_size+1-m,prod_ind,:)),...
                policy.postSuccessProb_f(mm.n_size+1,succ),policy.pi_f(prod_ind,:)',mm.net_size+1-m,mm.r,abs(diag(mm.Q_f)),mm.Q_f_d,squeeze(V(succ,mm.n_size+1,1,mm.net_size+2-m,prod_ind,:)),mm.cost_f,mm.l_opt_func_f);
        end
    end

for tri = 1:mm.n_size 
    for succ = 1:mm.n_size+1-tri

        new_v = val_func_iter(squeeze(V(succ,mm.n_size+2-tri,1,succ,prod_ind,:)),mm.b,policy.postSuccessProb_f,policy.pi_f(prod_ind,:)',mm.n_size+1-tri,mm.Q_f,mm.r,abs(diag(mm.Q_f)),mm.Q_f_d,squeeze(V(succ+1,mm.n_size-tri+2,1,succ+1,prod_ind,:)),squeeze(V(succ,mm.n_size-tri+2,1,succ,prod_ind,:)),mm.gam,1,mm.cs_f,mm.v_tolerance,mm.v_rel_tolerance,0,mm.n_size,tri,succ,mm.cost_f,mm.l_opt_func_f,policy.pi_f(prod_ind,:)');
        V(succ,mm.n_size+1-tri,1,succ,prod_ind,:) = new_v;
        [~,~,l_opt(succ,mm.n_size+1-tri,1,succ,prod_ind,:)] = iterateValueFunction_f(squeeze(V(succ,mm.n_size+1-tri,1,succ,prod_ind,:)),policy.postSuccessProb_f,policy.pi_f(prod_ind,:)',mm.n_size+1-tri,succ,mm.r,abs(diag(mm.Q_f)),mm.Q_f_d,squeeze(V(succ,mm.n_size-tri+2,1,succ,prod_ind,:)),squeeze(V(succ+1,mm.n_size-tri+2,1,succ+1,prod_ind,:)),mm.cost_f,mm.l_opt_func_f);
    end
end
end

end


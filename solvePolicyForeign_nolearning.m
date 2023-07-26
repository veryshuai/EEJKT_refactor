function [V,l_opt] = solvePolicyForeign_nolearning(policy,mm)

V = zeros(mm.n_size+1,mm.n_size+1,1,mm.net_size+1,2*mm.phi_size+1,2*mm.x_size+1);
l_opt = zeros(size(V));

    % lambda_f(succ, trial, common succ rate (defunct), network size, prod of firm, F macro shock) 

for prod_ind = 1:size(mm.Phi)
    for k = 1:mm.dim1
        for succ = 1:mm.n_size+1
    
            l_opt(succ,mm.n_size+1,1,mm.net_size+1,prod_ind,:) = mm.l_opt_func_h(mm.theta1(k),mm.net_size+1,policy.pi_f(prod_ind,:),0,0,0);
            Q0_diag = -mm.Q_f_d;
            Q0_diag(1:size(mm.Q_f_d,1)+1:end) = mm.r+abs(diag(mm.Q_f));
            perm_flow = mm.theta1(k)*squeeze(l_opt(succ,mm.n_size+1,1,mm.net_size+1,prod_ind,:)).*policy.pi_f(prod_ind,:)'-mm.cost_f(squeeze(l_opt(succ,mm.n_size+1,1,mm.net_size+1,prod_ind,:)),mm.net_size+1);
        % next eq. follows from V = (mm.r + abs(diag(Q0)))^-1 * (perm_flow + Q0_d * V)
        V(succ,mm.n_size+1,1,mm.net_size+1,prod_ind,:) = max(Q0_diag^-1*perm_flow,0);
    
        % learning maxed out
        for m = 1:mm.net_size+1-succ
    
            new_v = val_func_iter('no_learning',squeeze(V(succ,mm.n_size+1,1,mm.net_size+2-m,prod_ind,:)),squeeze(V(succ,mm.n_size+1,1,mm.net_size+2-m,prod_ind,:)),[],mm.theta1(k),policy.pi_f(prod_ind,:),mm.net_size+1-m,mm,policy,'foreign');
            V(succ,mm.n_size+1,1,mm.net_size+1-m,prod_ind,:) = new_v;
            [~,l_opt(succ,mm.n_size+1,1,mm.net_size+1-m,prod_ind,:)] = iterationError_h(squeeze(V(succ,mm.n_size+1,1,mm.net_size+1-m,prod_ind,:)),...
                mm.theta1(k),policy.pi_f(prod_ind,:),mm.net_size+1-m,squeeze(V(succ,mm.n_size+1,1,mm.net_size+2-m,prod_ind,:)),mm,policy,'foreign');
    
        end
    end
    
    for tri_iterator = 1:mm.n_size
        curr_trys = mm.n_size+1-tri_iterator;
        for succ = 1:curr_trys
    
            new_v = val_func_iter('no_learning',squeeze(V(succ,curr_trys+1,1,succ,prod_ind,:)),squeeze(V(succ+1,curr_trys+1,1,succ+1,prod_ind,:)),squeeze(V(succ,curr_trys+1,1,succ,prod_ind,:)),mm.theta1(k),policy.pi_f(prod_ind,:),succ,mm,policy,'foreign');
            V(succ,curr_trys,1,succ,prod_ind,:) = new_v;
            [~,l_opt(succ,curr_trys,1,succ,prod_ind,:)] = iterationError_h(squeeze(V(succ,curr_trys,1,succ,prod_ind,:)),mm.theta1(k),policy.pi_f(prod_ind,:),succ,squeeze(V(succ,curr_trys+1,1,succ,prod_ind,:)),mm,policy,'foreign');
    
        end
    end
end


end


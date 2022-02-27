function [V,l_opt] = solvePolicyForeign(policy,mm)

V = zeros(size(mm.Q_px_f,1),mm.n_size+1,mm.n_size+1,1,mm.net_size+1);
l_opt = zeros(size(V));

for succ = 1:mm.n_size+1

    l_opt(:,mm.n_size+1,succ,1,mm.net_size+1) = mm.l_opt_func_f(policy.postSuccessProb_f(mm.n_size+1,succ),mm.net_size+1,policy.pi_f,0,0,0);
    Q0_diag = -mm.Q_px_f_d;
    Q0_diag(1:size(mm.Q_px_f_d,1)+1:end) = mm.r+abs(diag(mm.Q_px_f));
    perm_flow = policy.postSuccessProb_f(mm.n_size+1,succ)*l_opt(:,mm.n_size+1,succ,1,mm.net_size+1).*policy.pi_f-mm.cost_f(l_opt(:,mm.n_size+1,succ,1,mm.net_size+1),mm.net_size+1);
    % next eq. follows from V = (mm.r + abs(diag(Q0)))^-1 * (perm_flow + Q0_d * V)
    V(:,mm.n_size+1,succ,1,mm.net_size+1) = max(Q0_diag^-1*perm_flow,0);

    % learning maxed out
    for m = 1:mm.net_size+1-succ

        new_v = val_func_iter(V(:,mm.n_size+1,succ,1,mm.net_size+2-m),mm.b,policy.postSuccessProb_f(mm.n_size+1,succ),policy.pi_f,mm.net_size+1-m,mm.Q_px_f,mm.r,abs(diag(mm.Q_px_f)),mm.Q_px_f_d,V(:,mm.n_size+1,succ,1,mm.net_size+2-m),1,mm.gam,1,mm.cs_f,mm.v_tolerance,mm.v_rel_tolerance,1,[],[],[],mm.cost_f,mm.l_opt_func_f,policy.pi_f);
        V(:,mm.n_size+1,succ,1,mm.net_size+1-m) = new_v;

        [~,~,l_opt(:,mm.n_size+1,succ,1,mm.net_size+1-m)] = iterateValueFunction_h(V(:,mm.n_size+1,succ,1,mm.net_size+1-m),...
            policy.postSuccessProb_f(mm.n_size+1,succ),policy.pi_f,mm.net_size+1-m,mm.r,abs(diag(mm.Q_px_f)),mm.Q_px_f_d,V(:,mm.n_size+1,succ,1,mm.net_size+2-m),mm.cost_f,mm.l_opt_func_f);
    end
end

for tri = 1:mm.n_size
    for succ = 1:mm.n_size+1-tri

        new_v = val_func_iter(V(:,mm.n_size+2-tri,succ,1,succ),mm.b,policy.postSuccessProb_f,policy.pi_f,mm.n_size+1-tri,mm.Q_px_f,mm.r,abs(diag(mm.Q_px_f)),mm.Q_px_f_d,V(:,mm.n_size-tri+2,succ+1,1,succ+1),V(:,mm.n_size-tri+2,succ,1,succ),mm.gam,1,mm.cs_f,mm.v_tolerance,mm.v_rel_tolerance,0,mm.n_size,tri,succ,mm.cost_f,mm.l_opt_func_f,policy.pi_f);
        V(:,mm.n_size+1-tri,succ,1,succ) = new_v;
        [~,~,l_opt(:,mm.n_size+1-tri,succ,1,succ)] = iterateValueFunction_f(V(:,mm.n_size+1-tri,succ,1,succ),policy.postSuccessProb_f,policy.pi_f,mm.n_size+1-tri,succ,mm.r,abs(diag(mm.Q_px_f)),mm.Q_px_f_d,V(:,mm.n_size-tri+2,succ,1,succ),V(:,mm.n_size-tri+2,succ+1,1,succ+1),mm.cost_f,mm.l_opt_func_f);
    end
end

end


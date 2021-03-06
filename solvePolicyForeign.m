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

        new_v = val_func_iter('no_learning',V(:,mm.n_size+1,succ,1,mm.net_size+2-m),V(:,mm.n_size+1,succ,1,mm.net_size+2-m),[],policy.postSuccessProb_f(mm.n_size+1,succ),mm.net_size+1-m,mm,policy,'foreign');
        V(:,mm.n_size+1,succ,1,mm.net_size+1-m) = new_v;
        [~,l_opt(:,mm.n_size+1,succ,1,mm.net_size+1-m)] = iterationError_h(V(:,mm.n_size+1,succ,1,mm.net_size+1-m),...
            policy.postSuccessProb_f(mm.n_size+1,succ),mm.net_size+1-m,V(:,mm.n_size+1,succ,1,mm.net_size+2-m),mm,policy,'foreign');

    end
end

for tri = 1:mm.n_size
    for succ = 1:mm.n_size+1-tri

        new_v = val_func_iter('learning',V(:,mm.n_size+2-tri,succ,1,succ),V(:,mm.n_size-tri+2,succ+1,1,succ+1),V(:,mm.n_size-tri+2,succ,1,succ),policy.postSuccessProb_f(mm.n_size+1-tri,succ),succ,mm,policy,'foreign');
        V(:,mm.n_size+1-tri,succ,1,succ) = new_v;
        [~,l_opt(:,mm.n_size+1-tri,succ,1,succ)] = iterationError_f(V(:,mm.n_size+1-tri,succ,1,succ),policy.postSuccessProb_f(mm.n_size+1-tri,succ),succ,V(:,mm.n_size-tri+2,succ,1,succ),V(:,mm.n_size-tri+2,succ+1,1,succ+1),mm,policy);

    end
end

end


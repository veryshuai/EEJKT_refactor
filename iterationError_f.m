function [err,search] = iterationError_f(V_curr,a,pi,net,V_fail,V_succ,mm)

    search = mm.l_opt_func_f(a,net,pi,V_succ,V_fail,V_curr);  
    denom = mm.r+search+abs(diag(mm.Q_f));
    meeting_payoff = search.*(a*(pi+V_succ)+(1-a)*V_fail);
    mac_state_change = mm.Q_f_d*V_curr;
    err = denom.^-1.*(-mm.cost_f(search,net)+mac_state_change+meeting_payoff)-V_curr;
    
end 

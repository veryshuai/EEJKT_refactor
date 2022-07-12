function [err,search] = iterationError_h(V_curr,a,pi,net,V_succ,mm,policy,country)  
                                         
    [~,l_opt_func,~,Q_m,Q_m_d,cost] = retrieveCountrySpecificParameters(country,mm,policy);

    search = l_opt_func(a,net,pi,V_succ,V_curr,V_curr);  
    denom = mm.r+search+abs(diag(Q_m));
    meeting_payoff = search.*(a*(pi+V_succ)+(1-a)*V_curr);
    mac_state_change = Q_m_d*V_curr;
    err = denom.^-1.*(-cost(search,net)+mac_state_change+meeting_payoff)-V_curr;

end

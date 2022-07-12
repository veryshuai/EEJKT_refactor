function [max_pi,l_opt_func,pi,Q_m,Q_m_d,cost] = retrieveCountrySpecificParameters(country,mm,policy)

    if country == "foreign"
        max_pi = max(policy.pi_f);
        l_opt_func = mm.l_opt_func_f;
        pi = policy.pi_f;
        Q_m = mm.Q_f;
        Q_m_d = mm.Q_f_d;
        cost = mm.cost_f;
    elseif country == "home"
        max_pi = max(policy.pi_h);
        l_opt_func = mm.l_opt_func_h;
        pi = policy.pi_h;
        Q_m = mm.Q_h;
        Q_m_d = mm.Q_h_d;
        cost = mm.cost_h;
    end

end
function [max_pi,l_opt_func,pi,Q_px,Q_px_d,cost] = retrieveCountrySpecificParameters(country,mm,policy)

    if country == "foreign"
        max_pi = max(policy.pi_f);
        l_opt_func = mm.l_opt_func_f;
        pi = policy.pi_f;
        Q_px = mm.Q_px_f;
        Q_px_d = mm.Q_px_f_d;
        cost = mm.cost_f;
    elseif country == "home"
        max_pi = max(policy.pi_h);
        l_opt_func = mm.l_opt_func_h;
        pi = policy.pi_h;
        Q_px = mm.Q_px_h;
        Q_px_d = mm.Q_px_h_d;
        cost = mm.cost_h;
    end

end
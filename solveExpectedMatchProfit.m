function [expected_match_profit,c_val] = solveExpectedMatchProfit(pie_scale,px,Q_px,Q_px_d,F,mm)

event_hazard = mm.r + mm.delta + mm.L_b + abs(Q_px(1,1)) + abs(mm.Q_z(1,1));

shipment_payoff_prod_x_no_z = (1/mm.eta) * exp(pie_scale) * exp((mm.eta-1)*px(1,:) + px(2,:)); %depends on prod and macro state x
shipment_payoff = repmat(exp(mm.Z)',size(shipment_payoff_prod_x_no_z,2),1) .* repmat(shipment_payoff_prod_x_no_z',1,size(mm.Z,1));

match_prof_z = 1 / (mm.delta + mm.r) * shipment_payoff;
c_val = match_prof_z;

myEps = 1e12;
it = 0;
max_iter = 50000;
while myEps > mm.pi_tolerance && it<=max_iter

    it = it + 1;
    if it == max_iter
        display('WARNING: makepie.m did not converge!');
    end

    c_val_gross_new = (Q_px_d * c_val + c_val * mm.Q_z_d' + mm.L_b * match_prof_z) / event_hazard;
    c_val_new = max(-F + c_val_gross_new,0);
    match_prof_z_new = shipment_payoff + c_val_new;

    myEps = max(max(abs(match_prof_z-match_prof_z_new)));
    match_prof_z = match_prof_z_new;
    c_val = c_val_new;

end

expected_match_profit = match_prof_z*mm.erg_pz;

end

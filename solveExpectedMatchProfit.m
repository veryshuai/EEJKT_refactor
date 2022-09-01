function [expected_match_profit,continuation_value] = solveExpectedMatchProfit(pie_scale,X,Q0,Q0_d,F,mm)
                                     
event_hazard = mm.r + mm.delta + mm.L_b + abs(Q0(1,1)) + abs(mm.Q_z(1,1)); 

shipment_payoff_except_z = exp(pie_scale + X + (mm.eta-1)*mm.Phi'); %dimension is macro shock by prod
shipment_payoff = repmat(repmat(exp(mm.Z),1,size(shipment_payoff_except_z,2)),1,1,size(shipment_payoff_except_z,1))...
    .* permute(repmat(shipment_payoff_except_z,1,1,size(mm.Z,1)),[3 2 1]); %dimension is demand shk by prod by macro

continuation_value = zeros(size(mm.Z,1),size(mm.Phi,1),size(X,1)); 
expected_match_profit = zeros(size(shipment_payoff_except_z))';

for j=1:size(mm.Phi)
    
    myEps = 1e12;
    it = 0;
    max_iter = 100000;
    c_val = zeros(size(mm.Z,1),size(X,1));
    pi_z = 1 / (mm.delta + mm.r) * squeeze(shipment_payoff(:,j,:)); 

    while myEps > mm.pi_tolerance && it<=max_iter

    it = it + 1;
    if it == max_iter
        display('WARNING: makepie.m did not converge!');
    end

    c_val_gross_new = (c_val * Q0_d' + mm.Q_z_d * c_val + mm.L_b * pi_z) / event_hazard ;          
    c_val_new = max(-F + c_val_gross_new,0);
    pi_z_new = squeeze(shipment_payoff(:,j,:)) + c_val_new;
    
    myEps = max(max(abs(pi_z-pi_z_new)));
    pi_z = pi_z_new; 
    c_val = c_val_new; 
    end 

    continuation_value(:,j,:) = max(c_val,0);
    expected_match_profit(j,:) = max(mm.erg_pz' * pi_z,0);
end

end 

function policy = solvePolicyExchangeRateShock(mm,policy)

    [policy.pi_f_exch_rate_shk,policy.c_val_f_exch_rate_shk] = solveExpectedMatchProfit(mm.exchange_rate_shock_multiplier + mm.scale_f,mm.L_bF,mm.X_f,mm.Q_f,mm.Q_f_d,mm.F_f,mm);

    [policy.value_f_exch_rate_shk,policy.lambda_f_exch_rate_shk] = solvePolicyForeign_exch_rate_shk(policy,mm); 

    %shouldHaveMonotonicPolicyFunctions(policy,mm);

end

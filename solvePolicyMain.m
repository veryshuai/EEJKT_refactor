function policy = solvePolicyMain(mm)

    policy = struct();

    [policy.pi_f,policy.c_val_f] = solveExpectedMatchProfit(mm.scale_f,mm.px_f,mm.Q_px_f,mm.Q_px_f_d,mm.F_f,mm); 
    [policy.pi_h,policy.c_val_h] = solveExpectedMatchProfit(mm.scale_h,mm.px_h,mm.Q_px_h,mm.Q_px_h_d,mm.F_h,mm); 

    policy.postSuccessProb_f = makeForeignSuccessPosteriors(mm); %[trials,successes]

    [policy.value_h,policy.lambda_h] = solvePolicyHome(policy,mm);  
    [policy.value_f,policy.lambda_f] = solvePolicyForeign(policy,mm); 

    policy = reshapePolicy(policy, mm);

end

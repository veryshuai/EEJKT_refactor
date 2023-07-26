function policy = solvePolicyMain(mm)

    policy = struct();

    [policy.pi_f,policy.c_val_f] = solveExpectedMatchProfit(mm.scale_f,mm.L_bF,mm.X_f,mm.Q_f,mm.Q_f_d,mm.F_f,mm);
    [policy.pi_h,policy.c_val_h] = solveExpectedMatchProfit(mm.scale_h,mm.L_bH,mm.X_h,mm.Q_h,mm.Q_h_d,mm.F_h,mm); 

    policy.postSuccessProb_f = makeForeignSuccessPosteriors(mm); %[trials,successes]

    tic
    [policy.value_h,policy.lambda_h] = solvePolicyHome(policy,mm);  
    [policy.value_f,policy.lambda_f] = solvePolicyForeign_nolearning(policy,mm); 
    toc

    %shouldHaveMonotonicPolicyFunctions(policy,mm);

end

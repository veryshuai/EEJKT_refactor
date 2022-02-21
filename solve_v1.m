function policy = solve_v1(mm)

    policy = struct();

    [policy.pi_f,policy.c_val_f] = solveExpectedMatchProfit(mm.scale_f,mm.X_f,mm.Q_f,mm.Q_f_d,mm.F_f,mm); %pi[prod,macro]
    [policy.pi_h,policy.c_val_h] = solveExpectedMatchProfit(mm.scale_h,mm.X_h,mm.Q_h,mm.Q_h_d,mm.F_h,mm); %c_val[demand shk, prod,macro]
    
    shouldMatchMoments(pi_f,c_val_f,"test","results/shouldMatchProfitsForeignData");
    shouldMatchMoments(pi_h,c_val_h,"test","results/shouldMatchProfitsHomeData");

    policy.postSuccessProb_f = makeForeignSuccessPosteriors(mm); %[trials,successes]

    [policy.value_h,policy.lambda_h] = solvePolicyHome(policy,mm);  
    [policy.value_f,policy.lambda_f] = solvePolicyForeign(policy,mm); 

    shouldMatchMoments(val_f,lf,"test","results/shouldMatchPolicyForeignData");
    shouldMatchMoments(val_h,lh,"test","results/shouldMatchPolicyHomeData");

end

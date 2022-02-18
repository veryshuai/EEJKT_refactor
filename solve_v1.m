function [lambda_f,lambda_h,pi_h,pi_f,c_val_h,c_val_f,my_flag,value_h,value_f] = solve_v1(mm)

    policy = struct();
        
    [policy.pi_f,policy.c_val_f] = solveExpectedMatchProfit(mm.scale_f,mm.X_f,mm.Q_f,mm.Q_f_d,mm.F_f,mm); %pi[prod,macro]
    [policy.pi_h,policy.c_val_h] = solveExpectedMatchProfit(mm.scale_h,mm.X_h,mm.Q_h,mm.Q_h_d,mm.F_h,mm); %c_val[prod,macro,demand shk]
    
    policy.postSuccessProb_f = makeForeignSuccessPosteriors(mm); %[trials,successes]
    
    [policy.value_h,policy.lambda_h] = solvePolicyHome(policy,mm);  
%   [~,lf,flag_f] = val_loop_f(Q0_f,Q0_f_d,a_f,pi_f,mm); 

end

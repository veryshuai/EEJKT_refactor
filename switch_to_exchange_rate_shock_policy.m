function policy = switch_to_exchange_rate_shock_policy(policy)

%keep pre-shock values
policy.pi_f_preshock = policy.pi_f;
policy.c_val_f_preshock = policy.c_val_f;
policy.value_f_preshock = policy.value_f;
policy.lambda_f_preshock = policy.lambda_f;
policy.pmat_cum_f_preshock = policy.pmat_cum_f;

%Implement shock
policy.pi_f = policy.pi_f_exch_rate_shk;
policy.c_val_f = policy.c_val_f_exch_rate_shk;
policy.value_f = policy.value_f_exch_rate_shk;
policy.lambda_f = policy.lambda_f_exch_rate_shk;
policy.pmat_cum_f = policy.pmat_cum_f_exch_rate_shk;

end
    

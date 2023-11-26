function check_FOC(mm,policy)

    max_abs_error = 0;
    for trials = 1:20
        for succs = 1:trials
            pol = policy.lambda_f(succs,trials,1,succs,17,8);
            val_base = policy.value_f(succs,trials,1,succs,17,8);
            val_succ = policy.value_f(succs+1,trials+1,1,succs+1,17,8);
            val_fail = policy.value_f(succs,trials+1,1,succs,17,8);
            %slope_cost = (mm.cost_f(pol + 1e-12,1) - mm.cost_f(pol,1))/1e-12;
            slope_cost = mm.cs_f * pol^(mm.kappa1 - 1);
            my_theta = (mm.af + succs - 1) / (mm.af + mm.bf + trials - 1);
            FOC_RHS = my_theta * (policy.pi_f(17,8) + val_succ) + (1 - my_theta) * val_fail - val_base;
            FOC_error = slope_cost - FOC_RHS;
            if abs(FOC_error) > max_abs_error
                max_abs_error = FOC_error;
                offender_succs = succs;
                offender_trials = trials;
                %rel error (from value function solver)
                %differ = norm(new_v - old_v)/norm(old_v);
                rel_error_at_max = max_abs_error / val_base;
            end
        end
    end
    display(max_abs_error)
    display(rel_error_at_max)
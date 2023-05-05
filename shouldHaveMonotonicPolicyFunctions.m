function shouldHaveMonotonicPolicyFunctions(policy,mm)

%this script tests for monotonicity in EEJKT policy functions

% reminder of what the arguments are
% lambda_f (succ, trial, common succ rate (defunct), network size, prod of firm, F macro shock) 
% lambda_h (common succ rate (defunct), known succ rate, network size, prod of firm, H macro shock)
% c_val* (demand shock, prod of firm, macro shock)

%conditional on trials, search should be increasing in success

%foreign
for prod_ind = 1:2*mm.phi_size + 1
    for mac_shk_ind = 1:2*mm.x_size + 1
        for trial_ind = 2:mm.n_size + 1
            for succ_ind = 1:trial_ind
                
                %monotonic in success
                if succ_ind > 1 && mm.gam >= 0
                    pol_diff = policy.lambda_f(succ_ind,trial_ind,1,succ_ind,prod_ind,mac_shk_ind)...
                        - policy.lambda_f(succ_ind-1,trial_ind,1,succ_ind,prod_ind,mac_shk_ind);
                    if pol_diff < 0
                        disp('Warning: Foreign policy function decreasing in successes conditional on trials')
                        disp(strcat('successes:' ,num2str(succ_ind)));
                        disp(strcat('trials:' ,num2str(trial_ind)));
                        disp(strcat('mac shock:' ,num2str(mac_shk_ind)));
                        disp(strcat('prod_ind:' ,num2str(prod_ind)));
                    end
                end
                
                %monotonic in macro shocks
                if mac_shk_ind > 1
                    pol_diff = policy.lambda_f(succ_ind,trial_ind,1,succ_ind,prod_ind,mac_shk_ind)...
                        - policy.lambda_f(succ_ind,trial_ind,1,succ_ind,prod_ind,mac_shk_ind-1);
                    if pol_diff < 0
                        disp('Warning: Foreign policy function decreasing in macro shocks')
                        disp(strcat('successes:' ,num2str(succ_ind)));
                        disp(strcat('trials:' ,num2str(trial_ind)));
                        disp(strcat('mac shock:' ,num2str(mac_shk_ind)));
                        disp(strcat('prod_ind:' ,num2str(prod_ind)));
                    end
                end

                %monotonic in productivity
                if prod_ind > 1
                    pol_diff = policy.lambda_f(succ_ind,trial_ind,1,succ_ind,prod_ind,mac_shk_ind)...
                        - policy.lambda_f(succ_ind,trial_ind,1,succ_ind,prod_ind-1,mac_shk_ind);
                    if pol_diff < 0
                        disp('Warning: Foreign policy function decreasing in productivity')
                        disp(strcat('successes:' ,num2str(succ_ind)));
                        disp(strcat('trials:' ,num2str(trial_ind)));
                        disp(strcat('mac shock:' ,num2str(mac_shk_ind)));
                        disp(strcat('prod_ind:' ,num2str(prod_ind)));
                    end
                end

                %monotonic in trials conditional on success
                if trial_ind > 1 && succ_ind < trial_ind
                    pol_diff = policy.lambda_f(succ_ind,trial_ind,1,succ_ind,prod_ind,mac_shk_ind)...
                        - policy.lambda_f(succ_ind,trial_ind-1,1,succ_ind,prod_ind,mac_shk_ind);
                    if pol_diff > 0
                        disp('Warning: Foreign policy function increasing in trials conditional on success')
                        disp(strcat('successes:',num2str(succ_ind)));
                        disp(strcat('trials:' ,num2str(trial_ind)));
                        disp(strcat('mac shock:' ,num2str(mac_shk_ind)));
                        disp(strcat('prod_ind:' ,num2str(prod_ind)));
                        disp(strcat('pct increase' , num2str(pol_diff/policy.lambda_f(succ_ind,trial_ind-1,1,succ_ind,prod_ind,mac_shk_ind))));
                    end
                end              
            end
        end
    end
end

%home
for prod_ind = 1:2*mm.phi_size + 1
    for mac_shk_ind = 1:2*mm.x_size + 1
        for succ_rt_ind = 1:mm.dim1
            for succ_ind = 1:succ_rt_ind
                
                %monotonic in success (gamma > 0)
                if succ_ind > 1 && mm.gam >=0
                    pol_diff = policy.lambda_h(1,succ_rt_ind,succ_ind,prod_ind,mac_shk_ind)...
                        - policy.lambda_h(1,succ_rt_ind,succ_ind-1,prod_ind,mac_shk_ind);
                    if pol_diff < 0
                        disp('Warning: Home policy function decreasing in successes (gamma > 0)')
                        disp(strcat('successes:' ,num2str(succ_ind)));
                        disp(strcat('success rate index:' ,num2str(succ_rt_ind)));
                        disp(strcat('mac shock:' ,num2str(mac_shk_ind)));
                        disp(strcat('prod_ind:' ,num2str(prod_ind)));
                    end
                end
                
                %monotonic in success (gamma < 0)
                if succ_ind > 1 && mm.gam < 0
                    pol_diff = policy.lambda_h(1,succ_rt_ind,succ_ind,prod_ind,mac_shk_ind)...
                        - policy.lambda_h(1,succ_rt_ind,succ_ind-1,prod_ind,mac_shk_ind);
                    if pol_diff > 0
                        disp('Warning: Home policy function increasing in successes (gamma < 0)')
                        disp(strcat('successes:' ,num2str(succ_ind)));
                        disp(strcat('success rate index:' ,num2str(succ_rt_ind)));
                        disp(strcat('mac shock:' ,num2str(mac_shk_ind)));
                        disp(strcat('prod_ind:' ,num2str(prod_ind)));
                    end
                end
                
                %monotonic in productivity
                if prod_ind > 1
                    pol_diff = policy.lambda_h(1,succ_rt_ind,succ_ind,prod_ind,mac_shk_ind)...
                        - policy.lambda_h(1,succ_rt_ind,succ_ind,prod_ind-1,mac_shk_ind);
                    if pol_diff < 0
                        disp('Warning: Home policy function decreasing in productivity')
                        disp(strcat('successes:' ,num2str(succ_ind)));
                        disp(strcat('success rate index:' ,num2str(succ_rt_ind)));
                        disp(strcat('mac shock:' ,num2str(mac_shk_ind)));
                        disp(strcat('prod_ind:' ,num2str(prod_ind)));
                    end
                end
                
                %monotonic in macro shock
                if mac_shk_ind > 1
                    pol_diff = policy.lambda_h(1,succ_rt_ind,succ_ind,prod_ind,mac_shk_ind)...
                        - policy.lambda_h(1,succ_rt_ind,succ_ind,prod_ind,mac_shk_ind-1);
                    if pol_diff < 0
                        disp('Warning: Home policy function decreasing in macro shock')
                        disp(strcat('successes:' ,num2str(succ_ind)));
                        disp(strcat('success rate index:' ,num2str(succ_rt_ind)));
                        disp(strcat('mac shock:' ,num2str(mac_shk_ind)));
                        disp(strcat('prod_ind:' ,num2str(prod_ind)));
                    end
                end

                %monotonic in success rate index
                if succ_rt_ind > 1
                    pol_diff = policy.lambda_h(1,succ_rt_ind,succ_ind,prod_ind,mac_shk_ind)...
                        - policy.lambda_h(1,succ_rt_ind-1,succ_ind,prod_ind,mac_shk_ind);
                    if pol_diff < 0
                        disp('Warning: Home policy function decreasing in success rate')
                        disp(strcat('successes:' ,num2str(succ_ind)));
                        disp(strcat('success rate index:' ,num2str(succ_rt_ind)));
                        disp(strcat('mac shock:' ,num2str(mac_shk_ind)));
                        disp(strcat('prod_ind:' ,num2str(prod_ind)));
                    end
                end

            end
        end
    end
end

end

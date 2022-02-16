% This script gives the user options over what policy function related task to perform

% Query user for desired task
acceptable = 'FALSE'; %Is user input interpretable?
while acceptable == 'FALSE'
    task = input(['There are two things you can currently do.  You can either plot the policy function, or you can estimate the fraction of exporters which would have searched had they known their product appeal, but due to bad luck give up searching.' char(10) 'est (estimate)' char(10) 'plt (plot)' char(10) 'Which task shall I perform?: '], 's');

    display([char(10) 'Your input was ' num2str(task) '.' char(10)]);

    %Check for validity
    if (task == 'est') | (task == 'plt') 
        acceptable = 'TRUE ';
    else
        display(['Sorry, I do not understand.  Try again.' char(10)]);
    end
end

if task == 'plt'

    % Get policy function
    [lambda_f_orig,lambda_f_new,pi_tilda_f_orig,pi_tilda_f_new,c_val_f_orig,c_val_f_new,punishment] = solve(mm);

    % Plot the policy function
    case_str = 'non';
    increase = 1; %no counterfactual
    policy_plot('baseline',lambda_f_orig,0,lambda_f_orig);

elseif task == 'est'

    % Estimate the missed opportunities
    case_str = 'pol';
    SetParams_noprod;

    % Calculate the full information foreign lambdas
    [lambda_f,lambda_f_full_info,pi_tilda_full_info,pi_tilda_f,c_val_f_full_info,c_val_f,punishment] = solve(mm);

    % Calculate the missed opportunities    
    calc_missed_ops;

    % Get no learning sales simulation
    lambda_f_new   =  lambda_f;
    lambda_h_new   =  lambda_f_full_info;
    pi_tilda_h_new =  pi_tilda_full_info;
    pi_tilda_f_new =  pi_tilda_f;
    c_val_h_new    =  c_val_f_full_info;
    c_val_f_new    =  c_val_f;
    punishment_new =  punishment;
    lambda_h_orig   =  lambda_f_full_info;
    pi_tilda_h_orig =  pi_tilda_full_info;
    c_val_h_orig    =  c_val_f_full_info;
    punishment_orig =  punishment;
    
end

%Inform user that an error after this script is run is normal
display('NOTE: AN ERROR FOLLOWING THIS MESSAGE IS NORMAL!')



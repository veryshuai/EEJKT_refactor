% Executed by inten_sim_v1, this script calls solve_v1 to generate policy
% functions for both markets.

% First half policy 

case_str = 'non';

% Get policy and value functions
if cf_num ~= 9 % Don't waste time if we are doing no learning version of the model
    [lambda_f_orig,lambda_h_orig,pi_tilda_h_orig,pi_tilda_f_orig,c_val_h_orig,c_val_f_orig,flag_orig,value_h,value_f] = solve_v1(mm);    
end

% Second half policy

switch cf_num
    case 0
        case_str = 'est';
        increase = 1; % no counterfactual

        % Get policy and value functions
        lambda_f_new   =  lambda_f_orig;
        lambda_h_new   =  lambda_h_orig;
        pi_tilda_h_new =  pi_tilda_h_orig;
        pi_tilda_f_new =  pi_tilda_f_orig;
        c_val_h_new    =  c_val_h_orig;
        c_val_f_new    =  c_val_f_orig;
        flag_new =  flag_orig;

    case 1
        case_str = 'sim';
        increase = 1; % no counterfactual

        % Get parameters (firm number)
        SetParams_noprod;

        % Get policy and value functions
        lambda_f_new   =  lambda_f_orig;
        lambda_h_new   =  lambda_h_orig;
        pi_tilda_h_new =  pi_tilda_h_orig;
        pi_tilda_f_new =  pi_tilda_f_orig;
        c_val_h_new    =  c_val_h_orig;
        c_val_f_new    =  c_val_f_orig;
        flag_new =  flag_orig;

    case 2
      %cost_dec_trans(1);

    case 3

        case_str = 'mac';

        % Set counterfactual
        increase = 1.2; % 20% jump in macro shock
        scale_f    =  scale_f + log(increase);

        % Get parameters
 %      SetParams_noprod; % commented this out. Don't know why it didn't throw an error.
        SetParams;       
 
        %reset random seed
        if seed == 1
            rng(80085);
        end
    
        % Get policy and value functions
        [lambda_f_new,lambda_h_new,pi_tilda_h_new,pi_tilda_f_new,c_val_h_new,c_val_f_new,flag_new] = solve_v1(mm);

    case 4
      %red_var_trans(1);

    case 5
      %search_dec_trans(1);

    case 6 % calculate value of network

        case_str = 'val';
        increase = 1; % no counterfactual

        % Get parameters
        SetParams_noprod;

        % Get policy and value functions
        lambda_f_new   =  lambda_f_orig;
        lambda_h_new   =  lambda_h_orig;
        pi_tilda_h_new =  pi_tilda_h_orig;
        pi_tilda_f_new =  pi_tilda_f_orig;
        c_val_h_new    =  c_val_h_orig;
        c_val_f_new    =  c_val_f_orig;
        flag_new =  flag_orig;

    case 7 % plot policies 

        increase = 1; % no counterfactual
        run_policy_stuff;

    case 8 % bootstrap for std err 

        case_str = 'boo';

        increase = 1; % no counterfactual

        % Get parameters
        SetParams_noprod;

        % Get policy and value functions
        lambda_f_new   =  lambda_f_orig;
        lambda_h_new   =  lambda_h_orig;
        pi_tilda_h_new =  pi_tilda_h_orig;
        pi_tilda_f_new =  pi_tilda_f_orig;
        c_val_h_new    =  c_val_h_orig;
        c_val_f_new    =  c_val_f_orig;
        flag_new =  flag_orig;

    case 9

        case_str = 'nln';

        % Get policy and value functions
        [lambda_f_orig,lambda_h_orig,pi_tilda_h_orig,pi_tilda_f_orig,c_val_h_orig,c_val_f_orig,flag_orig] = solve_nln(mm);

        increase = 1; % no counterfactual 

        % Get parameters
        SetParams_noprod;
        
        %reset random seed
        if seed == 1
            rng(80085);
        end
    
        % Get policy and value functions
        lambda_f_new   =  lambda_f_orig;
        lambda_h_new   =  lambda_h_orig;
        pi_tilda_h_new =  pi_tilda_h_orig;
        pi_tilda_f_new =  pi_tilda_f_orig;
        c_val_h_new    =  c_val_h_orig;
        c_val_f_new    =  c_val_f_orig;
        flag_new =  flag_orig;

    otherwise

        case_str = 'non';
        increase = 1; % no counterfactual

        % Get policy and value functions
        lambda_f_new   =  lambda_f_orig;
        lambda_h_new   =  lambda_h_orig;
        pi_tilda_h_new =  pi_tilda_h_orig;
        pi_tilda_f_new =  pi_tilda_f_orig;
        c_val_h_new    =  c_val_h_orig;
        c_val_f_new    =  c_val_f_orig;
        flag_new =  flag_orig;

end

% Make flag variable
try
    my_flag = max(flag_orig,flag_new);
catch e
    my_flag = flag_orig;
end

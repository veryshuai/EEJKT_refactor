function percent_missing = calc_single_missing_prob(prd,succ_prob,history_list,lambda_f)
% This function returns the fraction of firms with productivity prd and success prob succ_prob that stop searching due to negative signals 

%search at all?
if lambda_f{1,1,1,1}(prd,8) == 0
    percent_missing = 1;
    return;
end

%loop through events
run_sum_quitters = 0; %running sum of prob contribution of quitters
for hist_ind = 1:size(history_list)

    % Get history
    my_hist = history_list{hist_ind}';

    % Get partial Success sums
    succ = zeros(size(my_hist,1),1); %preallocate
    succ(1) = str2double(my_hist(1)); %get first trial outcome
    for trial = 2:size(my_hist)
        succ(trial) = succ(trial-1) + str2double(my_hist(trial));
    end

    % Get probability
    hist_prob = succ_prob^succ(end,1)  *  (1 - succ_prob)^(size(my_hist,1) - succ(end,1));

    % Check for firms giving up search
    gave_up = 0;
    for trials = 1:size(succ)
        if lambda_f{succ(trials)+1,trials+1,1,succ(trials)+1}(prd,8) == 0 
            gave_up = 1;
            break;
        end
    end

    % Add to quitters
    run_sum_quitters = run_sum_quitters + hist_prob * gave_up;
end

percent_missing = run_sum_quitters;


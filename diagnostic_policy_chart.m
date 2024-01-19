function diagnostic_policy_chart(policy,lab)
%policy chart

    initial_search = zeros(10,1);
    initial_search(1) = 0.0474; %learning model
    initial_search(2) = 0.0099; %one failure
    initial_search(3) = 0.1505; %one success
    for k = 4:10
        initial_search(k) = policy.lambda_f(1,1,k-3,1,17,15);
    end
    initial_search_labels = ["base", "baseFail", "baseSucc", "nol1", "nol2", "nol3", "nol4", "nol5","nol6", "nol7", "nol8", "nol9", "nol10"];
    bar(initial_search);
    xticklabels(initial_search_labels);
    saveas(gcf, append('results/diagnostic_policy_plot_', lab, '.png'));

end
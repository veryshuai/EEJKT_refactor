function [elasticities,standard_errors] = exch_rate_analysis_calc_elas(deflator,match_recs,sim_source,mm);
%This function calculates elasticities at the firm level, and bootstraps
%standard errors
%For elasticity matrix, type of elasticity is row (sales,matches,firms)
%column is length (1 year, 3 year, 10 year)

rng(12345);

% Calculate elasticities
option_firms_per_sim = 1;
[elasticities,firms_per_sim] = exch_rate_analysis_calc_elas_inner(deflator,match_recs,sim_source,mm,option_firms_per_sim);

%% Bootstrap standard errors
%Strategy: Draw simulation rounds at random, keep adding firm numbers until
%we get to 3300 as in the final year of our data (in year 25 of the
%simulation).  Do this 1000 times.

bootstrapped_elasticities = zeros(1000,3,3); %bootstrap round,elasticity length,elasticity type (sales,matches,firms)
option_firms_per_sim = 0;
for k=1:1000
    display("bootstrap_iteration " + k)
    match_recs_boot = [];
    boot_firms = 0;
    while boot_firms < 3300
        sim_no_boot = randi([1,2000]);
        match_recs_boot = [match_recs_boot;match_recs(sim_source == sim_no_boot,:)];
        boot_firms = boot_firms + firms_per_sim(sim_no_boot,1); %update with number of firms just before shock
    end
    [boot_elas,~] = exch_rate_analysis_calc_elas_inner(deflator,match_recs_boot,sim_source,mm,option_firms_per_sim);
    bootstrapped_elasticities(k,:,:) = boot_elas;
end

standard_errors = zeros(3,3);
for elas_indx = 1:3
    for time_indx = 1:3
        standard_errors(elas_indx,time_indx) = std(bootstrapped_elasticities(:,elas_indx,time_indx));
    end
end

end
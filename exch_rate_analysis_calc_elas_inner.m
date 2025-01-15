function varargout  = exch_rate_analysis_calc_elas_inner(deflator,match_recs,sim_source,mm,option_firms_per_sim)

%% Prepare data
firms_per_sim = zeros(2000,4);
elas_indx = 0;
for t = [25,26,28,35]
    elas_indx = elas_indx + 1;

    match_recs_one_year = match_recs(match_recs(:,11) == t,:);
    sim_source_one_year = sim_source(match_recs(:,11) == t,:);
    new_id_one_year = match_recs_one_year(:,12);
    [unique_firms,uniq_firm_rows,~] = unique(new_id_one_year);
    total_firms(t) = size(unique_firms,1);
    total_matches(t) = size(match_recs_one_year,1);
    
    if option_firms_per_sim == 1
        for k = 1:2000
            firms_per_sim(k,elas_indx) = sum(sim_source_one_year(uniq_firm_rows)==k);
        end
    end
    haircut = 1;
    if t>25
        haircut = 6/5;
    end
    total_sales(t) = deflator * sum(haircut * match_recs_one_year(:,4),1);
 
end

%% elasticities
%SALES
elasticities = [(log(total_sales(26)) - log(total_sales(25)))/(log(1.2) - log(1)),...
    (log(total_sales(28)) - log(total_sales(25)))/(log(1.2) - log(1)),...
    (log(total_sales(35)) - log(total_sales(25)))/(log(1.2) - log(1));...
... %MATCHES
(log(total_matches(26)) - log(total_matches(25)))/(log(1.2) - log(1)),...
    (log(total_matches(28)) - log(total_matches(25)))/(log(1.2) - log(1)),...
    (log(total_matches(35)) - log(total_matches(25)))/(log(1.2) - log(1));...
... %FIRMS
(log(total_firms(26)) - log(total_firms(25)))/(log(1.2) - log(1)),...
    (log(total_firms(28)) - log(total_firms(25)))/(log(1.2) - log(1)),...
    (log(total_firms(35)) - log(total_firms(25)))/(log(1.2) - log(1))];

    varargout = {elasticities,firms_per_sim};

end
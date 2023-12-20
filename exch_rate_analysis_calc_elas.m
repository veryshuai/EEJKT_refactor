function [sales_elas,match_elas,firm_elas] = exch_rate_analysis_calc_elas(deflator,match_recs,mm) 

    match_recs(:,11) = match_recs(:,1)/mm.pd_per_yr; %years rather than months
    new_id = 0.5*(match_recs(:,2)+match_recs(:,3)).*(match_recs(:,2)+match_recs(:,3)+1)+match_recs(:,3);
    match_recs = [match_recs,new_id];
    
    match_recs = exch_shock_analysis_alt_age_calc(match_recs);
    
    for t=[25,26,28,35]
        
        match_recs_one_year = match_recs(match_recs(:,11) == t,:);
        %match_recs_one_year = succ_matches(succ_matches(:,11) == t,:);

        new_id_one_year = match_recs_one_year(:,12);
        [unique_firms,~,~] = unique(new_id_one_year);
        total_firms(t) = size(unique_firms,1);
        
        total_matches(t) = size(match_recs_one_year,1);
        haircut = 1;
        if t > 25
            haircut = 6/5;
        end
        total_sales(t) = deflator * sum(haircut * match_recs_one_year(:,4),1);
        
        
    end
    %Elasticity calculations
    sales_elas = zeros(3,1);
    match_elas = zeros(3,1);
    firm_elas = zeros(3,1);
    for k = [26,28,35]
        sales_elas(k) = (log(total_sales(k)) - log(total_sales(25)))/(log(1.2) - log(1));
        match_elas(k) = (log(total_matches(k)) - log(total_matches(25)))/(log(1.2) - log(1));
        firm_elas(k) = (log(total_firms(k)) - log(total_firms(25)))/(log(1.2) - log(1));
    end

end
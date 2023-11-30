function data_sorted = exch_shock_analysis_alt_age_calc(data)

 %Alternative age calculation
    data_sorted = sortrows(data,[12 1]);
    firm_age_in_yrs_alt = zeros(size(data_sorted,1),1);
    for ind = 2:size(data_sorted,1)
        if data_sorted(ind-1,12) == data_sorted(ind,12)
            if data_sorted(ind-1,1) < data_sorted(ind,1)
                firm_age_in_yrs_alt(ind) = firm_age_in_yrs_alt(ind-1) + 1;
            else
                firm_age_in_yrs_alt(ind) = firm_age_in_yrs_alt(ind-1);
            end
        end
    end

    data_sorted(:,9) = firm_age_in_yrs_alt * 12 + 1; %to correspond to the other obs

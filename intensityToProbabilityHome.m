function pmat_cum = intensityToProbabilityHome(mm,policy)

Q_size_h = 70;
Q_index_h = [(1:Q_size_h)',zeros(Q_size_h,1),(1:Q_size_h)'];

%% construct intensity matrix for each firm type

qmat_type = cell(size(policy.firm_type_prod_succ_macro,1),1);
pmat_type = cell(size(policy.firm_type_prod_succ_macro,1),1);
pmat_cum  = cell(size(policy.firm_type_prod_succ_macro,1),1);

%% Create a matrix for taking out rounding errors in matrix exponentiation
nonzero_tran = zeros(Q_size_h);
for i = 1:Q_size_h
    for j = 1:Q_size_h
        cond1 = (Q_index_h(j,2) - Q_index_h(i,2)<0)*(Q_index_h(j,1)>1);
        cond2 = (Q_index_h(j,3) - Q_index_h(i,3)<0)*(Q_index_h(j,1)>1);
        imposs_tran = max(cond1,cond2);
        nonzero_tran(i,j) = 1 - imposs_tran;
    end
end
%% Now ready to loop over types
for typ_indx = 1:size(policy.firm_type_prod_succ_macro,1)
    %set type
    ms             = policy.firm_type_prod_succ_macro(typ_indx,2);
    succ_prob      = mm.theta2(policy.firm_type_prod_succ_macro(typ_indx,3));
    prod_lvl       = policy.firm_type_prod_succ_macro(typ_indx,4);

    counter = 0;

        q_index_list = zeros(mm.n_size + 1,3);
        for ss=1:1:Q_size_h % number of successes, plus 1
            if ss <= mm.n_size + 1
                search_inten = policy.lambda_h(1,policy.firm_type_prod_succ_macro(typ_indx,3),ss,prod_lvl,ms);
            else
                search_inten = policy.lambda_h(1,policy.firm_type_prod_succ_macro(typ_indx,3),mm.n_size + 1,prod_lvl,ms);
            end

            % deal separately with first and last nn1. It is impossible to learn more from last trial
            if ss ~= Q_size_h && ss ~= 1
                counter = counter + 1;
                q_index_list(counter,:) = [ss,ss,-(succ_prob*search_inten + mm.firm_death_haz)]; %diag
                counter = counter + 1;
                q_index_list(counter,:) = [ss,ss+1,succ_prob*search_inten]; %trial and success
                counter = counter + 1;
                q_index_list(counter,:) = [ss,1,mm.firm_death_haz]; %firm death
            elseif ss == Q_size_h
                counter = counter + 1;
                q_index_list(counter,:) = [ss,ss,-(mm.firm_death_haz)];
                counter = counter + 1;
                q_index_list(counter,:) = [ss,1,mm.firm_death_haz];
            elseif ss == 1
                counter = counter + 1;
                q_index_list(counter,:) = [1,1,-succ_prob*search_inten]; %no firm death hazard, no success
                counter = counter + 1;
                q_index_list(counter,:) = [1,2,succ_prob*search_inten]; %trial and success
            end
        end
        qmat_type{typ_indx} = sparse(q_index_list(:,1),q_index_list(:,2),q_index_list(:,3));
    %% construct transition probabilities for discrete time intervals

    pmat_type{typ_indx} = nonzero_tran.*expm(qmat_type{typ_indx});
    pmat_type{typ_indx} = pmat_type{typ_indx}./sum(pmat_type{typ_indx},2);
    pmat_cum{typ_indx}  = cumsum(pmat_type{typ_indx}')';
end

end





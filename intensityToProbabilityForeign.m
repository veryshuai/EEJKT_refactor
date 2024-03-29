function policy = intensityToProbabilityForeign(mm,policy)


%% construct intensity matrix for each firm type

qmat_type = cell(size(policy.firm_type_macro_succ_prod,1),1);
pmat_type = cell(size(policy.firm_type_macro_succ_prod,1),1);
pmat_cum  = cell(size(policy.firm_type_macro_succ_prod,1),1);

for typ_indx = 1:size(policy.firm_type_macro_succ_prod,1)
    ms             = policy.firm_type_macro_succ_prod(typ_indx,2);
    succ_prob      = mm.theta2(policy.firm_type_macro_succ_prod(typ_indx,3));
    prod_lvl       = policy.firm_type_macro_succ_prod(typ_indx,4);

    counter = 0;

    q_index_list = zeros((mm.n_size + 1)*(mm.n_size + 2)/2,3);
    for i=1:mm.n_size + 1  % number of meetings, plus 1
        pos_i = sum(1:i-1); %advance to the correct number of meetings
        for ss=1:i % number of successes, plus 1
            r_ind = pos_i+ss;
            search_inten = policy.lambda_f(ss,i,1,ss,prod_lvl,ms);
            if i ~= mm.n_size + 1 && i ~= 1 % deal separately with first and last nn1.  It is impossible to learn more from last trial
                counter = counter + 1;
%                q_index_list(counter,:) = [r_ind,r_ind,-(search_inten + mm.firm_death_haz)]; %diag
                q_index_list(counter,:) = [r_ind,r_ind,-search_inten]; %diag
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,r_ind+i,(1-succ_prob)*search_inten]; %trial and failure
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,r_ind+i+1,succ_prob*search_inten]; %trial and success
%                counter = counter + 1;
%                q_index_list(counter,:) = [r_ind,1,mm.firm_death_haz]; %firm death
            elseif i == mm.n_size + 1
                counter = counter + 1;
%                q_index_list(counter,:) = [r_ind,r_ind,-(mm.firm_death_haz)]; %no learning (or searching for now)
                q_index_list(counter,:) = [r_ind,r_ind,-(mm.firm_death_haz)]; %no learning (or searching for now)
%                counter = counter + 1;
%                q_index_list(counter,:) = [r_ind,1,mm.firm_death_haz]; %firm death hazard
%                q_index_list(counter,:) = [r_ind,1,mm.firm_death_haz]; %firm death hazard
            elseif i == 1
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,r_ind,-(search_inten)]; %no firm death hazard, no clients to lose
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,r_ind+i,(1-succ_prob)*search_inten]; %trial and failure
                counter = counter + 1;
                q_index_list(counter,:) = [r_ind,r_ind+i+1,succ_prob*search_inten]; %trial and success
            end
        end
    end
    qmat_type{typ_indx} = sparse(q_index_list(:,1),q_index_list(:,2),q_index_list(:,3));

    %% construct transition probabilities for discrete time intervals


    policy.pmat_to_meets_succs = zeros((mm.n_size + 1) *(mm.n_size+2)/2,3); % [index,trials,successes]
    counter = 1;
    for i=1:1:mm.n_size+1    % number of meetings, plus 1
        for ss=1:1:i % number of successes, plus 1
            policy.pmat_to_meets_succs(counter,:) = [counter,i,ss];
            counter = counter + 1;
        end
    end

%    nonzero_tran = zeros(size(policy.pmat_to_meets_succs,1));%% Create a matrix for taking out rounding errors in matrix exponentiation
%    for i = 1:size(policy.pmat_to_meets_succs,1)
%        for j = 1:size(policy.pmat_to_meets_succs,1)
%            cond1 = (policy.pmat_to_meets_succs(j,2) - policy.pmat_to_meets_succs(i,2)<0)*(policy.pmat_to_meets_succs(j,1)>1);
%            cond2 = (policy.pmat_to_meets_succs(j,3) - policy.pmat_to_meets_succs(i,3)<0)*(policy.pmat_to_meets_succs(j,1)>1);
%            imposs_tran = max(cond1,cond2);
%            % when cum matches or cum successes fall, they must go all the way to zero
%            cond3 = ((policy.pmat_to_meets_succs(j,3) - policy.pmat_to_meets_succs(i,3)) > (policy.pmat_to_meets_succs(j,2) - policy.pmat_to_meets_succs(i,2)))*(policy.pmat_to_meets_succs(j,1)>1);
%            % # new successes cannot be greater than # new meetings for continuing firms
%            imposs_tran = max(imposs_tran,cond3);
%            nonzero_tran(i,j) = 1 - imposs_tran;
%        end
%    end

%    pmat_type{typ_indx} = nonzero_tran.*expm(qmat_type{typ_indx}); 
    pmat_type{typ_indx} = expm(qmat_type{typ_indx}); 
    pmat_type{typ_indx} = pmat_type{typ_indx}./sum(pmat_type{typ_indx},2);
    pmat_cum{typ_indx}  = cumsum(pmat_type{typ_indx}')';

end

policy.pmat_cum_f = pmat_cum;

end





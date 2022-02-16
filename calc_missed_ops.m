% This script takes two inputs -- the full information policy function lambda_f_full_info and the regular learning policy function lambda_f, then calculates which productivity product appeal firms will export under full information, and then the fraction of each of those firms that will not be exporting after 20 signals.

% Which firms export?  Average macro shock, prod appeal by prod

%read in actual success probabilities
th_f = mm.theta2;

%Preallocate
full_info_grid = zeros(size(lambda_f_full_info{1,1,1},1),size(lambda_f_full_info,2)); %rows prod, columns, prod_appeal

%Populate the full_info_grid from lambda_f
for k=1:size(lambda_f_full_info,2)
    temp = lambda_f_full_info{1,k,1};
    full_info_grid(:,k) = temp(:,7);
end

full_info_grid

%Find productivity rows that begin with zeros and end with positive numbers
%These are the firms that will only search if they believe they have high
%enough product appeal
zero_in_front = find(full_info_grid(:,1) == 0);
pos_in_back = find(full_info_grid(:,end) > 0);
affected_rows = intersect(zero_in_front,pos_in_back);
affected_rows

%Find the positive columns in the affected rows
affected_columns = cell(size(affected_rows));
for k=1:numel(affected_rows)
    ind = affected_rows(k);
    affected_cols{k} = find(full_info_grid(ind,:) > 0)';
end

%How many generations out should we look?
gens = 10;

%Preallocate results
result_grid = zeros(size(full_info_grid));

%Create list of possible meeting results
num_of_binaries = 2^gens;
history_list = cell(num_of_binaries,1);
for k = 1:num_of_binaries
    history_list{k} = sprintf(strcat('%0',num2str(gens),'d'),str2num(dec2bin(k-1)));
end

%Now we need to start looping through the affected firms
for prd_ind=1:numel(affected_rows)
    cols_for_this_row = affected_cols{prd_ind};
    for app_ind=1:numel(cols_for_this_row)
        succ_prob = th_f(cols_for_this_row(app_ind));
        missing_frac =  calc_single_missing_prob(affected_rows(prd_ind),succ_prob,history_list,lambda_f);
        result_grid(affected_rows(prd_ind), cols_for_this_row(app_ind)) = missing_frac;
    end
end

result_grid


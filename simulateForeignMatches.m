function [iter_out] = simulateForeignMatches(pt_ndx,macro_state_f,mm,policy)

[iter_in, iter_out] = simulateForeignInnerInitialize(mm, pt_ndx, macro_state_f);

iterF_check.seas_tran     = cell(mm.tot_yrs,1); 
iterF_check.match_mat     = cell(mm.tot_yrs,1);
iterF_check.Zcut_H        = cell(mm.tot_yrs,1);
iterF_check.mat_yr_sales  = cell(mm.tot_yrs,1);
iterF_check.firm_yr_sales = cell(mm.tot_yrs,1);

for t = 2:1:mm.periods

% try
%     assert(sum(abs(iter_in.cur_cli_cnt(:,t-1)-sum(iter_in.lag_cli_zst,2)))==0)
% catch
%     'cur_cli_cnt and lag_cli_zst inconsistent'
%         [pt_ndx,t]
%         [iter_in.cur_cli_cnt(:,t-1),sum(iter_in.lag_cli_zst,2)]
% end   
   
    iter_in.t = t;
    if mod(iter_in.t-1,mm.pd_per_yr) == 0
        iter_in.season = 1; % reset season when previous period completes a year
    end

    % update year, productivity type, and transition probability matrix
%   iter_in.year = floor((iter_in.t-1)/mm.pd_per_yr);
    iter_in.year = floor((iter_in.t-1)/mm.pd_per_yr)+1;
    iter_in.mic_type = find(policy.firm_type_macro_succ_prod(:,3) == mm.pt_type(pt_ndx,2) ...
        & policy.firm_type_macro_succ_prod(:,4) == mm.pt_type(pt_ndx,1),1,'first');
    p_mat_type = find(policy.firm_type_macro_succ_prod(:,2) == macro_state_f(iter_in.t)...
        & policy.firm_type_macro_succ_prod(:,3) == mm.pt_type(pt_ndx,2)...
        & policy.firm_type_macro_succ_prod(:,4) == mm.pt_type(pt_ndx,1),1,'first');
    iter_in.pmat_cum_t = policy.pmat_cum_f{p_mat_type}; % holds cumulative transition probs across (#success, #meeting) pairs.

    iter_in = simulateForeignMatchesInnerSim(iter_in,mm,policy);
    
    iter_in.keep_cli_lag = iter_in.keep_cli;
    
    if iter_in.season == mm.pd_per_yr
        
        [iter_in,iter_out] = simulateForeignMatchesInnerAnnualize(iter_in,iter_out,mm);
        
        [iter_in,iter_out] = simulateForeignMatchesInnerMoments(iter_in,iter_out,mm);
        

    if pt_ndx == mm.check_type    
    yr_ndx = iter_in.year;
    iterF_check.seas_tran{yr_ndx}     = iter_in.seas_tran;
    iterF_check.trans_count{yr_ndx}   = iter_in.trans_count;
    iterF_check.Zcut_H{yr_ndx}        = iter_in.seas_Zcut;
    iterF_check.mat_yr_sales{yr_ndx}  = iter_in.mat_yr_sales;
    iterF_check.firm_yr_sales{yr_ndx} = iter_in.firm_yr_sales;
    end
    
    end

    iter_in.season = iter_in.season + 1;

    iter_in.lag_cli_zst  = iter_in.cur_cli_zst;
    iter_in.new_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iter_in.die_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iter_in.trans_zst    = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iter_in.trans_count  = zeros(size(mm.Z,1)+1,size(mm.Z,1)+1,mm.sim_firm_num_by_prod_succ_type(pt_ndx));

if iter_in.t == mm.periods
%   for checking only: collect the count matrices in transF, by firm type
    find_xcli = find(sum(iter_in.cur_cli_cnt,2)>0);
    iter_out.transF{pt_ndx,1} = find_xcli;
    iter_out.transF{pt_ndx,2} = iter_in.cur_cli_cnt(find_xcli,:);
    iter_out.transF{pt_ndx,3} = iter_in.cum_succ(find_xcli,:);
    iter_out.transF{pt_ndx,4}  = iter_in.cumage(find_xcli,:);  
    iter_out.transF{pt_ndx,5}  = iter_in.new_firm(find_xcli,:);  
    iter_out.transF{pt_ndx,6} = iter_in.cum_meets(find_xcli,:);
    
% for checking only: rearrange count matrices in blocks, by firm # & type
    rooms =  iter_out.transF{pt_ndx,1};
    stackF = zeros(length(rooms),mm.periods+1); 
    for i=1:length(rooms)
    lb = (i-1)*4 + 1;
    ub = i*4;
    stackF(lb:ub,:) = ...  % for debugging only
    [rooms(i),iter_out.transF{pt_ndx,2}(i,:);...
    rooms(i),iter_out.transF{pt_ndx,3}(i,:);...
    rooms(i),iter_out.transF{pt_ndx,4}(i,:);... 
    rooms(i),iter_out.transF{pt_ndx,5}(i,:)];
    iter_out.stackF = stackF;
    end

    iter_out.iterF_check = iterF_check;   
end



end     
end

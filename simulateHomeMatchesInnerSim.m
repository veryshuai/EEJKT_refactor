function [iter_out] = simulateHomeMatchesInnerSim(iter_out, mm, iterH_in, pt_ndx, policy)

tic
iterH_in.pt_ndx = pt_ndx;
iterH_in.season = 1;

iterH_check.seas_tran       = cell(mm.tot_yrs,1); 
iterH_check.match_mat       = cell(mm.tot_yrs,1);
iterH_check.Zcut_eoyH       = cell(mm.tot_yrs,1);
iterH_check.mat_h_yr_sales  = cell(mm.tot_yrs,1);
iterH_check.firm_h_yr_sales = cell(mm.tot_yrs,1);

for t = 2:1:mm.periods
    iterH_in.t = t;
    if mod(t,mm.pd_per_yr) == 0; iterH_in.season = 12;
    else
    iterH_in.season = mod(t,mm.pd_per_yr);
    end       

% The following block labels trans_count output, which are written to
% diagnostics.txt for debugging in simulateMatchesInnerSimMatchTrans
% fileID3 = fopen('results/diagnostics.txt','a');
%    fprintf(fileID3,'\r\n t =%4.0f',t);
%    fprintf(fileID3, '\r\n  ');
% fclose(fileID3);

iterH_in.year = floor((iterH_in.t-1)/mm.pd_per_yr);

[iterH_in] = simulateHomeMatchesInnerSimClientCounts(iterH_in, mm, policy);
[iterH_in] = simulateHomeMatchesInnerSimUpdZHotel(iterH_in, mm, policy);
[iterH_in] = simulateHomeMatchesInnerSimKickDormant(iterH_in, mm);
[iterH_in] = simulateHomeMatchesInnerSimFirmAge(iterH_in, mm);
[iterH_in] = simulateHomeMatchesInnerSimMatchLevelData(iterH_in, mm);


    
 mat_tran = iterH_in.mat_tran;
 iterH_in.N_match = iterH_in.N_match + size(mat_tran,1);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% construct annualized variables
    if iterH_in.season == mm.pd_per_yr

  [iterH_in.mat_h_yr_sales,iterH_in.firm_h_yr_sales,iterH_in.Zcut_eoy] = season_merge(iterH_in,mm);

        % firm_h_yr_sales:[firm ID, total dom. sales, total dom. shipments, firm age in domestic market]
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% the following matrices accumulate annualized values over time and firm types

        theta_h_firm = iterH_in.theta_h(floor(iterH_in.firm_h_yr_sales(:,1)));
        ttt = ones(size(iterH_in.firm_h_yr_sales,1),1).*[t,iterH_in.ptm_type];
        iter_out.firm_h_yr_sales = [iter_out.firm_h_yr_sales;[ttt,iterH_in.firm_h_yr_sales]];
        iter_out.theta_h_firm  = [iter_out.theta_h_firm;theta_h_firm]; % keep track of domestic thetas for each firm
        % iter_out.firm_h_yr_sales: [t,type,firm ID, total sales, total shipments,firm age]

        if pt_ndx == mm.check_type
        yr_ndx = iterH_in.year+1;
        iterH_check.seas_tran{yr_ndx}       = iterH_in.seas_tran;
        iterH_check.match_mat{yr_ndx}       = iterH_in.mat_tran;
        iterH_check.Zcut_eoyH{yr_ndx}       = iterH_in.Zcut_eoy; 
        iterH_check.mat_h_yr_sales{yr_ndx}  = iterH_in.mat_h_yr_sales;
        iterH_check.firm_h_yr_sales{yr_ndx} = iterH_in.firm_h_yr_sales;
        end
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Construct and cumulate moments
           
        
        if iterH_in.year > mm.burn  % cosntruct moments for firm domestic sales regressions
            % autoregressions and degree distribution
            
         iterH_in.mat_yr_sales = iterH_in.mat_h_yr_sales;
         iterH_in.mat_yr_sales_lag = iterH_in.mat_h_yr_sales_lag;

         [iterH_in.mat_h_cont_2yr,iterH_in.mat_h_yr_sales,iterH_in.mat_h_yr_sales_lag,~] =...
            mat_yr_splice_v2(iterH_in,mm,iterH_in.year);  
   
            [x,y,fmoms_h_xx,fmoms_h_xy,fysum_h,fn_obs_h] = firm_reg_h_moms(iterH_in,mm);

            iter_out.x_fsales_h   = [iter_out.x_fsales_h;x];
            iter_out.y_fsales_h   = [iter_out.y_fsales_h;y];
            iter_out.fmoms_h_xx = iter_out.fmoms_h_xx + fmoms_h_xx; % cumulate moments for home sales AR1
            iter_out.fmoms_h_xy = iter_out.fmoms_h_xy + fmoms_h_xy; % cumulate moments for home sales AR1
            iter_out.fysum_h    = iter_out.fysum_h + fysum_h;   
            iter_out.fnobs_h    = iter_out.fnobs_h + fn_obs_h ;

        end   % year > mm.burn if statement
        iterH_in.firm_h_yr_sales_lag = iterH_in.firm_h_yr_sales; % stack data for firm regression
        iterH_in.mat_h_yr_sales_lag  = iterH_in.mat_h_yr_sales;   
        iterH_in.Zcut_eoy_lag        = iterH_in.Zcut_eoy;        
        
        
    end   % season == mm.pd_per_yr if statement
    
    iterH_in.keep_cli_lag        = iterH_in.keep_cli;
 
    %% load lagged client state matrix and re-initialize objects

    iterH_in.lag_cli_zst  = iterH_in.cur_cli_zst;
    iterH_in.new_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iterH_in.die_cli_zst  = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iterH_in.trans_zst    = zeros(mm.sim_firm_num_by_prod_succ_type(pt_ndx),size(mm.Z,1));
    iterH_in.trans_count  = zeros(size(mm.Z,1)+1,size(mm.Z,1)+1,mm.sim_firm_num_by_prod_succ_type(pt_ndx));

 %     fprintf('\r passing line 106 simulateHomeMatchesInnerSim, period %.0f\n', t);
if t == mm.periods   
%   for checking only: collect count series for each firm type
    find_hcli        = find(sum(iterH_in.cur_cli_cnt,2)>0);
    iter_out.transH{pt_ndx,1} = find_hcli;
    iter_out.transH{pt_ndx,2} = iterH_in.cur_cli_cnt(find_hcli,:);
    iter_out.transH{pt_ndx,3} = iterH_in.cum_succ(find_hcli,:);
    iter_out.transH{pt_ndx,4}  = iterH_in.cumage(find_hcli,:);  
    iter_out.transH{pt_ndx,5}  = iterH_in.new_firm(find_hcli,:);  
    
    
%% for checking only: rearrange count series in blocks, by firm t and firm #
%     rooms =  iter_out.transH{pt_ndx,1};
%     stackH = zeros(length(rooms),mm.periods+1); 
%     for i=1:length(rooms)
%     lb = (i-1)*4 + 1;
%     ub = i*4;
%     stackH(lb:ub,:) = ... 
%     [rooms(i),iter_out.transH{pt_ndx,2}(i,:);...
%     rooms(i),iter_out.transH{pt_ndx,3}(i,:);...
%     rooms(i),iter_out.transH{pt_ndx,4}(i,:);... 
%     rooms(i),iter_out.transH{pt_ndx,5}(i,:)];
%     iter_out.stackH = stackH;
%           
%     end
% 
%    iter_out.domfirm_count = length(rooms);
 %% end checking block 
 
    iter_out.iterH_check = iterH_check;

end


end


end

function [iter_in,iter_out] = simulateForeignMatchesInnerMoments(iter_in,iter_out,mm)
if iter_in.year > mm.burn  % don't start building simulated data set until burn-in finished
    if size(iter_in.mat_cont_2yr)>0

        tt =  ones(size(iter_in.mat_yr_sales,1),1).*[iter_in.t,iter_in.mic_type]; % add cols 1 and 2: t, firm type
        iter_out.mat_yr_sales  = [iter_out.mat_yr_sales;[tt,iter_in.mat_yr_sales]];
       % agg_mat_yr_sales: [t,type,firm ID, match sales, shipments, boy Z, eoy Z, w/in yr. match age, firm age]

         ttt = ones(size(iter_in.firm_yr_sales,1),1).*[iter_in.t,iter_in.mic_type];
         iter_out.firm_f_yr_sales = [iter_out.firm_f_yr_sales;[ttt,iter_in.firm_yr_sales]]; % add cols 1 and 2: t, firm type
       % agg_firm_yr_sales: [t,type,firm ID, total exports,total shipments,firm age]

        iter_out.mat_matur =  [iter_out.mat_matur; iter_in.mat_matur_dat];
%       iter_in.mat_matur_dat: [sales, boy Z, eoy Z, match age, firm age, firm_ID, yr]
%       iter_in.mat_cont_2yr:  [firm_ID, sales, shipments, boy Z, eoy Z, match age, firm age] x 2 (lagged, then current) 

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  match regression moments

        [mat_ar1_x,mat_ar1_y,moms_xx,moms_xy,ysum,n_obs] = match_reg_moms(iter_in.mat_cont_2yr,iter_in.ncols);

        iter_out.moms_xx = iter_out.moms_xx + moms_xx; % cumulate moments for match regression
        iter_out.moms_xy = iter_out.moms_xy + moms_xy; % cumulate moments for match regression
        iter_out.ysum    = iter_out.ysum + ysum;
        iter_out.nobs    = iter_out.nobs + n_obs;

        iter_out.mat_ar1_x = [iter_out.mat_ar1_x;mat_ar1_x];
        iter_out.mat_ar1_y = [iter_out.mat_ar1_y;mat_ar1_y];
    end

end

if iter_in.year >= mm.burn  
    
    % autoregressions and degree distribution

    [fmoms_xx,fmoms_xy,fysum,fn_obs] = firm_reg_moms(iter_in,mm);
% JT: I don't think we're using this regression and in any case it is 
% done using data at the match level (iter_in.mat_cont_2yr), not the firm level

    iter_out.fmoms_xx = iter_out.fmoms_xx + fmoms_xx; % cumulate moments for firm regression
    iter_out.fmoms_xy = iter_out.fmoms_xy + fmoms_xy; % cumulate moments for firm regression
    iter_out.fysum    = iter_out.fysum + fysum;
    iter_out.fnobs    = iter_out.fnobs + fn_obs ;

    % foreign market exit regression moments

    ff_exit = iter_in.mkt_exit(:,2)>0;  % column 2 of mkt_exit is number of meetings
    if sum(ff_exit,1)>0
        [exit_moms_xx,exit_moms_xy,sum_succ_rate,sum_exits,exit_obs] = mkt_exit_moms(iter_in.mkt_exit);

        iter_out.exit_xx = iter_out.exit_xx + exit_moms_xx;
        iter_out.exit_xy = iter_out.exit_xy + exit_moms_xy';
        iter_out.sum_succ_rate = iter_out.sum_succ_rate + sum_succ_rate;
        iter_out.exit_obs = iter_out.exit_obs + exit_obs;
        iter_out.sum_exits = iter_out.sum_exits + sum_exits;
    end
    % match exit regression moments
    if iter_in.year_lag == iter_in.year
        ff_mexit = iter_in.mat_yr_sales(:,2)>0; % positive lagged exports
        if sum(ff_mexit,1)>0 
            
            [mat_exit_x,mat_exit_y,mat_exit_moms_xx,mat_exit_moms_xy,mat_obs,nmat_exit]...
                = match_exit_moms(iter_in.mat_yr_sales(ff_mexit,:),mm.pd_per_yr);

            % Notes on variables:

            %         mat_exit_y  = matches(ff,5)==0;         % match dead by end of year
            %         mat_exit_x = [x0,x1,x2,x3,x4];
            %           x0 = ones(size(ff,1),1);              % ff picks off non-missing obs.
            %           x1 = matches(ff,4)==0;                % first year dummy
            %           x2 = log(matches(ff,2));              % sales during year
            %           x3 = log(1+matches(ff,6)./mm.pd_per_yr); % age of match
            %           x4 = log(1+matches(ff,7));            % age of exporter

            iter_out.mat_exit_moms_xx = iter_out.mat_exit_moms_xx + mat_exit_moms_xx;
            iter_out.mat_exit_moms_xy = iter_out.mat_exit_moms_xy + mat_exit_moms_xy;
            iter_out.mat_obs      = iter_out.mat_obs + mat_obs;
            iter_out.nmat_exit    = iter_out.nmat_exit + nmat_exit;
            iter_out.mat_exit_x   = [iter_out.mat_exit_x;mat_exit_x];
            iter_out.mat_exit_y   = [iter_out.mat_exit_y;mat_exit_y];
        end
    end
    % shipment and match counter
    
    [nship_obs,ln_ships,match_count] = match_shpt_cntr(iter_in.mat_yr_sales,mm.max_match);

    iter_out.ship_obs    = iter_out.ship_obs + nship_obs ;
    iter_out.ln_ships    = iter_out.ln_ships + ln_ships ;
    
   try
   if size(match_count,1)>0
     match_histogram = sum(match_count*ones(1,mm.max_match) - ones(size(match_count,1),1)*(1:mm.max_match)==0);
     iter_out.match_count = iter_out.match_count + match_histogram ;
   end
   catch 
        'pause in simulateForeignMatchesInnerSim line 96'
   end

end  

iter_in.mat_yr_sales_lag = iter_in.mat_yr_sales;   % update lags
iter_in.firm_yr_sales_lag = iter_in.firm_yr_sales;

end
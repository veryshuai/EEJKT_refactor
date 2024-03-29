function [iter_in,iter_out] = simulateForeignMatchesInnerMoments(iter_in,iter_out,mm)
if iter_in.year > mm.burn  % don't start building simulated data set until burn-in finished
    
        tt =  ones(size(iter_in.mat_yr_sales,1),1).*[iter_in.t,iter_in.mic_type]; % add cols 1 and 2: t, firm type
        iter_out.mat_yr_sales  = [iter_out.mat_yr_sales;[tt,iter_in.mat_yr_sales]];    
    
    if size(iter_in.mat_cont_2yr)>0

%        tt =  ones(size(iter_in.mat_yr_sales,1),1).*[iter_in.t,iter_in.mic_type]; % add cols 1 and 2: t, firm type
%        iter_out.mat_yr_sales  = [iter_out.mat_yr_sales;[tt,iter_in.mat_yr_sales]];
       % mat_yr_sales: [t,type,firm ID, match sales, shipments, boy Z, eoy Z, match age, firm age]

         ttt = ones(size(iter_in.firm_yr_sales,1),1).*[iter_in.t,iter_in.mic_type];
         iter_out.firm_f_yr_sales = [iter_out.firm_f_yr_sales;[ttt,iter_in.firm_yr_sales]]; % add cols 1 and 2: t, firm type
       % firm_f_yr_sales: [t,type,firm ID, total exports,total shipments,firm age]

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%  match AR1 regression moments

        [mat_ar1_x,mat_ar1_y,moms_xx,moms_xy,ysum,n_obs] = match_reg_moms(iter_in.mat_cont_2yr,iter_in.ncols,mm);
%       iter_in.mat_cont_2yr:  [firm_ID, sales, shipments, boy Z, eoy Z, match age, firm age] x 2 (lagged, then current) 

        iter_out.moms_xx = iter_out.moms_xx + moms_xx; % cumulate moments for match regression
        iter_out.moms_xy = iter_out.moms_xy + moms_xy; % cumulate moments for match regression
        iter_out.ysum    = iter_out.ysum + ysum;
        iter_out.nobs    = iter_out.nobs + n_obs;

        iter_out.mat_ar1_x = [iter_out.mat_ar1_x;mat_ar1_x]; % cumulate explanatory variable matrix
        iter_out.mat_ar1_y = [iter_out.mat_ar1_y;mat_ar1_y]; % cumulate dependent variable
        
        % NOTE: we don't really need to cumulate both the moments and the raw data
    end

end

if iter_in.year >= mm.burn  
    
    
    %% match exit regression moments
    if iter_in.year_lag == iter_in.year
        ff_mexit = iter_in.mat_yr_sales(:,2)>0; % positive exports
        if sum(ff_mexit,1)>0 
            
            [mat_exit_x,mat_exit_y,mat_exit_moms_xx,mat_exit_moms_xy,mat_obs,nmat_exit]...
                = match_exit_moms(iter_in.mat_yr_sales(ff_mexit,:),mm.pd_per_yr);

            % Notes on variables:

            %  mat_exit_y  = matches(ff,5)==0;         % match dead by end of year
            %  mat_exit_x = [x0,x1,x2,x3,x4];
            %    x0 = ones(size(ff,1),1);                 % ff picks off non-missing obs.
            %    x1 = matches(ff,4)==0;                   % match new this year
            %    x2 = log(matches(ff,2));                 % log sales during year
            %    x3 = log(1+matches(ff,6)./mm.pd_per_yr); % log age of match (in years)
            %    x4 = log(1+matches(ff,7)./mm.pd_per_yr); % log age of exporter (in years)
                       

            iter_out.mat_exit_moms_xx = iter_out.mat_exit_moms_xx + mat_exit_moms_xx;
            iter_out.mat_exit_moms_xy = iter_out.mat_exit_moms_xy + mat_exit_moms_xy;
            iter_out.mat_obs      = iter_out.mat_obs + mat_obs;
            iter_out.nmat_exit    = iter_out.nmat_exit + nmat_exit;
            iter_out.mat_exit_x   = [iter_out.mat_exit_x;mat_exit_x];
            iter_out.mat_exit_y   = [iter_out.mat_exit_y;mat_exit_y];
        end
    end
    %% shipment and match counter
    
    [nship_obs,ln_ships,match_count,match_countD,dud_matches] = match_shpt_cntr(iter_in,mm);
  % dud matches: [firm_ID, sales at Zcut, shipments (1), bop Z = eop Z = match_age = 0, firm age]
  % match_count and match_countD are vectors of match counts for active firms w/out % w/ duds
  % ln_ships is the sum of ln # shipments across matches; nships_obs is
  % number of matches with positive shipments.
  
    iter_out.ship_obs    = iter_out.ship_obs + nship_obs ;
    iter_out.ln_ships    = iter_out.ln_ships + ln_ships ;
    
    ttt = ones(size(dud_matches,1),1).*[iter_in.t,iter_in.mic_type];
    iter_out.dud_matches = [iter_out.dud_matches; [ttt,dud_matches]] ;
  % iter_out.dud matches: [t, mic_type, firm_ID, sales at Zcut, shipment=1,  bop Z = eop Z = match_age = 0, firm age]
    
    % topcode match counts for each firm_ID 
    match_count  = min(mm.max_match, match_count);
    match_countD = min(mm.max_match, match_countD);
    
    % convert to frequency counts (# firms with each possible match count)
   try
   if size(match_count,1)>0 % counts excluding duds
     match_histogram = sum(match_count*ones(1,mm.max_match) - ones(size(match_count,1),1)*(1:mm.max_match)==0,1);
     iter_out.match_hist = iter_out.match_hist + match_histogram ;
   end
   
   if size(match_countD,1)>0 % counts including duds (alternative meaasure)
     match_histogramD = sum(match_countD*ones(1,mm.max_match) - ones(size(match_countD,1),1)*(1:mm.max_match)==0,1);
     iter_out.match_histD = iter_out.match_histD + match_histogramD ;
   end  
   
%   iter_out.duds = iter_out.duds + sum(duds);
   
   catch 
        'problem in simulateForeignMatchesInnerMoments lines 93-101'
   end

end  
%% update lags

iter_in.mat_yr_sales_lag  = iter_in.mat_yr_sales;    
iter_in.firm_yr_sales_lag = iter_in.firm_yr_sales;
iter_in.Zcut_eoy_lag      = iter_in.Zcut_eoy;

end
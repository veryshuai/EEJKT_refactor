% clear
  seed = 1; % reset random number generator
    if seed == 1
        rng(80085,'twister');
    end

% control size of simulated dataset
tot_sims  = mm.S; % total number of firms to simulate
tot_yrs   = mm.tot_yrs; % total number of years to simulate
pd_per_yr = mm.pd_per_yr;       % periods per year
frac_of_year = 1/pd_per_yr;
periods   = round(tot_yrs*pd_per_yr); % number of periods to simulate

max_macro_shk = size(mm.X_f,1);
N_types       = size(typemat,1);
N_mic_types   = N_types/max_macro_shk; % number of exog. micro types (N_phi x N_theta)
N_mic_state   = size(pmat_cum_f{1},2); % number of micro states for any given type

% adjust all hazards for period length
% prob_mdeath = 1-exp(-frac_of_year * mm.delta);  % per period probability of match death 
% prob_fdeath = 1-exp(-frac_of_year * mm.d);      % per period probability of firm death 
% haz_ship    = frac_of_year * mm.L_b;            % per period shipment hazard

% convert hazards to per-period probabilities
prob_mdeath = 1-exp(-mm.delta);  % per period probability of match death 
prob_fdeath = 1-exp(-mm.d);      % per period probability of firm death 

% construct CDF for within-period shipments of active matches
haz_ship  =  mm.L_b;           % shipment hazard
max_ships = 3*round(haz_ship); % maximum within-period shipments is triple expected number
poisCDF   = poisscdf(1:1:max_ships,haz_ship);

% create CDF for match shocks (for later use inside loop)
cum_pz    = cumsum(mm.erg_pz); 

% create macro trajectory
macro_state_f    = zeros(periods,1);
macro_state_f(1) = 8; %  start macro trajectory at midpoint of distribution
macro_state_h    = zeros(periods,1);
macro_state_h(1) = 8; %  start macro trajectory at midpoint of distribution
for t = 2:periods
    macro_state_f(t) = find(pmat_cum_msf(macro_state_f(t-1),:)>rand(1,1),1,'first'); % update macro state
    macro_state_h(t) = find(pmat_cum_msh(macro_state_h(t-1),:)>rand(1,1),1,'first'); % update macro state
end

%% Create objects that will store firm-type specific simulated data
% pt_type = [kron((1:N_Phi)',ones(N_theta2,1)),kron(ones(N_Phi,1),(1:N_theta2)')];
N_pt    = N_Phi*N_theta2;

s_mat_yr_sales    = cell(N_pt,1); % for analysis of match dynamics
s_mat_yr_sales_adj= cell(N_pt,1); % for analysis of match exits
s_mat_ar1_x       = cell(N_pt,1); % for match ar1 regression residuals
s_mat_ar1_y       = cell(N_pt,1); % for match ar1 regression residuals
s_mat_exit_x      = cell(N_pt,1); % for match exit regression residuals
s_mat_exit_y      = cell(N_pt,1); % for match exit regression residuals  
s_x_hf            = cell(N_pt,1); % for home-foreign sales regression residuals
s_y_hf            = cell(N_pt,1); % for home-foreign sales regression residuals
s_x_fsales_h      = cell(N_pt,1); % for home sales AR1 residuals
s_y_fsales_h      = cell(N_pt,1); % for home sales AR1 residuals
s_time_gaps       = cell(N_pt,1); % for match hazard analysis

firm_cntr = zeros(N_Phi,N_theta2);  % simulated firm counter, by type

% match level moment aggregators
s_moms_xx   = zeros(N_pt,4,4);
s_moms_xy   = zeros(N_pt,4);
s_ysum      = zeros(N_pt,1);
s_nobs      = zeros(N_pt,1);
s_ship_obs  = zeros(N_pt,1);
s_ln_ships = zeros(N_pt,1);

% firm level moment aggregators
s_fmoms_xx = zeros(N_pt,4,4);
s_fmoms_xy = zeros(N_pt,4);
s_fmoms_h_xx = zeros(N_pt,2,2);
s_fmoms_h_xy = zeros(N_pt,2);
s_fysum    = zeros(N_pt,1);
s_fnobs    = zeros(N_pt,1);
s_fysum_h  = zeros(N_pt,1);
s_fnobs_h  = zeros(N_pt,1); 

s_exit_xx       = zeros(N_pt,6,6);
s_exit_xy       = zeros(N_pt,6);
s_sum_succ_rate = zeros(N_pt,1);
s_exit_obs      = zeros(N_pt,1);
s_sum_exits     = zeros(N_pt,1);

% home-foreign firm level moment aggregators
s_hfmoms_xx = zeros(N_pt,2,2);
s_hfmoms_xy = zeros(N_pt,2);
s_hfysum    = zeros(N_pt,1); 
s_hf_nobs   = zeros(N_pt,1);  
s_nfirm     = zeros(N_pt,1); 
s_nexptr    = zeros(N_pt,1); 
s_expt_rate = cell(N_pt);  

% match death regressions
s_mat_exit_moms_xx  = zeros(N_pt,5,5);
s_mat_exit_moms_xy  = zeros(N_pt,5);
s_mat_obs       = zeros(N_pt,1); 
s_nmat_exit     = zeros(N_pt,1); 

% match counter for degree distribution
 max_match = 50; % upper bound on number of matches to be counted 
 s_match_count = zeros(N_pt,max_match);
 s_singletons  = zeros(N_pt,1);
 NN            = zeros(N_pt,1);

 %% Load storage objects by looping over firm types

N_theta1 = size(mm.theta1,2);    % number of thetas in home market
th1_cdf  = betacdf(mm.theta1,mm.ah,mm.bh); % cdf for home theta draws
th1_cdf(N_theta1) = 1; % to deal with rounding problem

%seed_cntr = 1:1:N_pt;
seeds = randi(1e6,N_Phi,2);
tic3 = tic;
too_slow = zeros(N_pt,1); % flags when jumped out of matdat_gen loops

pt_type = [kron((1:N_Phi)',ones(N_theta2,1)),kron(ones(N_Phi,1),(1:N_theta2)')];

parfor pt_ndx = 1:1:N_pt
%   for pt_ndx = 1:1:N_pt % use this for loop for debugging only

prod_ndx  = pt_type(pt_ndx,1);
theta_ndx = pt_type(pt_ndx,2);

if seed == 1
 rng(seeds(prod_ndx,1),'twister');
 seed_crand(seeds(prod_ndx,2));
end

  prod_lvl  = mm.Phi(prod_ndx);
  succ_prob = mm.theta2(theta_ndx);
  
  % number of firms to simulate of this particular typec_val_f_origc_val_f_origc_val_f_origc_val_f_origc_val_f_orig
  N_firms = round(mm.erg_pp(prod_ndx)*th2_pdf(theta_ndx)*tot_sims);
  NN(pt_ndx) = N_firms;
    
   if N_firms>0
       
%   if pt_ndx == 102
%       'pause here'
%   end
%%                         
  [s_singletons(pt_ndx),s_time_gaps{pt_ndx},s_mat_yr_sales{pt_ndx},s_mat_yr_sales_adj{pt_ndx},...
   firm_f_yr_sales,s_mat_ar1_x{pt_ndx},s_mat_ar1_y{pt_ndx},s_moms_xx(pt_ndx,:,:),...
   moms_xy,s_ysum(pt_ndx),s_nobs(pt_ndx),s_fmoms_xx(pt_ndx,:,:),fmoms_xy,... 
   s_fysum(pt_ndx),s_fnobs(pt_ndx),s_exit_xx(pt_ndx,:,:),exit_xy,s_sum_succ_rate(pt_ndx),s_exit_obs(pt_ndx),...
   s_sum_exits(pt_ndx),s_mat_exit_moms_xx(pt_ndx,:,:),mat_exit_moms_xy,...
     s_mat_exit_x{pt_ndx} ,s_mat_exit_y{pt_ndx},s_mat_obs(pt_ndx),s_nmat_exit(pt_ndx),...
    s_ship_obs(pt_ndx),s_ln_ships(pt_ndx), s_match_count(pt_ndx,:),abort_flag_f]...
      = matchdat_gen_f(N_firms,N_types,N_Z,N_Phi,N_mic_state,pd_per_yr,frac_of_year,periods,max_ships,typemat,macro_state_f,...
              theta_ndx,prod_ndx,mm,pmat_cum_f,c_val_f_orig,succ_prob,prod_lvl,cum_pz,poisCDF,...
              prob_mdeath,prob_fdeath,haz_ship,pmat_cum_z,nn2,nn4,lambda_f,max_match);
          
   s_fmoms_xy(pt_ndx,:) = fmoms_xy';
   s_exit_xy(pt_ndx,:,:) = exit_xy';     
   s_mat_exit_moms_xy(pt_ndx,:) = mat_exit_moms_xy';
   s_moms_xy(pt_ndx,:) = moms_xy';

 % s_mat_yr_sales{}: [firm ID, match-specific sales, shipments, boy Z, eoy Z, match age in periods (w/in year), firm age in periods]   
 % firm_f_yr_sales: [t,type,firm ID, total exports,total foreign shipments,firm age in export mkt.]      
 % s_mat_matur{}: [sales, boy Z, eoy Z, match age, firm age, firm_ID, yr]  
 
   [firm_h_yr_sales,s_x_fsales_h{pt_ndx},s_y_fsales_h{pt_ndx},theta_h_firm,s_fmoms_h_xx(pt_ndx,:,:),...
      fmoms_h_xy,s_fysum_h(pt_ndx), s_fnobs_h(pt_ndx),abort_flag_h]...
      = matchdat_gen_h(N_firms,N_Z,N_Phi,nn_h,pd_per_yr,frac_of_year,periods,max_ships,typemat,macro_state_h,...
              theta_ndx,prod_ndx,mm,pmat_cum_h,c_val_h_orig,cum_pz,poisCDF,...
              prob_mdeath,prob_fdeath,haz_ship,pmat_cum_z,N_theta1,th1_cdf);
 %%         
 
% Alternative code for domestic sales statistics doesn't use intensity
% matrix. It's much slower, and the results are pretty similar:
% 
%   [firm_h_yr_sales,s_x_fsales_h{pt_ndx},s_y_fsales_h{pt_ndx},theta_h_firm,s_fmoms_h_xx(pt_ndx,:,:),...
%    fmoms_h_xy,s_fysum_h(pt_ndx), s_fnobs_h(pt_ndx),abort_flag_h]...
%       = matchdat_gen_h2(N_firms,N_Z,pd_per_yr,frac_of_year,periods,max_ships,typemat,macro_state_h,...                     
%               theta_ndx,prod_ndx,mm,c_val_h_orig,succ_prob,cum_pz,poisCDF,...
%               prob_mdeath,prob_fdeath,pmat_cum_z,nn2,lambda_h,N_theta1,th1_cdf,max_match);             
      
  s_fmoms_h_xy(pt_ndx,:) = fmoms_h_xy';        
 
  
 % firm_h_yr_sales: [t,type,firm ID, total dom. sales, total dom. shipments,firm age in dom. mkt.]

  too_slow(pt_ndx,1) = abort_flag_f + abort_flag_h;

% splice home and foreign sales and create moments for home-foreign regression

% if N_firms>=20  && sum(sum(firm_f_yr_sales(:,4)>0))>20 
%     pt_ndx
%     N_firms
%     'pause here'    
% end

  [s_x_hf{pt_ndx}, s_y_hf{pt_ndx},s_expt_rate{pt_ndx},s_nfirm(pt_ndx),...
      s_nexptr(pt_ndx),nhfirms,hfmoms_xy,s_hfmoms_xx(pt_ndx,:,:),...
      s_hfysum(pt_ndx),s_hf_nobs(pt_ndx)]...
     = splice_hf(firm_h_yr_sales,firm_f_yr_sales,typemat,theta_h_firm,periods);
 
 s_hfmoms_xy(pt_ndx,:) = hfmoms_xy';
 
   end
 
end        % end of pt_ndx loop 
 
time3 = toc(tic3);
fprintf(' discrete_sim big loop time: %.2f\n', time3);
%% Initialize objects that will hold cumulated values


% create first observation on firm-year level aggregates (will concatenate below)
agg_mat_yr_sales    = zeros(0,9); % for analysis of match dynamics
agg_mat_yr_sales_adj= zeros(0,9); % for analysis of match exit
agg_mat_ar1_x       = zeros(0,4); % for match ar1 regression residuals
agg_mat_ar1_y       = zeros(0,1); % for match ar1 regression residuals
agg_mat_exit_x      = zeros(0,5); % for match exit regression residuals
agg_mat_exit_y      = zeros(0,1); % for match exit regression residuals
agg_mat_matur       = zeros(0,8); % for match maturation analysis
agg_x_hf            = zeros(0,2); % for home-foreign sales regression residuals
agg_y_hf            = zeros(0,1); % for home-foreign sales regression residuals
agg_x_fsales_h      = zeros(0,2); % for home sales AR1 residuals
agg_y_fsales_h      = zeros(0,1); % for home sales AR1 residuals
agg_time_gaps       = zeros(0,7); % for match hazard analysis

% firm_cntr = 0;  % simulated firm counter, all types combined
chk1 = 0;
chk2 = 0;

% match level moment aggregators
agg_moms_xx   = zeros(4,4);
agg_moms_xy   = zeros(4,1);
agg_ysum      = 0;
agg_nobs      = 0;
agg_ship_obs  = 0;
agg_ln_ships = 0;

% firm level moment aggregators
agg_fmoms_xx = zeros(4,4);
agg_fmoms_xy = zeros(4,1);
agg_fmoms_h_xx = zeros(2,2);
agg_fmoms_h_xy = zeros(2,1);
agg_fysum    = 0;
agg_fnobs    = 0; 
agg_fysum_h  = 0;
agg_fnobs_h  = 0; 

agg_exit_xx = zeros(6,6);
agg_exit_xy = zeros(6,1);
agg_sum_succ_rate = 0;
agg_exit_obs = 0;
agg_sum_exits = 0;

% home-foreign firm level moment aggregators
agg_hfmoms_xx = zeros(2,2);
agg_hfmoms_xy = zeros(2,1);
agg_hfysum    = 0;
agg_hf_nobs    = 0; 
agg_nfirm     = 0;
agg_nexptr    = 0;
agg_expt_rate = zeros(0,1);  

% match death regressions
agg_mat_exit_moms_xx  = zeros(5,5);
agg_mat_exit_moms_xy  = zeros(5,1);
agg_mat_obs       = 0;
agg_nmat_exit     = 0;
 
agg_match_count  = zeros(max_match,1);
singletons = 0;

%% Cumulate over firm types

for pt_ndx = 1:1:N_pt
 
  if NN(pt_ndx) > 0
     % cumulate time gaps for match hazard regressions

        agg_time_gaps = [agg_time_gaps;s_time_gaps{pt_ndx}];
      
% ccumulate annualized values over time and firm types 

         agg_mat_yr_sales  = [agg_mat_yr_sales;s_mat_yr_sales{pt_ndx}]; 
         agg_mat_yr_sales_adj  = [agg_mat_yr_sales_adj;s_mat_yr_sales_adj{pt_ndx}]; 
         
       % agg_mat_yr_sales: [t,type,firm ID, match sales, shipments, boy Z, eoy Z, match age, firm age] 
         agg_x_hf = [agg_x_hf;s_x_hf{pt_ndx}];
         agg_y_hf = [agg_y_hf;s_y_hf{pt_ndx}];
       % agg_x_hf and agg_y_hf are data for the home-foreign regression 

         agg_x_fsales_h = [agg_x_fsales_h;s_x_fsales_h{pt_ndx}];
         agg_y_fsales_h = [agg_y_fsales_h;s_y_fsales_h{pt_ndx}];
       % agg_x_fsales_h and agg_y_fsales_h are obs. for the home sales AR1   
      

   %  match regression moments
         agg_moms_xx = agg_moms_xx + squeeze(s_moms_xx(pt_ndx,:,:)); % cumulate moments for match regression
%        agg_moms_xy = agg_moms_xy + squeeze(s_moms_xy(pt_ndx,:)); % cumulate moments for match regression
         agg_moms_xy = agg_moms_xy + s_moms_xy(pt_ndx,:)'; % cumulate moments for match regression
         agg_ysum    = agg_ysum    + s_ysum(pt_ndx);
         agg_nobs    = agg_nobs    + s_nobs(pt_ndx);   
                
         agg_mat_ar1_x = [agg_mat_ar1_x;s_mat_ar1_x{pt_ndx}];
         agg_mat_ar1_y = [agg_mat_ar1_y;s_mat_ar1_y{pt_ndx}];
         
       % agg_mat_ar1_x and agg_mat_ar1_y are data for the home-foreign regression
 
   % firm level autoregressions
         agg_fmoms_xx = agg_fmoms_xx + squeeze(s_fmoms_xx(pt_ndx,:,:)); % cumulate moments for firm regression, foreign
  %      agg_fmoms_xy = agg_fmoms_xy + squeeze(s_fmoms_xy(pt_ndx,:)); % cumulate moments for firm regression, foreign
         agg_fmoms_xy = agg_fmoms_xy + squeeze(s_fmoms_xy(pt_ndx,:)'); % cumulate moments for firm regression, foreign
         
         agg_fysum    = agg_fysum + s_fysum(pt_ndx);
         agg_fnobs    = agg_fnobs + s_fnobs(pt_ndx) ; 
         
         agg_fmoms_h_xx = agg_fmoms_h_xx + squeeze(s_fmoms_h_xx(pt_ndx,:,:)); % cumulate moments for firm regression, home
%        agg_fmoms_h_xy = agg_fmoms_h_xy + squeeze(s_fmoms_h_xy(pt_ndx,:)); % cumulate moments for firm regression, home
         agg_fmoms_h_xy = agg_fmoms_h_xy + squeeze(s_fmoms_h_xy(pt_ndx,:)'); % cumulate moments for firm regression, home
 

         agg_fysum_h    = agg_fysum_h + s_fysum_h(pt_ndx);
         agg_fnobs_h    = agg_fnobs_h + s_fnobs_h(pt_ndx) ; 
            
    % foreign-home firm level autoregressions and export rates
         agg_hfmoms_xx = agg_hfmoms_xx + squeeze(s_hfmoms_xx(pt_ndx,:,:));
%        agg_hfmoms_xy = agg_hfmoms_xy + squeeze(s_hfmoms_xy(pt_ndx,:));
         agg_hfmoms_xy = agg_hfmoms_xy + squeeze(s_hfmoms_xy(pt_ndx,:)');
         agg_hfysum    = agg_hfysum    + s_hfysum(pt_ndx);
         agg_hf_nobs   = agg_hf_nobs   + s_hf_nobs(pt_ndx);
         
         agg_nfirm     = agg_nfirm + s_nfirm(pt_ndx);
         agg_nexptr    = agg_nexptr + s_nexptr(pt_ndx);
         agg_expt_rate = [agg_expt_rate;s_expt_rate{pt_ndx}];        
         
    % foreign market exit regression moments
         agg_exit_xx = agg_exit_xx + squeeze(s_exit_xx(pt_ndx,:,:));
   %     agg_exit_xy = agg_exit_xy + squeeze(s_exit_xy(pt_ndx,:));
         agg_exit_xy = agg_exit_xy + squeeze(s_exit_xy(pt_ndx,:)');
         agg_sum_succ_rate = agg_sum_succ_rate + s_sum_succ_rate(pt_ndx);
         agg_exit_obs = agg_exit_obs + s_exit_obs(pt_ndx);
         agg_sum_exits = agg_sum_exits + s_sum_exits(pt_ndx);
         
    % match exit regression moments
         agg_mat_exit_moms_xx = agg_mat_exit_moms_xx + squeeze(s_mat_exit_moms_xx(pt_ndx,:,:));
         agg_mat_exit_moms_xy = agg_mat_exit_moms_xy + squeeze(s_mat_exit_moms_xy(pt_ndx,:)');
         agg_mat_obs          = agg_mat_obs + s_mat_obs(pt_ndx);
         
         agg_nmat_exit        = agg_nmat_exit + s_nmat_exit(pt_ndx);   
         agg_mat_exit_x       = [agg_mat_exit_x;s_mat_exit_x{pt_ndx}];
         agg_mat_exit_y       = [agg_mat_exit_y;s_mat_exit_y{pt_ndx}];        
%        Add firm type identifier as first column of s_mat_matur before aggregating across firm types:
        
    % shipment and match counter
        agg_ship_obs    = agg_ship_obs    + s_ship_obs(pt_ndx) ;
        agg_ln_ships    = agg_ln_ships    + s_ln_ships(pt_ndx) ;
        agg_match_count = agg_match_count + squeeze(s_match_count(pt_ndx,:)') ;
        singletons      = singletons + s_singletons(pt_ndx);
        
 
  end

end        % end of prod_ndx loop

%% Construct simulated statistics

% save ('match_count.mat','s_match_count','pi_tilda_f_new','pt_type','NN')
% active_client_value

    simMoms = struct; %container for all simulated moments
        
        % Some numbers from above
        simMoms.agg_nexptr = agg_nexptr;
        simMoms.agg_nfirm = agg_nfirm;

        % match-level autoregression
        if rank(agg_moms_xx) == size(agg_moms_xx,2)
            inv_agg_moms_xx = inv(agg_moms_xx);
%         elseif rank(agg_moms_xx) == size(agg_moms_xx,2) - 1
%             'singular matrix for match-level AR1: rank=3'
%             inv_agg_moms_xx = [inv(agg_moms_xx(1:3,1:3)), zeros(3,1); zeros(1,4)];
        else
            rank_xx = rank(agg_moms_xx);
            fprintf('\r Warning: singular matrix for match-level AR1. Rank: %.1f\n', rank_xx); 
            inv_agg_moms_xx = [inv(agg_moms_xx(1:2,1:2)), zeros(2,2); zeros(2,4)];
        end         
        simMoms.beta_match = inv_agg_moms_xx*agg_moms_xy;
        simMoms.ybar_match = agg_ysum/agg_nobs;
        
        simMoms.mse_match_ar1 = (agg_mat_ar1_y - agg_mat_ar1_x*simMoms.beta_match)'*...
           (agg_mat_ar1_y - agg_mat_ar1_x*simMoms.beta_match)/size(agg_mat_ar1_x,1);
          
        
        % firm-level autoregression
        beta_fsales = inv(agg_fmoms_xx)*agg_fmoms_xy;
        ybar_fsales = agg_fysum/agg_fnobs;
        
        simMoms.beta_fsales_h = inv(agg_fmoms_h_xx)*agg_fmoms_h_xy;
        simMoms.mse_h         = (agg_y_fsales_h - agg_x_fsales_h*simMoms.beta_fsales_h)'*...
                        (agg_y_fsales_h - agg_x_fsales_h*simMoms.beta_fsales_h)/size(agg_y_fsales_h,1);
        simMoms.ybar_fsales_h = agg_fysum_h/agg_fnobs_h ;  
        
        % home-foreign regression and export rates
        simMoms.beta_hfsales = inv(agg_hfmoms_xx)*agg_hfmoms_xy;
        simMoms.mse_hf       = (agg_y_hf - agg_x_hf*simMoms.beta_hfsales)'*...
                       (agg_y_hf - agg_x_hf*simMoms.beta_hfsales)/size(agg_x_hf,1);
        simMoms.ybar_hfsales  = agg_hfysum/agg_hf_nobs;

        simMoms.avg_expt_rate = mean(agg_expt_rate);
        simMoms.share_exptr   = agg_nexptr/agg_nfirm;
        
        % market exit regression        
     if rank(agg_exit_xx) == size(agg_exit_xx,2)
        inv_agg_exit_xx = inv(agg_exit_xx);
     else
        rank_xx = rank(agg_exit_xx);
        fprintf('\r Warning: singular matrix for market exit. Rank: %.1f\n', rank_xx); 
        inv_agg_exit_xx = [inv(agg_exit_xx(1:3,1:3)), zeros(3,3); zeros(3,6)];
     end  
        simMoms.beta_mkt_exit   = inv_agg_exit_xx*agg_exit_xy;
        simMoms.mkt_exit_rate   = agg_sum_exits/agg_exit_obs;
        match_succ_rate = agg_sum_succ_rate/agg_exit_obs;
        
        % match exit regression       
        simMoms.beta_match_exit = inv(agg_mat_exit_moms_xx)*agg_mat_exit_moms_xy;
        simMoms.match_exit_rate = agg_nmat_exit/agg_mat_obs;
        mse_match_exit  = (agg_mat_exit_y - agg_mat_exit_x*simMoms.beta_match_exit)'*...
             (agg_mat_exit_y - agg_mat_exit_x*simMoms.beta_match_exit)/size(agg_mat_exit_y,1);
        
 
% average log #shipments
        simMoms.avg_ln_ships = agg_ln_ships/agg_ship_obs;
        
%       % create variables for analysis of degree distribution
%    
        simMoms.ff_sim_max      = find(cumsum(agg_match_count)./sum(agg_match_count)<1);
        log_compCDF     = log(1 - cumsum(agg_match_count(simMoms.ff_sim_max))./sum(agg_match_count));
        log_matches     = log(1:1:size(simMoms.ff_sim_max,1))';
        xmat            = [ones(size(simMoms.ff_sim_max)),log_matches,log_matches.^2];
        
       % quadratic regression approximating degree distribution
        simMoms.b_degree        = regress(log_compCDF,xmat);
       % linear regression approximating degree distribution 
        xmat_linear     = [ones(size(simMoms.ff_sim_max)),log_matches];
        b_degree_linear = regress(log_compCDF,xmat_linear);
       % nonparametric plot of degree distribution 
%       scatter(log_matches,log_compCDF)
        
        % match hazards
        % agg_time_gaps: (1) firm_ID, (2) periods into interval, (3) time gap,  
        %        (4) # new meetings,(5) t, (6) cum. meetings, (7) cum succeses 
        
      % plot histogram of frequencies for meeting hazards
         agg_time_gaps = agg_time_gaps(2:size(agg_time_gaps,1),:);
%         histogram(agg_time_gaps(:,3))
        
       % create variables for hazard regressions     
        ln_haz = log(1./agg_time_gaps(:,3));
        ln_csucc = log(1+agg_time_gaps(:,7)); 
        ln_meet = log(agg_time_gaps(:,6));
        const = ones(size(ln_haz,1),1);
        ln_succ_rate = log(1+(agg_time_gaps(:,7)./agg_time_gaps(:,6))); 
        
       % success rate regression 
        succ_rate = agg_time_gaps(:,7)./agg_time_gaps(:,6);
        [simMoms.b_succ_rate,~,uu] = regress(succ_rate,[const, ln_meet]);
        usq_succ = uu.^2;
        simMoms.b_usq_succ = regress(usq_succ,[const, ln_meet]);
        
      % translog meeting hazard regression 
        X_haz = [const, ln_csucc, ln_csucc.^2, ln_succ_rate, ln_succ_rate.^2, ln_succ_rate.*ln_csucc];
        simMoms.b_haz = regress(ln_haz,X_haz);
      % regression explaining squared residuals of meeting hazard equation  
        usq_haz = (ln_haz - X_haz*simMoms.b_haz).^2;
        b_haz_usq = regress(usq_haz,[const ln_meet]);
      % means of log hazard rate, success rate, and squared residuals from succ rate and hazard regression 
        means_vec = mean([ln_haz,succ_rate,usq_succ,usq_haz]);
        simMoms.mean_ln_haz    = means_vec(1);
        simMoms.mean_succ_rate = means_vec(2);
        simMoms.mean_usq_succ  = means_vec(3);
        mean_usq_haz   = means_vec(4);









     







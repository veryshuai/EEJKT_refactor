function mm = setModelParameters(X)

mm = struct();

%% technology parameters

mm.param_vec = X;      % carry along parameter vector for diagnostic checks

mm.pd_per_yr = 12;        % number of periods per year
mm.r         = 0.13/mm.pd_per_yr;   % Rate of time preference per period
mm.firm_death_haz = 0.08/mm.pd_per_yr;   % Component of time preference due to exogenous death
mm.delta     = 0.326/mm.pd_per_yr;  % Exogenous match separation rate 
mm.scale_h   = X(2);       % Domestic profit function scale parameter
mm.scale_f   = X(2);       % Export profit function scale parameter (same as home)
mm.eta       = 5;          % Demand elasticity 
mm.gam       = X(7);       % Network effect parameter
mm.cs_h      = exp(X(8));  % Cost scaling parameter, home market
mm.cs_f      = exp(X(11)); % Cost scaling parameter, foreign market
mm.sig_p     = X(9);       %standard deviation of productivity distribution
mm.optimism  = X(12);      %parameter on prior distribution (positive means optimistic, negative pessamisitic)


%% Discretization of state-space
%mm.grid_length   = 2.5;   % number of standard deviations from mean used for discretization
mm.n_size        = 20;    % Maximum number of informative signals per firm 
mm.net_size      = 40;    % maximum number of network effects
mm.z_size        = 7;     % Number of discretized demand shock states (2*n+1) 
mm.phi_size      = 8;     % number of different discretized profit shocks (2*n+1)
mm.x_size        = 7;     % Number of different discretized macro shocks; same for home and foreign (2*n+1)
mm.theta_size    = 20;    % Number of possible market potential values (WAS 51)
mm.dim1          = 7;     % Number of possible theta1 values (specific to home market);
mm.dim2          = 7;     % Number pf possible theta2 values (specific to foreign market);


%% Theta distributions

mm.ah        = X(4)*X(3); % Beta function, home (theta1) success parameter
mm.bh        = X(4)*(1-X(3));% Beta function, home (theta1) failure parameter
mm.af        = mm.ah;     % Beta function, foreign (theta2) success parameter (assume same as home)
mm.bf        = mm.bh;     % Beta function, foreign (theta2) success parameter (assume same as home)
mm.F_f       = exp(X(10));       % cost of maintaining a client- foreign
mm.F_h       = exp(X(1)); % cost of maintaining a client- home 

mm.theta1           = linspace((1/mm.dim1)*0.5,(1 - 1/mm.dim1*0.5),mm.dim1); %was 1/mm.dim1:1/mm.dim1:1;
mm.theta2           = linspace((1/mm.dim2)*0.5,(1 - 1/mm.dim2*0.5),mm.dim2); % was 1/mm.dim2:1/mm.dim2:1;
mm.theta1(mm.dim1)  =  mm.theta1(mm.dim1) - 0.0001;
mm.theta2(mm.dim2)  =  mm.theta2(mm.dim2) - 0.0001;

mm.th1_cdf  = betacdf(mm.theta1,mm.ah,mm.bh); % cdf for home theta draws
mm.th1_cdf(size(mm.theta1,2)) = 1; % higher draws assigned to highest success type
mm.th2_cdf = betacdf(mm.theta2,mm.af,mm.bf); % cdf for foreign theta draws
mm.th2_pdf = [mm.th2_cdf(1),mm.th2_cdf(2:mm.dim1)-mm.th2_cdf(1:mm.dim1-1)]; 


%% Solution parameters
mm.v_tolerance   = 1e-3;    % convergence tolerance, value function iterations (WAS 1e-3)
mm.max_iter      = 5e4; % maximum number of value function iterations

mm.pi_tolerance  = 1e-8;    % convergence tolerance, profit function (WAS .001)
mm.T             = 50;      % horizon for calculating profit function
mm.tot_yrs       = 50;     % years to simulate, including burn-in (mm.burn)
mm.periods       = round(mm.tot_yrs*mm.pd_per_yr); % number of periods to simulate


mm.S         = 10000;    % number of potential exporting firms to simulate 
mm.burn      = 10;       %number of burn-in years
mm.max_match = 50;    % upper bound on number of matches to be counted for foreign market
mm.max_match_h = 70;    % Number of possible matches for domestic market

%% Cost function

mm.kappa1 = 2; %quadratic search costs 

mm.cost_h = @(x,net) (mm.cs_h * ((1+x).^mm.kappa1-(1+mm.kappa1*x))) /(mm.kappa1*(1 + log(net))^mm.gam);
mm.cost_f = @(x,net) (mm.cs_f * ((1+x).^mm.kappa1-(1+mm.kappa1*x))) /(mm.kappa1*(1 + log(net))^mm.gam);

mm.l_opt_func_h = @(a,net,pi,V_succ,V_fail,V_orig)... 
    max(max(((1+log(net))^mm.gam*(a*(pi+V_succ) + (1-a)*V_fail - V_orig)/mm.cs_h)+1,0).^(1/(mm.kappa1-1))-1,0);            
mm.l_opt_func_f = @(a,net,pi,V_succ,V_fail,V_orig)... 
    max(max(((1+log(net))^mm.gam*(a*(pi+V_succ) + (1-a)*V_fail - V_orig)/mm.cs_f)+1,0).^(1/(mm.kappa1-1))-1,0);            

%% Exogenous Jump Process Parameters
      
gam_h = 1 - 0.8412593; %DJ: calculations from COL_AR1.do, reversion coef, use Euler-Maruyama discretization of OU  
sig_h = 0.04688; %DJ: calculations from COL_AR1.do, root MSE, Euler-Maruyama discretization of OU 
gam_f = 1 - 0.948018; %DJ: from Jim's email 2023/1/11, AR1 reversion coef, Euler-Maruyama discretization of OU
sig_f = 0.11338; %DJ: from Jim's email 2023/1/11, AR1 root MSE, Euler-Maruyama discretization of OU

L_h = gam_h * mm.x_size; %lambda, arrival rate of shock
D_h = sig_h*L_h^(-.5); %delta, size of jump states
[Q_h,X_h] = makeq(L_h,D_h,mm.x_size);

L_f = gam_f * mm.x_size; %lambda, arrival rate of shock
D_f = sig_f*L_f^(-.5);   % delta, size of jump states
[Q_f,X_f] = makeq(L_f,D_f,mm.x_size);

%load exog_est/exog.mat  % old estimates for checking

% Shipment orders in home market are twice as frequent 
% (Alessandria, Kaboski, and Midrigan, AER, 2010)
mm.L_bF = X(6)/mm.pd_per_yr;
mm.L_bH = 1*mm.L_bF;

mm.max_shipsF = 3*round(mm.L_bF); % maximum within-period shipments is triple expected number
mm.poisCDF_shipmentsF   = poisscdf(1:1:mm.max_shipsF,mm.L_bF);

mm.max_shipsH = 3*round(mm.L_bH); % maximum within-period shipments is triple expected number
mm.poisCDF_shipmentsH   = poisscdf(1:1:mm.max_shipsH,mm.L_bH);

L_z = 4/mm.pd_per_yr; % four demand shock jumps per year (where is this from?)
D_z = X(5)/mm.pd_per_yr;
[Q_z,Z] = makeq(L_z,D_z,mm.z_size);
erg_pz = make_erg(L_z,D_z,Z); 

% Normal around zero (lognormal ultimately)
erg_pp = zeros(2 * mm.phi_size,1);
for k = 1:2 * mm.phi_size + 1
    erg_pp(k) = normpdf(-3 + 3/mm.phi_size * (k-1));
end
erg_pp = erg_pp./sum(erg_pp);
Phi = (-3:3/mm.phi_size:3)' * mm.sig_p;

%create Q_z with zeros on the diagonal (will use this later)
Q_z_d = Q_z;
Q_z_d(1:size(Q_z,1)+1:end) = 0;

%mm.actual_h     = actual_h; %actual indexes of home macro shocks for available years
%mm.actual_f     = actual_f; %actual indexes of foreign macro shocks for available years    

mm.Z            = Z;        %buyer productivities
mm.Phi          = Phi;      %seller productivies
mm.X_f          = X_f;      %foreign macro shocks
mm.X_h          = X_h;      %home macro shocks

mm.D_h          = D_h;      %size of jump in home macro shock
mm.Q_h          = Q_h;      %intensity matrix for home macro shock
mm.Q_h_d = Q_h;
mm.Q_h_d(1:size(Q_h,1)+1:end) = 0;

mm.D_f          = D_f;      %size of jump in foreign macro shock
mm.Q_f          = Q_f;      %intensity matrix for foreign macro shock
mm.Q_f_d = Q_f;
mm.Q_f_d(1:size(Q_f,1)+1:end) = 0;

mm.erg_pp       = erg_pp;   %ergodic distribution of seller productivities

mm.L_z          = L_z;      %arrival rate for jumps in other firms productivity
mm.D_z          = D_z;      %size of jump in other firms productivity
mm.Q_z          = Q_z;      %intensity matrix for demand shocks 
mm.Q_z_d        = Q_z_d;    %with zeros on the diagonal
mm.erg_pz       = erg_pz;   %ergodic distribution of demand shocks

mm.N_pt          = size(mm.Phi,1)*size(mm.theta2,2);
mm.pt_type = [kron((1:size(mm.Phi,1))',ones(size(mm.theta2,2),1)),kron(ones(size(mm.Phi,1),1),(1:size(mm.theta2,2))')];
mm.sim_firm_num_by_prod_succ_type = round(mm.erg_pp(mm.pt_type(:,1)).*mm.th2_pdf(mm.pt_type(:,2))'*mm.S);


%function TBA = generatePolicyAndValueFunctions(mm) 
 
    format long;

%    cf_num = 0;
%    type = 1; % initialize exporter type

 % Get policy functions
  
    policy = solvePolicyMain(mm);    

    lambda_h = policy.lambda_h;
    lambda_f = policy.lambda_f;
    c_val_h_orig = policy.c_val_h;
    c_val_f_orig = policy.c_val_f;
    value_f = policy.value_f;
    value_h = policy.value_h;
    % lambda_f (succ, trial, common succ rate (defunct), network size, prod of firm, F macro shock) 
    % lambda_h (common succ rate (defunct), known succ rate, network size, prod of firm, H macro shock)
    % c_val* (demand shock, prod of firm, macro shock)
        
    Q_size_f  = (mm.n_size + 1) *(mm.n_size+2)/2;    % was nn1 instead of nn2--both have same value, but need to check arguments of lambda_f
    Q_index_f = zeros(Q_size_f,3); % [index,trials,successes]
    
    counter = 1;
    for i=1:1:mm.n_size+1    % number of meetings, plus 1
        for ss=1:1:i % number of successes, plus 1
            Q_index_f(counter,:) = [counter,i,ss];
            counter = counter + 1;
        end
    end 

    Q_size_h = 70;
    Q_index_h = [(1:Q_size_h)',zeros(Q_size_h,1),(1:Q_size_h)']; 
    % (1) micro state (2) # matches (3) # successes 
    % structure chosen to be symmetric with Q_index_f. Only one column
    % really needed: we only keep track of successes in home market


%% calculate probabilities of generic firm types (theta draw by Z draw)
                
    th2_cdf = betacdf(mm.theta2,mm.af,mm.bf); % cdf for foreign theta draws
    th2_pdf = [th2_cdf(1),th2_cdf(2:mm.dim1)-th2_cdf(1:mm.dim1-1)]; 
         
    N_theta2 = size(mm.theta2',1);    % number of success rates
    N_Z      = size(mm.Z,1);          % number of buyer productivities
    N_Phi    = size(mm.Phi,1);        % number of exporter productivities

    % steady state macro shock dist
    mac_dist = normpdf(mm.X_f,0,0.0584);
    mac_dist = mac_dist./sum(mac_dist);   % (normalization patch for now)
    mac_dist_h = normpdf(mm.X_h,0,0.047);
    mac_dist_h = mac_dist_h./sum(mac_dist_h);   % (normalization patch for now)
     
    type_ndx = (1:1:size(mm.Phi,1)*size(mm.theta2,2)*size(mm.X_f,1))'; % # types is prod * succ prob * macro shocks
    
    pr_non_macro = kron(ones(size(mm.Phi,1),1),th2_pdf').*kron(mm.erg_pp,ones(size(mm.theta2,2),1)); % (prod, succ prob) joint pdf
    % NOTE: I changed erg_pz to erg_pp, i.e., probability of exporter (phi), not match shocks (z)
    
    pr_type = kron(mac_dist,ones(size(pr_non_macro,1),1)) .* repmat(pr_non_macro,size(mm.X_f,1),1);   % (prod, succ prob, macro shock) joint pdf
%    typemat = [type_ndx, kron((1:N_Xf)',ones(N_Phi * N_theta2,1)), repmat(kron(ones(N_Phi,1),(1:N_theta2)'),size(mm.Z,1),1), repmat(kron((1:N_Phi)',ones(N_theta2,1)), N_Phi,1)];
    typemat = [type_ndx, kron((1:size(mm.X_f,1))',ones(size(mm.Phi,1) * size(mm.theta2,2),1)),...
                         repmat(kron(ones(size(mm.Phi,1),1),(1:size(mm.theta2,2))'),size(mm.X_f,1),1),...
                         repmat(kron((1:size(mm.Phi,1))',ones(size(mm.theta2,2),1)), size(mm.X_f,1),1)];
    % 1: type indx, 2: macro shock indx, col 3: theta2 indx, 4: seller prod (phi) indx
        
    N_cli = 30; %maximum number of clients

%% construct intensity matrix for given firm type k

 [pmat_type_f,pmat_cum_f] = Q2pmat_v2(mm,mm.n_size+1,typemat,lambda_f,Q_size_f,Q_index_f,1);
 % last argument is cum. foreign market transition probability 
 

% not needed if using matchdat_gen_h2:
[pmat_type_h,pmat_cum_h,nn_h] = Q2pmat_v2(mm,mm.n_size+1,typemat,lambda_h,Q_size_h,Q_index_h,2);
% last output argument indicates maximum number of clients
  
%% construct transition probabilities for time interval Delta
% JT: division by sum() is needed only because of rounding error in expm()

pmat_msf = expm((1/mm.pd_per_yr).*mm.Q_f); % JT: added this (foreign macro state transitions)
% Q_f based on annual data, so need frac_of_year
pmat_msf = pmat_msf./sum(pmat_msf,2);
pmat_cum_msf = cumsum(pmat_msf')';   

pmat_msh = expm((1/mm.pd_per_yr).*mm.Q_h);
% Q_h based on annual data, so need frac_of_year
pmat_msh = pmat_msh./sum(pmat_msh,2); % JT: added this (home macro state transitions)
pmat_cum_msh = cumsum(pmat_msh')';   

% pmat_z = expm(frac_of_year.*mm.Q_z);  % JT: added this (buyer state transitions)
pmat_z = expm(mm.Q_z); 
pmat_z = pmat_z./sum(pmat_z,2);
pmat_cum_z = cumsum(pmat_z')';    

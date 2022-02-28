function policy = makeExporterTransitionProbabilities(mm,policy)

Q_size_f  = (mm.n_size + 1) *(mm.n_size+2)/2;
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
policy.th2_pdf = [th2_cdf(1),th2_cdf(2:mm.dim1)-th2_cdf(1:mm.dim1-1)];

% steady state macro shock dist
mac_dist = normpdf(mm.X_f,0,0.0584);
mac_dist = mac_dist./sum(mac_dist);   
mac_dist_h = normpdf(mm.X_h,0,0.047);
mac_dist_h = mac_dist_h./sum(mac_dist_h);   

policy.firm_type_prod_succ_macro = [(1:1:size(mm.Phi,1)*size(mm.theta2,2)*size(mm.X_f,1))', kron((1:size(mm.X_f,1))',ones(size(mm.Phi,1) * size(mm.theta2,2),1)),...
    repmat(kron(ones(size(mm.Phi,1),1),(1:size(mm.theta2,2))'),size(mm.X_f,1),1),...
    repmat(kron((1:size(mm.Phi,1))',ones(size(mm.theta2,2),1)), size(mm.X_f,1),1)];

%N_cli = 30; Not right...hope this isn't carried through

%% construct intensity matrix for given firm type k

[~,policy.pmat_cum_f] = Q2pmat_v2(mm,mm.n_size+1,policy.firm_type_prod_succ_macro,policy.lambda_f,Q_size_f,Q_index_f,1);

[~,policy.pmat_cum_h] = Q2pmat_v2(mm,mm.n_size+1,policy.firm_type_prod_succ_macro,policy.lambda_h,Q_size_h,Q_index_h,2);

%% construct transition probabilities for time interval Delta

pmat_msf = expm((1/mm.pd_per_yr).*mm.Q_f);
pmat_msf = pmat_msf./sum(pmat_msf,2);
policy.pmat_cum_msf = cumsum(pmat_msf')';

pmat_msh = expm((1/mm.pd_per_yr).*mm.Q_h);
pmat_msh = pmat_msh./sum(pmat_msh,2);
policy.pmat_cum_msh = cumsum(pmat_msh')';

pmat_z = expm(mm.Q_z);
pmat_z = pmat_z./sum(pmat_z,2);
policy.pmat_cum_z = cumsum(pmat_z')';
end
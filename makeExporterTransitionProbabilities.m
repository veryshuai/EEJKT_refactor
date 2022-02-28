function policy = makeExporterTransitionProbabilities(mm,policy)



% (1) micro state (2) # matches (3) # successes
% structure chosen to be symmetric with Q_index_f. Only one column
% really needed: we only keep track of successes in home market


%% calculate probabilities of generic firm types (theta draw by Z draw)

th2_cdf = betacdf(mm.theta2,mm.af,mm.bf); % cdf for foreign theta draws
policy.th2_pdf = [th2_cdf(1),th2_cdf(2:mm.dim1)-th2_cdf(1:mm.dim1-1)]; 

policy.firm_type_prod_succ_macro = [(1:1:size(mm.Phi,1)*size(mm.theta2,2)*size(mm.X_f,1))',...
    kron((1:size(mm.X_f,1))',ones(size(mm.Phi,1) * size(mm.theta2,2),1)),...
    repmat(kron(ones(size(mm.Phi,1),1),(1:size(mm.theta2,2))'),size(mm.X_f,1),1),...
    repmat(kron((1:size(mm.Phi,1))',ones(size(mm.theta2,2),1)), size(mm.X_f,1),1)];

%% construct intensity matrix for given firm type k

policy.pmat_cum_f = intensityToProbabilityForeign(mm,policy);

policy.pmat_cum_h = intensityToProbabilityHome(mm,policy);

%% construct transition probabilities for time interval Delta

policy = makeExogenousExogenousToFirmTransitionProbabilities(mm, policy);

end
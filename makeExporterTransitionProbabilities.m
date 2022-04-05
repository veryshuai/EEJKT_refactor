function policy = makeExporterTransitionProbabilities(mm,policy)

policy.firm_type_prod_succ_macro = [(1:1:size(mm.Phi,1)*size(mm.theta2,2)*size(mm.X_f,1))',...
    kron((1:size(mm.X_f,1))',ones(size(mm.Phi,1) * size(mm.theta2,2),1)),...
    repmat(kron(ones(size(mm.Phi,1),1),(1:size(mm.theta2,2))'),size(mm.X_f,1),1),...
    repmat(kron((1:size(mm.Phi,1))',ones(size(mm.theta2,2),1)), size(mm.X_f,1),1)];

%% construct intensity matrix for given firm type k

policy = intensityToProbabilityForeign(mm,policy);

policy = intensityToProbabilityHome(mm,policy);

%% construct transition probabilities for time interval Delta

policy = makeExogenousFirmTransitionProbabilities(mm, policy);

end
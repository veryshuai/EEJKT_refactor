function policy = makeExogenousFirmTransitionProbabilities(mm, policy)

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
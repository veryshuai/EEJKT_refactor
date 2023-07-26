function policy = makeExporterTransitionProbabilitiesExchangeRateShock(mm,policy)

%% construct intensity matrix for given firm type k

policy = intensityToProbabilityForeign_exch_rate_shk(mm,policy);

end

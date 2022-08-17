function D = calculateDistanceAndPrint(simMoms,mm,X)

[Data, W, Model] = read_in_and_organize_moments(simMoms);

error = Data'-Model;

[Data, error, Model] = do_not_match_some_moments(Data, error, Model);

D = log((error'/W)*error);

inv_W = W^-1;
err_comp = @(sta,fin) error(sta:fin)' * inv_W(sta:fin,sta:fin) * error(sta:fin); 

real_moms_and_sim_moms = cat(2,Data',Model);
print_diagnostics_to_standard_output(D, X, real_moms_and_sim_moms, err_comp, simMoms,mm);
%plots; 
%summary_tables_v2; 
%estimate_summary;

shouldMatchMoments(real_moms_and_sim_moms,D,"overwrite","results/shouldMatchMomentsData");
% shouldMatchMoments(real_moms_and_sim_moms,D,"test","results/shouldMatchMomentsData");
 

end

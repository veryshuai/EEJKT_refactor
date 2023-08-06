function [D,D_alt,real_moms_and_sim_moms] = calculateDistanceAndPrint(simMoms,mm,X)

[Data, Data_alt, W, W_alt, Model, Model_alt] = read_in_and_organize_moments(simMoms);

%% Statistic using coefficients from degree distribution regression
error = Data'-Model;
[Data, error, Model] = do_not_match_some_moments(Data, error, Model);
D = log((error'/W)*error);

inv_W = W^-1;
err_comp = @(sta,fin) error(sta:fin)' * inv_W(sta:fin,sta:fin) * error(sta:fin); 
real_moms_and_sim_moms = cat(2,Data',Model);
print_diagnostics_to_standard_output(D, X, real_moms_and_sim_moms, err_comp, simMoms,mm);

% plots_v2(simMoms); 
 summary_tables(simMoms,mm); 

%% Statistics using shares distribution of buyers per seller 
error_alt = Data_alt' - Model_alt;
D_alt     = log((error_alt'/W_alt)*error_alt);
[Data_alt, error_alt, Model_alt] = do_not_match_some_moments_alt(Data_alt, error_alt, Model_alt);

inv_W_alt = W_alt^-1;
err_comp_alt = @(sta,fin) error_alt(sta:fin)' * inv_W_alt(sta:fin,sta:fin) * error_alt(sta:fin); 

real_moms_and_sim_moms_alt = cat(2,Data_alt',Model_alt);
print_diagnostics_to_standard_output_alt(D_alt, D, X, real_moms_and_sim_moms_alt, err_comp_alt, simMoms,mm);

%%
% shouldMatchMoments(real_moms_and_sim_moms,D,"overwrite","results/shouldMatchMomentsData");
% shouldMatchMoments(real_moms_and_sim_moms,D,"test","results/shouldMatchMomentsData");

end

function [D,W,error] = distance(X)

rng(80085,'twister');
seed_crand(80085);

format long;

X = [-3.60529  -3.87941  0.23156  2.46487  1.88176  15.42634 0.38246  11.72211  1.38618  -1.21819  13.00238  -6.13506];
    
mm = setModelParameters(X);
policy = generatePolicyAndValueFunctions(mm);
simulateMomentsMain;

[Data, W, Model] = read_in_and_organize_moments(simMoms);

error = Data'-Model;

[Data, error, Model] = do_not_match_some_moments(Data, error, Model);

D = log((error'/W)*error);

dat_mod_moms = cat(2,Data',Model);  
inv_W = W^-1;
err_comp = @(sta,fin) error(sta:fin)' * inv_W(sta:fin,sta:fin) * error(sta:fin); 

print_diagnostics_to_standard_output(D, X, dat_mod_moms, err_comp, simMoms,mm);
plots; 
summary_tables_v2; 
estimate_summary;

shouldMatchMoments(dat_mod_moms,D,"test","results/shouldMatchMomentsData");
   
end

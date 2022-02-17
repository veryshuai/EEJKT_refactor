function [D,W,error] = distance(X)

rng(80085,'twister');
seed_crand(80085);

format long;

X = [-3.60529  -3.87941  0.23156  2.46487  1.88176  15.42634 0.38246  11.72211  1.38618  -1.21819  13.00238  -6.13506];
    
tic;
X2params;
SetParams;
inten_sim_v1;
time1 = toc;

tic2 = tic;
discrete_sim_parfor3;
time2 = toc(tic2);

fprintf('\r\n inten_sim run time: %.2f\n', time1); 
fprintf(' discrete_sim run time: %.2f\n', time2);

[Data, W, Model] = read_in_and_organize_moments(simMoms);

error = Data'-Model;

[Data, error, Model] = do_not_match_some_moments(Data, error, Model);

W_D   = log((error'/W)*error);
Old_D = norm(error)/norm(Data');

D = W_D;

%%  data/model comparison  
mmm = cat(2,Data',Model);  
inv_W = W^-1;
err_comp = @(sta,fin) error(sta:fin)' * inv_W(sta:fin,sta:fin) * error(sta:fin); 

print_diagnostics_to_standard_output(W_D, Old_D, X, mmm, err_comp, simMoms,mm);
plots; 
summary_tables_v2; 
estimate_summary;

shouldMatchMoments(mmm,W_D,"test");
   
end

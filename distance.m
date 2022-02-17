function [D,W,error] = distance(X)

rng(80085,'twister');
seed_crand(80085);

format long;

X = [-3.60529  -3.87941  0.23156  2.46487  1.88176  15.42634 0.38246  11.72211  1.38618  -1.21819  13.00238  -6.13506];

try
    
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


%% Targets and weights

[Data, W] = target_stats();

%% Realizations

    match_death_coefsSIM = [match_exit_rate;beta_match_exit(2:5)]; % [match exit rate, 1st yr. dummy, lnXf(ijt), ln(match age), ln(exporter age),mse]
    match_ar1_coefsSIM   = [ybar_match;beta_match(2:4);mse_match_ar1]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]   
    loglog_coefsSIM      = [b_degree]; % [intercept, slope, quadratic term]
    mavshipSIM           = [avg_ln_ships]; % average ln(# shipments) 
    exp_dom_coefsSIM     = [ybar_hfsales;beta_hfsales(2);mse_hf]; % [mean dep var.,coef,MSE]  
    dom_ar1_coefsSIM     = [ybar_fsales_h;beta_fsales_h(2);mse_h]; % [mean dep var.,coef,MSE] 
    ln_haz_coefsSIM      = [mean_ln_haz;b_haz(2:6)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)] 
    last_match_coefsSIM  = [mkt_exit_rate;beta_mkt_exit(2:6)]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)]            
    succ_rate_coefsSIM   = [mean_succ_rate;b_succ_rate(2)]; % [mean succ rate, ln(1+meetings)]
    sr_var_coefsSIM      = [mean_usq_succ;b_usq_succ(2)]; % [mean dep. var, ln(1+meetings)]
    for_sales_shrSIM     = [avg_expt_rate]; % mean share of exports to U.S. in total sales 
    exp_fracSIM          = [share_exptr]; % fraction of firms exporting to U.S.  

 Model = cat(1,match_death_coefsSIM,match_ar1_coefsSIM,loglog_coefsSIM,...
    mavshipSIM,exp_dom_coefsSIM,dom_ar1_coefsSIM,ln_haz_coefsSIM,...   
    last_match_coefsSIM,succ_rate_coefsSIM,sr_var_coefsSIM,for_sales_shrSIM,...    
    exp_fracSIM);   


error = Data'-Model;

[Data, error, Model] = do_not_match_some_moments(Data, error, Model);

W_D   = log((error'/W)*error);
Old_D = norm(error)/norm(Data');

D = W_D;

nanflag = isnan(D); 
if nanflag>0
        D = D * 10; 
end

%%  data/model comparison  
mmm = cat(2,Data',Model);  

    % Create Diagnostics
    fprintf('\r\n weighted metric:   %.15f\n', W_D); 
        
    %Simple unweighted loss

    fprintf(' unweighted metric: %.15f\n', Old_D); 
    
   fprintf('\r\n params = ');
   fprintf('\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',X(1:6));
   fprintf('\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',X(7:12));
   fprintf( '\r\n  ');   
    
    format shortG
    fprintf('\r\n moments: ');
    cat(2,mmm(1:10,:),mmm(11:20,:),mmm(21:30,:),[mmm(31:38,:);zeros(2,2)])
    format long
    

    fprintf('\r\n Fit diagnostics: ');  
    inv_W = W^-1;
    err_comp = @(sta,fin) error(sta:fin)' * inv_W(sta:fin,sta:fin) * error(sta:fin);    
    fprintf('\r\n match_death_coefs  = %.3f\n',err_comp(1,5));
    fprintf(' match_ar1_coefs    = %.3f\n',err_comp(6,10));
    fprintf(' log_log_coefs      = %.3f\n',err_comp(11,13));
    fprintf(' av_shipments       = %.3f\n',err_comp(14,14));
    fprintf(' exp_dom            = %.3f\n',err_comp(15,17));
    fprintf(' dom_ar1            = %.3f\n',err_comp(18,20));
    fprintf(' match_lag_coef     = %.3f\n',err_comp(21,26));
    fprintf(' last_match_coef    = %.3f\n',err_comp(27,32));   
    fprintf(' succ_rate_coef     = %.3f\n',err_comp(33,34));
    fprintf(' sr_var_coef        = %.3f\n',err_comp(35,36));
    fprintf(' for_sales_shr_coef = %.3f\n',err_comp(37,37));
    fprintf(' exp_frac_coef      = %.3f\n',err_comp(38,38));
    
    fprintf('\r\n number of exporters per yr = %.3f\n',agg_nexptr/(mm.tot_yrs - mm.burn));
    fprintf(' maximum number of clients  = %.3f\n',size(ff_sim_max,1));
    fprintf(' number of firms per yr     = %.3f\n',agg_nfirm/(mm.tot_yrs - mm.burn));
    fprintf( '\r\n  '); 
    
%%   write results to text files
      fitvec = [W_D Old_D];
      fileID2 = fopen('results/ga_fitlog.txt','a');
      fprintf(fileID2,'\r\n fit metric = ');
      dlmwrite('results/ga_fitlog.txt',fitvec,'-append','precision',12);
      fclose(fileID2);
 

      fileID1 = fopen('results/ga_running_output.txt','a');
      fprintf(fileID1,'\r\n fit metrics (weighted and unweighted): ');
      dlmwrite('results/ga_running_output.txt',fitvec,'-append','precision',12);
    
      fprintf(fileID1, '\r\n parameters: ');
      fprintf(fileID1, '\r\n%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f',X(1:6));
      fprintf(fileID1, '\r\n%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f',X(7:12));
      fprintf(fileID1, '\r\n  ');
  
      fprintf(fileID1, '\r\n moments: ');   
      fprintf(fileID1, '\r\n%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f',full(mmm(1:10,:)));  
      fprintf(fileID1, '\r\n  ');
      fprintf(fileID1, '\r\n%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f',full(mmm(11:20,:)));  
      fprintf(fileID1, '\r\n  ');
      fprintf(fileID1, '\r\n%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f',full(mmm(21:30,:)));  
      fprintf(fileID1, '\r\n  ');
      fprintf(fileID1, '\r\n%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f',full(mmm(31:38,:)));  
      fprintf(fileID1, '\r\n  ');               
      fclose(fileID1);

 
plots
summary_tables_v2 
estimate_summary

shouldMatchMoments(mmm,W_D,"test");

catch err
    
        D = 1e12;
        W = 1; 
        error = 1;
end 
   
end

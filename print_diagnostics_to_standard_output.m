function print_diagnostics_to_standard_output(W_D, Old_D, X, mmm, err_comp,simMoms,mm)
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
    
    fprintf('\r\n number of exporters per yr = %.3f\n',simMoms.agg_nexptr/(mm.tot_yrs - mm.burn));
    fprintf(' maximum number of clients  = %.3f\n',size(simMoms.ff_sim_max,1));
    fprintf(' number of firms per yr     = %.3f\n',simMoms.agg_nfirm/(mm.tot_yrs - mm.burn));
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
end
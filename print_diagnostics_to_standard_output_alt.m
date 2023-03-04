function print_diagnostics_to_standard_output_alt(D, D0, X, mmm, err_comp,simMoms,mm)

fprintf('\r\n alternative fit metric: %.15f\n', D); 
    
   fprintf('\r\n params = ');
   fprintf('\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',X(1:6));
   fprintf('\r%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f',X(7:end));
   fprintf( '\r\n  ');   
    
    format shortG
    fprintf('\r\n alternative moments: ');
    cat(2,mmm(1:10,:),mmm(11:20,:),mmm(21:30,:),[mmm(31:36,:);zeros(4,2)])
    format long
    

    fprintf('\r\n Fit diagnostics: ');  
   
    fprintf('\r\n match_death_coefs  = %.3f\n',err_comp(1,5));
    fprintf(' match_ar1_coefs    = %.3f\n',err_comp(6,10));   
    fprintf(' av_shipments       = %.3f\n',err_comp(11,11));
    fprintf(' exp_dom            = %.3f\n',err_comp(12,14));
    fprintf(' dom_ar1            = %.3f\n',err_comp(15,17));
    fprintf(' match_haz_coef     = %.3f\n',err_comp(18,23));   
    fprintf(' succ_rate_coef     = %.3f\n',err_comp(24,25));
    fprintf(' sr_var_coef        = %.3f\n',err_comp(26,27));
    fprintf(' for_sales_shr_coef = %.3f\n',err_comp(28,28));
    fprintf(' exp_frac_coef      = %.3f\n',err_comp(29,29));
    fprintf(' BPS distribution   = %.3f\n',err_comp(30,36));
    
    fprintf('\r\n number of exporters per yr = %.3f\n',simMoms.agg_nexptr/(mm.tot_yrs - mm.burn));
    fprintf(' maximum number of clients  = %.3f\n',max(simMoms.max_countD));
    fprintf(' number of firms per yr     = %.3f\n',simMoms.agg_nfirm/(mm.tot_yrs - mm.burn));
    fprintf( '\r\n  '); 


      fileID2 = fopen('results/ga_fitlog_2Arv_1Scl_opt.txt','a');
      fprintf(fileID2,'\r\n fit metrics, original and alternative:  ');
      dlmwrite('results/ga_fitlog_2Arv_1Scl_opt.txt',[D0 D],'-append','precision',12);
      fclose(fileID2);
    
%       fileID2 = fopen('results/ga_running_output_2Arv_2Scl.txt','a');
%       fprintf(fileID2,'\r\n alternative fit metric = ');
%       dlmwrite('results/ga_running_output_2Arv_2Scl.txt',D,'-append','precision',12);
%       fclose(fileID2);
 

      fileID1 = fopen('results/ga_running_output_2Arv_1Scl_opt.txt','a');
      fprintf(fileID1,'\r\n alternative fit metric: ');
      dlmwrite('results/ga_running_output_2Arv_1Scl_opt.txt',D,'-append','precision',12);
    
      fprintf(fileID1, '\r\n parameters: ');
      fprintf(fileID1, '\r\n%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f',X(1:6));
      fprintf(fileID1, '\r\n%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f',X(7:end));
      fprintf(fileID1, '\r\n  ');
  
      fprintf(fileID1, '\r\n alternative moments: ');   
      fprintf(fileID1, '\r\n%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f',full(mmm(1:10,:)));  
      fprintf(fileID1, '\r\n  ');
      fprintf(fileID1, '\r\n%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f',full(mmm(11:20,:)));  
      fprintf(fileID1, '\r\n  ');
      fprintf(fileID1, '\r\n%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f',full(mmm(21:30,:)));  
      fprintf(fileID1, '\r\n  ');
      fprintf(fileID1, '\r\n%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f',full([mmm(31:36,:);zeros(4,2)])); 
      fprintf(fileID1, '\r\n  ');               
      fclose(fileID1);
end
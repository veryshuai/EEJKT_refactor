function [Data, W] = target_stats()
% returns the data statistics used in the loss function

    %% TARGET SOURCES: 
    % match_death_coefs from CJ 5/12/15 p. 11, also in CJ spreadsheet of 3/29/16
    % match_ar1_coefs   from CJ 5/12/15 p. 11, also in CJ spreadsheet of 3/29/16
    % loglog_coefs      from CJ 5/12/15 p. 11, also in CJ spreadsheet of 3/29/16
    % mavship           from CJ 6/11/15
    % ex_dom_coefs      from Marcela 6/22/17 except for mean of dep. var, which is in $ and from C.J.
    % dom_ar1_coefs     from Marcela 6/22/17 except for mean of dep. var*
    % ln_haz_coefs      from US customs records, disclosed 12-12-16 & 05-07-19
    % last_match_coefs  from US customs records, disclosed 12-12-16 & 05-07-19 
    % succ_rate_coefs   from RDC disclosure 2-3-17 & 05-07-19
    % sr_var_coefs      from RDC disclosure 2-3-17 & 05-07-19
    % for_sales_shr     from Marcela 6/22/17
    % exp_frac          from Marcela 6/22/17
    
    % *avg. ln(exports|#matches>0) = 
    %  avg. ln(match sales) + avg. ln(#matches) +var(match sales)/2 
    %  = 10.6653 + 0.37463 + 1.62678 = 12.66671
    % *In constant pesos, difference between avg. ln(domsal|domsal>0) and   
    %  avg. ln(exports|export>0,domsal>0) is 14.46852 -12.84174 = 1.62678.
    %  Add to avg. ln(exports)in $ from US Census to get $-equivalent 
    %  avg. ln(domsal):  12.66671 + 1.62678 = 14.29349
    
 % coefficient vector with means of dep. var. instead of intercepts, revised 2-24-17 
 
    match_death_coefsDAT = [0.3951 0.03421 -0.03163 -0.05370 -0.02778]; % [mean D, R(ijt), lnXf(ijt), ln(match age), ln(exporter age)]
    match_ar1_coefsDAT   = [10.6653 0.82637 0.32834 0.06312 1.20786^2]; % [mean ln Xf(ijt), ln Xf(ijt-1), R(ijt-1), ln(exporter age)]   
    loglog_coefsDAT      = [-5.9731 -1.88130 -0.05446];  % [intercept, slope, quadratic term]
    mavshipDAT           = 0.9706402010; % average ln(# shipments) 
%   exp_dom_coefsDAT     = [11.0399,0.3228142,2.606^2]; % [mean dep var.,coef,MSE]  
    exp_dom_coefsDAT     = [11.0399,0.3726,2.4679^2]; % [mean dep var.,coef,MSE]  REVISED EXCEPT MEAN 9/6/22
%   dom_ar1_coefsDAT     = [14.29349,0.9764422,0.46207^2]; % [mean dep var.,coef,MSE] 
    dom_ar1_coefsDAT     = [14.29349,0.97689,0.45104^2]; % [mean dep var.,coef,MSE] REVISED EXCEPT MEAN 9/6/22
    ln_haz_coefsDAT      = [-0.7188,-0.8181,0.3117,-1.1323,2.4514,-0.7082]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)] 
    last_match_coefsDAT  = [0.635,-0.6290,0.1205,0.0211,0.5976,-0.1290]; % [mean dep. var, ln(1+a), ln(1+a)^2, ln(1+r), ln(1+r)^2, ln(1+a)*ln(1+r)]            
    succ_rate_coefsDAT   = [0.413,0.093];  % [mean succ rate, ln(1+meetings)]
    sr_var_coefsDAT      = [0.0912,-0.060]; % [mean succ rate, ln(1+meetings)]
%   for_sales_shrDAT     =  0.1270; % mean share of exports to U.S. in total sales
    for_sales_shrDAT     =  0.1291; % mean share of exports to U.S. in total sales REVISED 9/6/22
%   exp_fracDAT          =  0.1023; % fraction of firms that export to U.S. 
    exp_fracDAT          =  0.0954; % fraction of firms that export to U.S. REVISED 9/6/22
    
   Data = [match_death_coefsDAT,match_ar1_coefsDAT,loglog_coefsDAT,mavshipDAT,exp_dom_coefsDAT,...
       dom_ar1_coefsDAT,ln_haz_coefsDAT,last_match_coefsDAT,succ_rate_coefsDAT,sr_var_coefsDAT,...
       for_sales_shrDAT,exp_fracDAT];
     
%      Data = [match_death_coefsDAT (1-5),match_ar1_coefsDAT (6-10),loglog_coefsDAT (11-13),mavshipDAT (14),exp_dom_coefsDAT (15-17),...
%       dom_ar1_coefsDAT (18-20),match_lag_coefsDAT (21-26),last_match_coefsDAT (27-32),succ_rate_coefsDAT (33-34),sr_var_coefsDAT (35-36)];
  
    %% covariance matrices--covariances with dep. var. mean set to zero
 
match_death_coefsCOV = ...    
    [0.000427699832   0.000000000000   0.000000000000   0.000000000000  0.000000000000; 
     0.000000000000   0.000137342677   0.000001753036   0.000081273853 -0.000005694758; 
     0.000000000000   0.000001753036   0.000002501209  -0.000000672050 -0.000000605574;
     0.000000000000   0.000081273853  -0.000000672050   0.000081302304 -0.000027033495;
     0.000000000000  -0.000005694758  -0.000000605574  -0.000027033495  0.000042624529];
 
match_ar1_coefsCOV = ...
     [0.002362705419  0.000000000000  0.000000000000  0.000000000000  0       ;
      0.000000000000  0.000014721579  0.000014300618 -0.000004772804  0       ;
      0.000000000000  0.000014300618  0.000332671396  0.000124115947  0       ;
      0.000000000000 -0.000004772804  0.000124115947  0.000193084788  0       ;
      0               0               0               0               0.00012];

loglog_coefsCOV = ...
  [0.022689638548   0.000000000000   0.000000000000;
   0.000000000000   0.012614810554  -0.002268716169;             
   0.000000000000  -0.002268716169   0.000443232645].*0.01;  % giving extra weight           

    mavshipCOV = 0.00415553^2;  % note that I squared the standard error here

 % The following two matrices come from Marcela's printout (DIAN only 21 Feb. 2013)

%     exp_dom_coefsCOV = ...
%       [0.0133969  0             0       ;     % first element is std. dev. of mean of dep. var
%        0           0.0120728    0       ;     % 2.537503/sqrt(35877) = 0.0133969
%        0           0            0.089016].^2; % see notes on printout for se(RSME)
% 
%  dom_ar1_coefsCOV = ...
%        [0.0050684   0           0        ;    % first element is std. dev. of mean of dep. var
%         0           0.0008622   0        ;    % 1.718554/sqrt(114968) = 0.0050684
%         0           0           0.000958].^2; % see notes on printouts for se(RSME)
    


 % The following two matrices WILL come from Marcela's revised printout (Outputs_EAM_moms_2022) 

       exp_dom_coefsCOV = ...
      [0.0133969  0             0       ;     % first element is std. dev. of mean of dep. var
       0           0.0120728    0       ;     % 2.537503/sqrt(35877) = 0.0133969
       0           0            0.089016].^2; % see notes on printout for se(RSME)
   

    dom_ar1_coefsCOV = ...
       [0.0050684   0           0        ;    % first element is std. dev. of mean of dep. var
        0           0.0008622   0        ;    % 1.718554/sqrt(114968) = 0.0050684
        0           0           0.000958].^2; % see notes on printouts for se(RSME)
    
    
% The following matrices are based on U.S. customs records. 

  match_lag_coefsCOV = ...
   [0.00621 0      0      0      0      0        ;
    0       0.1128 0      0      0      0        ;
    0       0      0.0168 0      0      0        ;
    0       0      0      0.2962 0      0        ;
    0       0      0      0      0.3956 0        ;
    0       0      0      0      0      0.1341]^2; %Results disclosed 12-12-16 & 5-7-19


  last_match_coefsCOV = ...
   [0.00253 0      0      0      0      0        ;
    0      0.0422 0      0      0      0        ;
     0      0      0.0063 0      0      0        ;
     0      0      0      0.1101 0      0        ;
     0      0      0      0      0.1470 0        ;
     0      0      0      0      0      0.0501]^2; % Results disclosed 12-12-16 & 5-07-19

succ_rate_coefsCOV = ...  % based on RDC disclosure 2-3-17
    [0.00153^2   0.00000000 ;
     0.00000000  0.00000683];
 
 sr_var_coefsCOV = ...    % based on RDC disclosure 2-3-17 & 5-07-19
     [0.000265^2  0.00000000;
      0.00000000  0.00000012];
  
 for_sales_shrCOV = (.2409052^2)/10838;   % from Marcela's eam_moms_out.pdf  
 exp_fracCOV      = (.2937100^2)/113656;  % from Marcela's eam_moms_out.pdf
%%

    
 % weighting matrix for moment vector revised 12-17-16  
   
        W = blkdiag(match_death_coefsCOV,...
       match_ar1_coefsCOV,loglog_coefsCOV,mavshipCOV,exp_dom_coefsCOV,...
       dom_ar1_coefsCOV,match_lag_coefsCOV,last_match_coefsCOV,...
       succ_rate_coefsCOV,sr_var_coefsCOV,for_sales_shrCOV,exp_fracCOV);
   
end

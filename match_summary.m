function [FirmCount,exit_by_age,brooks] =  match_summary_v3(simMoms,mm)

% This function generates failure rates by initial sales and match age.
% More precisely, by type group and match age, where types are classified
% according to their initial sales.

%%  Preliminary data preparation
     succ_matches = simMoms.agg_mat_yr_sales;
     dud_matches  = simMoms.agg_dud_matches;
     match_recs = [succ_matches;dud_matches]; % stack successful and dud matches       
%    match_recs = succ_matches; % exclude duds as a test
      
  %  match_recs: [year, type, firm_ID, sales, shipments, boy Z, eoy Z, match age, firm age]       

     ff_ship         = match_recs(:,5) > 0;    % positive shipmments
     match_recs      = match_recs(ff_ship,:);  % select matches with shipments>0
     match_recs(:,1) = floor(match_recs(:,1)/mm.pd_per_yr);   % restate time in years
     match_recs(:,8) = floor(match_recs(:,8)/(mm.pd_per_yr+1e-6))+1; % restate match age in years
     match_recs(:,9) = floor(match_recs(:,9)/mm.pd_per_yr+1e-6)+1; % restate firm age in years
     max_age         = max(match_recs(:,8));
 
     all_matches = match_recs;        
     match_TF    = all_matches(:,2) + 0.0001*floor(all_matches(:,3)); % firm type by firmID identifier
     TF_list     = unique(match_TF);
     NumTF       = size(TF_list,1);    
     TF_matdat   = cell(NumTF,1);
     firstYrCount = zeros(NumTF,1);
     
     % pick off new matches for year-type-firm jj 
    for firmID=1:NumTF
      ff = match_TF == ones(size(all_matches,1),1).*TF_list(firmID);
      TF_matdat{firmID} = all_matches(ff,:);
      firstYrCount(firmID,1) = sum(TF_matdat{firmID}(:,4)==1);
    end
   
    agg_mat_lifecycle = double.empty(0,5*max_age+1);
    
%    for firmID = 4
      for firmID=1:NumTF
      fprintf('\r firm ID = %0.0f\n',firmID);  
      all_age = TF_matdat{firmID};
      new_mat = all_age(:,8)<=1;
      max_age_jj = max(all_age(:,8));
      NN = sum(new_mat);
      mat_lifecycle = zeros(NN,5*max_age);
      mat_lifecycle(:,1:5) = [all_age(new_mat,1),all_age(new_mat,8),...
               all_age(new_mat,6:7),all_age(new_mat,4)];         
    % all_age: 1) year, 2) type, 3) firm_ID, 4) sales, 5) shipments, 6) boy Z, 7) eoy Z, 8) match age, 9) firm age] 
    % mat_lifecycle(:,1:4): [year, match age, boy Z, eoy Z, sales]
     
      for aa = 2:max_age_jj
          % grab the aa-year-old matches  
           aa_yr_old = all_age(:,8) == aa;
           aa_cohort = all_age(aa_yr_old,:);

     exitloop = 1;      
          while exitloop > 0 && size(aa_cohort,1) > 0
              
           if firstYrCount(firmID,1)==0 
           fprintf('\r firm ID = %0.0f, number of 1st year matches = %0.0f\n',[firmID,firstYrCount(firmID,1)]);
           exitloop = 0;               
           end
              
              ii =1; 
              while ii <= NN
              fprintf('\r firm ID = %0.0f, match age = %0.0f, match number = %0.0f\n',[firmID,aa,ii]);
              lag_yr  = mat_lifecycle(ii,5*(aa-2)+1);
              lag_age = mat_lifecycle(ii,5*(aa-2)+2);
              lag_eopZ = mat_lifecycle(ii,5*(aa-2)+4);

              % find a current year match to splice with last year's match ii
                flg = (lag_eopZ>0).*(aa_cohort(:,6)==lag_eopZ).*...
                      (aa_cohort(:,8)==lag_age + 1).*(aa_cohort(:,1)==lag_yr + 1);
                  

              if sum(flg)>0
                  mat_cont = find(flg==1,1);
                  lb = 5*(aa-1)+1;
                  ub = 5*(aa-1)+5;
                  mat_lifecycle(ii,lb:ub) = [aa_cohort(mat_cont,1),aa_cohort(mat_cont,8),...
                       aa_cohort(mat_cont,6:7),aa_cohort(mat_cont,4)]; 
              
                  % remove matched record from matchable subset 
                  keepers = ones(size(aa_cohort,1),1);
                  keepers(mat_cont,1) = 0;
                  exitloop = sum(keepers)>0;
                  aa_cohort = aa_cohort(logical(keepers),:); 
              else
                  exitloop = 0;
              end
 %          end
              ii = ii + 1;
              end
 %            aa
         end
      end
       agg_mat_lifecycle = [agg_mat_lifecycle; [firmID*ones(size(mat_lifecycle,1),1), mat_lifecycle]];
    end       
  
%% Construct match survival table

% Find sales quartiles for new matches
     yr1 = all_matches(:,8)<=1; % pick matches in first year   
     new_matches = all_matches(yr1,:); % matches in their 1st yr.    
     new_sales   = sort(new_matches(:,4),1); % sort year-type-firms by mean sales of new matches
     Nrows       = size(new_sales,1);
%    cum_cnt     = cumsum(new_sales./sum(new_sales)); % cum distrib. of # matches, sorted by avg. sales
     cum_cnt = cumsum(ones(Nrows,1)./Nrows);    
     ndx_q = cell(4,1); % set of dummy variables for match count quartiles
     mean_q = zeros(4,1);
     match_cnt_q = zeros(4,1);
     size_cutoff = zeros(4,1);
     upb = 0.25;
     lowb = 0;
     for ii = 1:4
        ndx_q{ii}  = cum_cnt>lowb & cum_cnt<=upb; % vector of row identifiers for types in sales quartile ii
        lowb = upb;
        upb = upb + 0.25;
        % quartile-wide avg. sales per match:
        mean_q(ii) =mean(new_sales(ndx_q{ii},1));
        match_cnt_q(ii) = size(new_sales(ndx_q{ii},1),1);
        size_cutoff(ii) = max(new_sales(ndx_q{ii}));  
     end
     
% label matches by initial size quartiles   
     initial_size = 1 ...
     + 1*((agg_mat_lifecycle(:,6) > size_cutoff(1)).*(agg_mat_lifecycle(:,6) <= size_cutoff(2)))...    
     + 2*((agg_mat_lifecycle(:,6) > size_cutoff(2)).*(agg_mat_lifecycle(:,6) <= size_cutoff(3)))...
     + 3*(agg_mat_lifecycle(:,6) > size_cutoff(3));
 
 
% match survival and age 
   maxAge = 10;   
   survive = zeros(maxAge,4);
   meanSale = zeros(maxAge,4);
   matchCount  = zeros(maxAge,4);
      for qq=1:4
          group = initial_size==qq;
          matchBySize = agg_mat_lifecycle(group,:);

          for aa = 1:maxAge
              currSale = (aa-1)*5 + 5;
              nextYrSale =  aa*5 + 5;
              meanSale(aa,qq)   = mean(matchBySize(matchBySize(:,1+currSale)>0,1+currSale));
              matchCount(aa,qq) = sum(matchBySize(:,1+currSale)>0);
              if aa > 1
              survive(aa,qq) = matchCount(aa,qq)/matchCount(aa-1,qq);
              end
          end
          
      end
 
      exit_by_age = 1 - survive;
      
%%  Construct brooks table
BrooksYrs =10;
maxAge   = 50;
firmMatchCount = zeros(maxAge,NumTF); 
firmCount      = zeros(maxAge,NumTF); 
firmTotSales   = zeros(maxAge,NumTF); 

   for firmID = 1:NumTF 
   firmMatches = TF_matdat{firmID};  
      for aa=1:maxAge
      firmCount(aa,firmID)      = size(firmMatches(firmMatches(:,9)==aa,:),1)>0;
      firmMatchCount(aa,firmID) = sum(firmMatches(firmMatches(:,9)==aa,4)>0);
      firmTotSales(aa,firmID)   = sum(firmMatches(firmMatches(:,9)==aa,4));
      end    
   end
   
MatchCount = sum(firmMatchCount,2);
TotFirms   = sum(firmCount,2);
TotSales   = sum(firmTotSales,2); 
AvgSales   = TotSales./MatchCount;

brooks = [TotFirms,TotSales,AvgSales];
brooks = brooks(1:BrooksYrs,:);
%% Construct degree distribution (This version not used)

maxMatches=50;
FirmFreq   = zeros(maxMatches,maxAge);
maxMatches = max(max(firmMatchCount));
    for Nmatch = 1:maxMatches
       for aa = 1:maxAge
       ff = find(firmMatchCount(aa,:) == Nmatch);
       FirmFreq(Nmatch,aa) = sum(firmMatchCount(aa,ff)>0);
    end
    end         
 FirmCount = sum(FirmFreq,2);
 %FirmFreq  = FirmCount./FirmCount(1);


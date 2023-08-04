function [DegreeDistCount,exit_by_age,brooks] =  match_summary(simMoms,mm)

% This function generates failure rates by initial sales and match age.
% It also constructs a summary (Brooks) table for the cohort maturation process.

%%  Preliminary data preparation
     succ_matches = simMoms.agg_mat_yr_sales;
     dud_matches  = simMoms.agg_dud_matches;
     match_recs   = [succ_matches;dud_matches]; % stack successful and dud matches            
  %  match_recs: [period, type, firm_ID, sales, shipments, boy Z, eoy Z, match age, firm age]       

     ff_ship         = match_recs(:,5) > 0;    % positive shipmments
     match_recs      = match_recs(ff_ship,:);  % select matches with shipments>0
     match_recs(:,1) = floor(match_recs(:,1)/(mm.pd_per_yr+1e-6)) + 1;          % restate time in years
     match_recs(:,8) = floor(match_recs(:,8)/(mm.pd_per_yr+1e-6))+1; % restate match age in years
     match_recs(:,9) = floor(match_recs(:,9)/(mm.pd_per_yr+1e-6))+1;   % restate firm age in years
     max_age         = max(match_recs(:,8));
 
     all_matches = match_recs;        
     match_TF    = all_matches(:,2) + 0.0001*floor(all_matches(:,3)); % firm type by firmID identifier
     TF_list     = unique(match_TF);
     NumTF       = size(TF_list,1);    
     TF_matdat   = cell(NumTF,1);
     mat_lifecycle_TF = cell(NumTF,1);
     firstYrCount = zeros(NumTF,1);
 
%% Create matrix of match life cycles, one row per match

  % Pick off first-year matches for each type-firm 
    for TF_id = 1:NumTF
      ff = match_TF == ones(size(all_matches,1),1).*TF_list(TF_id);
      TF_matdat{TF_id} = all_matches(ff,:);
      firstYrCount(TF_id,1) = sum(TF_matdat{TF_id}(:,8)==1);
    end
  
  % Build life cycle for each new match
    orphan_matches = double.empty(0,9);
  parfor TF_id=1:NumTF
    % fprintf('\r firm ID = %0.0f\n',TF_id);  
      all_age = TF_matdat{TF_id};
    % all_age: (1) year, (2) type, (3) firm_ID, (4) sales, (5) shipments, 
    %          (6) boy Z,(7) eoy Z,(8) match age,(9) firm age 
      new_match_TF = all_age(:,8)<=1;
      max_age_TF = max(all_age(:,8));
      N_TF = sum(new_match_TF); % number of new matches for this firm
      mat_lifecycle = zeros(N_TF,5*max_age);
      mat_lifecycle(:,1:5) = [all_age(new_match_TF,1),all_age(new_match_TF,8),...
               all_age(new_match_TF,6:7),all_age(new_match_TF,4)]; 
    % load one-year-olds into first 5 cols of mat_lifecycle       
    % mat_lifecycle(:,1:5): [year, match age (<=1), boy Z, eoy Z, sales]
     
      for aa = 2:max_age_TF
          % Grab the aa-year-old matches  
           aa_yr_old = all_age(:,8) == aa;
           aa_cohort = all_age(aa_yr_old,:);

      stayInLoop = 1;      
          while stayInLoop > 0 && size(aa_cohort,1) > 0
              
           if firstYrCount(TF_id,1)==0 
           fprintf('\r type-firm ID = %0.0f, number of 1st year matches = %0.0f\n',[TF_id,firstYrCount(TF_id,1)]);
           stayInLoop = 0;               
           end
              
              ii =1; 
              while ii <= N_TF %looping over matches of age aa for firm-type TF_id
              fprintf('\r type-firm ID = %0.0f, match age = %0.0f, match number = %0.0f',[TF_id,aa,ii]);
              lag_yr   = mat_lifecycle(ii,5*(aa-2)+1);
              lag_age  = mat_lifecycle(ii,5*(aa-2)+2);
              lag_eoyZ = mat_lifecycle(ii,5*(aa-2)+4);

            % Find compatible matches to splice with last year's match ii
              flg = (lag_eoyZ>0).*(aa_cohort(:,6)==lag_eoyZ).*...
                    (aa_cohort(:,8)==lag_age + 1).*(aa_cohort(:,1)==lag_yr + 1);
                  
             % If compatible matches are identified, put first one in the 
             % relevant bloc of mat_lifecycle and remove it from aa_cohort
              if sum(flg)>0
                  mat_cont = find(flg==1,1); % first compatible match in aa_cohort
                  lb = 5*(aa-1)+1;
                  ub = 5*(aa-1)+5;
                  mat_lifecycle(ii,lb:ub) = [aa_cohort(mat_cont,1),aa_cohort(mat_cont,8),...
                       aa_cohort(mat_cont,6:7),aa_cohort(mat_cont,4)]; 
              
                % remove matched record from aa_cohort 
                  keepers = ones(size(aa_cohort,1),1);
                  keepers(mat_cont,1) = 0;
                  stayInLoop = sum(keepers)>0;
                  aa_cohort = aa_cohort(logical(keepers),:); 
              else
                  stayInLoop = 0;
              end
              ii = ii + 1;
              end
          end
      if size(aa_cohort,1) > 0 % document leftover matches
        orphans = size(aa_cohort,1) > 0;
        orphan_matches = [orphan_matches;aa_cohort];
        fileID5 = fopen('results/EEJKT_orphan_log.txt','a');
        fprintf(fileID5,'\r %0.0f unmatched record(s), %0.0f year-old firm, type %3.6f',[orphans,aa,TF_id]);
        fclose(fileID5);
      end             
      end     
      mat_lifecycle_TF{TF_id} = [TF_id*ones(size(mat_lifecycle,1),1), mat_lifecycle];
   end       
  
  % Stack match histories for firm-types, putting firm-type ID in col. 1
    agg_mat_lifecycle = double.empty(0,5*max_age+1);    
    for TF_id = 1:NumTF
        agg_mat_lifecycle = [agg_mat_lifecycle; mat_lifecycle_TF{TF_id}];
    end
    
%% Construct match survival table

% Find sales quartiles for new matches
     yr1         = all_matches(:,8)<=1; % pick matches in first year   
     new_matches = all_matches(yr1,:);  % matches in their 1st yr.    
     new_sales   = sort(new_matches(:,4),1); % sort by sales of new matches
     Nrows       = size(new_sales,1);
     cum_cnt     = cumsum(ones(Nrows,1)./Nrows);    
     mean_q      = zeros(4,1);
     match_cnt_q = zeros(4,1);
     size_cutoff = zeros(4,1);
     
     upb = 0.25;
     lowb = 0;
     for ii = 1:4
     % vector of row identifiers for types in sales quartile ii:
        ndx_q  = cum_cnt>lowb & cum_cnt<=upb; 
        lowb = upb;
        upb = upb + 0.25;
     % quartile-wide avg. sales per match:
        mean_q(ii)      = mean(new_sales(ndx_q,1));
        match_cnt_q(ii) = size(new_sales(ndx_q,1),1);
        size_cutoff(ii) = max(new_sales(ndx_q));  
     end
     
   % create dummies for initial size quartiles   
     agg_mat_lifecycle = sortrows(agg_mat_lifecycle,6); % sort by 1st yr sales
     Nquartile         = floor(size(agg_mat_lifecycle,1)/4);
     sizedum           = logical(kron(eye(4),ones(Nquartile,1)));
   % group leftover big firms in largest size category, if necessary
     resid             = size(agg_mat_lifecycle,1) - 4*Nquartile;
     sizedum           = [sizedum;logical(ones(resid,1)*[0 0 0 1])]; 

% Find match survival by initial size quartile and age 
     maxAge = 10;   
     survive    = zeros(maxAge,4);
     meanSale   = zeros(maxAge,4);
     matchCount = zeros(maxAge,4);
     for qq=1:4
          matchBySize = agg_mat_lifecycle(sizedum(:,qq),:);
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
firmCount      = zeros(maxAge,1); 
cohortSales    = zeros(maxAge,1); 
cohortAvgSales = zeros(maxAge,1); 
 
   for aa=1:maxAge
  
      matchByAge = all_matches(all_matches(:,9)==aa,:); % matches for aa-yr-old firms
      match_YTF  = 0.0001*matchByAge(:,1) + matchByAge(:,2)... 
                 + 0.0000001*matchByAge(:,3); % year-type-firmID identifier
      YTF_list      = unique(match_YTF);
      firmCount(aa) = size(YTF_list,1); % # year-type-firms w/ aa-yr old matches
      
      YTFtotSales   = zeros(firmCount(aa),1);
      for YTF_id = 1:size(YTF_list) % aggregate aa-year-old match sales by year-type-firmID
         selectYTF =  match_YTF == YTF_list(YTF_id);
         YTFtotSales(YTF_id) = sum(matchByAge(selectYTF,4),1);    
      end
    % sum and average all aa-yr-old firm sales over years and firm-types 
      cohortSales(aa)    = sum(YTFtotSales);
      cohortAvgSales(aa) = cohortSales(aa)/firmCount(aa); % avg. aa-yr-old firm sales    
    end    
      
brooks = [firmCount,cohortSales,cohortAvgSales];

brooks = brooks(1:BrooksYrs,:);
%% Construct degree distribution (used for graph by summary_table)

DegreeDistCount = degree_dist(match_recs,mm);

save('match_summary_out_8-3-23b');

end

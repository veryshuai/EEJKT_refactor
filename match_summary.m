function [DegreeDistCount,exit_by_age,brooks] =  match_summary(simMoms,mm)

% This function generates failure rates by initial sales and match age.
% It also constructs a summary (Brooks) table for the cohort maturation process.

%%  Preliminary data preparation
     succ_matches = simMoms.agg_mat_yr_sales;
     dud_matches  = simMoms.agg_dud_matches;
     dud_matches(:,8) = 1; % set all dud matches to 1 month old
     match_recs   = [succ_matches;dud_matches]; % stack successful and dud matches            
  %  match_recs: [period, type, firm_ID, sales, shipments, boy Z, eoy Z, match age, firm age]       

%      after_burn      = match_recs(:,8) <= match_recs(:,1) - (mm.burn*mm.pd_per_yr); 
%     match_recs      = match_recs(after_burn,:); % drop matches that started during burn-in period
%     max_age         = max(match_recs(:,8));

%     match_TF    = match_recs(:,2) + 0.0001*floor(match_recs(:,3)); % firm type by firmID identifier
%     TF_list     = unique(match_TF);
%     NumTF       = size(TF_list,1);    
%     TF_matdat   = cell(NumTF,1);
%  
% %% Create matrix of match life cycles, one row per match
% 
%   % collect matches into cells by firm type and firmID 
%     for TF_id = 1:NumTF
%       ff = match_TF == ones(size(match_TF,1),1).*TF_list(TF_id);
%       TF_matdat{TF_id} = match_recs(ff,:);
%     end
% 
% %     [agg_mat_lifecycle,agg_orphan_matches] = lifecycle(NumTF,TF_matdat,mm,max_age);
%       
%  %  Adjust multi-year life cycles if zero shipments in first year 
%     shift_lifecycle =  agg_mat_lifecycle;
%     dormant = logical((agg_mat_lifecycle(:,6)==0).*(agg_mat_lifecycle(:,7)>0)); 
%     shift_lifecycle(dormant,2:end-5) = agg_mat_lifecycle(dormant,7:end); 
% 
%     for ndx = 3:5:(5*ceil(max_age/mm.pd_per_yr)+1)
%           shift_lifecycle(:,ndx) = max(shift_lifecycle(:,ndx) - dormant, 0);
%     end
%      
%  % When next year has 0 sales, 0 zero out remaining match-years 
%    for ndx = 6:5:(5*ceil(max_age/mm.pd_per_yr))
%      lastYr = logical((shift_lifecycle(:,ndx)>0).*(shift_lifecycle(:,ndx+5)==0));
%      shift_lifecycle(lastYr,ndx+1:end) = 0;
%    end
%  
%  % Drop unpopulated rows. (Some matches died before generating shipments.)   
%    all_matches = shift_lifecycle(shift_lifecycle(:,6)>0,:);                            
%     
% %% Construct match survival table
% 
% % Find sales quartiles for new matches
%      new_matches = all_matches(:,1:6);       % matches in their 1st yr.    
%      new_sales   = sort(new_matches(:,6),1); % sort by sales of new matches
%      Nrows       = size(new_sales,1);
%      cum_cnt     = cumsum(ones(Nrows,1)./Nrows);    
%      mean_q      = zeros(4,1);
%      match_cnt_q = zeros(4,1);
%      size_cutoff = zeros(4,1);
%      
%      upb = 0.25;
%      lowb = 0;
%      for ii = 1:4
%      % vector of row identifiers for types in sales quartile ii:
%         ndx_q  = cum_cnt>lowb & cum_cnt<=upb; 
%         lowb = upb;
%         upb = upb + 0.25;
%      % quartile-wide avg. sales per match:
%         mean_q(ii)      = mean(new_sales(ndx_q,1));
%         match_cnt_q(ii) = size(new_sales(ndx_q,1),1);
%         size_cutoff(ii) = max(new_sales(ndx_q));  
%      end
%      
%    % create dummies for initial size quartiles   
%      new_matches = sortrows(new_matches,6); % sort by 1st yr sales
%      Nquartile   = floor(size(new_matches,1)/4);
%      sizedum     = logical(kron(eye(4),ones(Nquartile,1)));
%    % group leftover big firms in largest size category, if necessary
%      resid       = size(new_matches,1) - 4*Nquartile;
%      sizedum     = [sizedum;logical(ones(resid,1)*[0 0 0 1])]; 
% 
% % Find match survival by initial size quartile and age 
%      maxAge = 11; % now expressed in years  
%      survive    = zeros(maxAge,4);
%      meanSale   = zeros(maxAge,4);
%      matchCount = zeros(maxAge,4);
%      all_matches = sortrows(all_matches,6);
%      for qq=1:4
%           matchBySize = all_matches(sizedum(:,qq),:);
%           age = 1;
%           for aa = 6:5:maxAge*5
%               meanSale(age,qq)   = mean(matchBySize(matchBySize(:,aa)>0,aa));
%               matchCount(age,qq) = sum(matchBySize(matchBySize(:,aa)>0));
%               if age > 1
%                 survive(age,qq) = matchCount(age,qq)/matchCount(age-1,qq);
%               end
%               age = age + 1;
%           end          
%      end 
%      exit_by_age = 1 - survive;
%       
% %%  Construct brooks table
% BrooksYrs =10;
% maxAge   = 50;
% firmCount      = zeros(maxAge,1); 
% cohortSales    = zeros(maxAge,1); 
% cohortAvgSales = zeros(maxAge,1); 
% % restate match age in years
% match_recs(:,9) = floor(match_recs(:,9)/(mm.pd_per_yr+1e-6))+1;
% 
%     for aa=1:maxAge
% 
%     % identify and count aa-yr-old firms in export mkt., by year
%       matchByAge = match_recs(match_recs(:,9)==aa,:); % matches for aa-yr-old firms
%     % Note: need to keep track of year because same hotel room can be
%     % occupied by different firms in different years  
%       match_YTF  = 0.0001*matchByAge(:,1) + matchByAge(:,2)... 
%                  + 0.0000001*matchByAge(:,3); % year-type-firmID identifier (in future, suggest Cantor pairing function)                    
%       YTF_list      = unique(match_YTF);
%       firmCount(aa) = size(YTF_list,1); 
%       
%     % aggregate aa-yr-old firms' match sales by year-type-firmID
%       YTFtotSales   = zeros(firmCount(aa),1);
%       for YTF_id = 1:size(YTF_list) 
%          selectYTF =  match_YTF == YTF_list(YTF_id);
%          YTFtotSales(YTF_id) = sum(matchByAge(selectYTF,4),1);    
%       end
%       
%     % sum and average all aa-yr-old firm sales over all years and firm-types 
%       cohortSales(aa)    = sum(YTFtotSales); % sales of aa-yr-olds, all years
%       cohortAvgSales(aa) = cohortSales(aa)/firmCount(aa); % avg. aa-yr-old firm sales    
%       
%     end    
%       
% brooks = [firmCount,cohortSales,cohortAvgSales];
% 
% brooks = brooks(1:BrooksYrs,:);
% %% Construct degree distribution (used for graph by summary_table)
% DegreeDistCount = degree_dist(match_recs,mm);

%save results/val_dat;
full_save_path = "results/exch_shock_plots/" + mm.save_name + mm.boot_iter_num;
save(char(full_save_path));
% save('match_summary_out_8-14-23');

DegreeDistCount = [];
exit_by_age = [];
brooks = [];
end
